// Filename       kf.cc
// Description    his file implements kf.h
// Author         Long Xingyu
// Date           July 8, 2020

#include "kf.h"

namespace nstd {
namespace navigation {
namespace {
void Inter(const nstd::navigation::Pva pre, const nstd::navigation::Pva aft,
           const double time, double* inter) {
  // the ratio of interpolation
  double ratio = (time - pre.time) / (aft.time - pre.time);
  // pos_dif = aft.pos - pre.pos
  double pos_dif[3]{0};
  cblas_dcopy(3, aft.pos, 1, pos_dif, 1);
  cblas_daxpy(3, -1, pre.pos, 1, pos_dif, 1);
  // pos_inter = pos_before + ratio * pos_dif;
  cblas_dcopy(3, pre.pos, 1, inter, 1);
  cblas_daxpy(3, ratio, pos_dif, 1, inter, 1);
  // vel_dif = aft.vel - pre.vel
  double vel_dif[3]{0};
  cblas_dcopy(3, aft.vel, 1, vel_dif, 1);
  cblas_daxpy(3, -1, pre.vel, 1, vel_dif, 1);
  // vel_inter = vel_before + ratio * vel_dif
  cblas_dcopy(3, pre.vel, 1, &inter[3], 1);
  cblas_daxpy(3, ratio, vel_dif, 1, &inter[3], 1);
}
}  // namespace
using std::shared_ptr;
using std::vector;

Kf::Kf(const ImuErr& err, const double* observeNoise, const SensorErrMod mod,
       const vector<double> mod_arg)
    : mod_(mod),
      dim_(mod == kNoCompensation
               ? 9
               : (mod == kGaussMarkov) || (mod == kRandomWalk) ||
                         (mod == kArma11) || (mod == kAr1)
                     ? 15
                     : mod == (kAr2 | kArma22) ? 21 : 27) {
  Resize();
  Iniq(err);
  Inir(observeNoise);
  IniSensorModArg(mod_arg);
}
void Kf::Resize() {
  stm_.resize(dim_ * dim_);
  x_.resize(dim_ * 1);
  g_.resize(dim_ * 12);
  qk_.resize(dim_ * dim_);
  h_.resize(6 * dim_);
  k_.resize(dim_ * 6);
  var_.resize(dim_ * dim_);
  switch (mod_) {
    case kGaussMarkov:
    case kRandomWalk:
    case kAr1:
      mod_arg_.resize(6);
      break;
    case kAr2:
    case kArma11:
      mod_arg_.resize(12);
      break;
    case kAr3:
      mod_arg_.resize(18);
      break;
    case kArma22:
      mod_arg_.resize(24);
      break;
    case kArma33:
      mod_arg_.resize(36);
      break;
    default:
      mod_arg_.resize(0);
      break;
  }
}
void Kf::Iniq(const ImuErr& err) {
  for (int i = 0; i < 3; i++) {
    q_[i * 12 + i] = pow(err.vrw[i], 2);
  }
  for (int i = 3; i < 6; i++) {
    q_[i * 12 + i] = pow(err.arw[i - 3], 2);
  }
  for (int i = 6; i < 9; i++) {
    q_[i * 12 + i] = pow(err.gyro_bais[i - 6], 2);
  }
  for (int i = 9; i < 12; i++) {
    q_[i * 12 + i] = pow(err.acce_bais[i - 9], 2);
  }
}
void Kf::Inir(const double* observeNoise) {
  for (size_t i = 0; i < 6; i++) {
    r_[i + i * 6] = pow(observeNoise[i], 2);
  }
}
void Kf::IniSensorModArg(const vector<double>& arg) {
  if (mod_arg_.size() != arg.size())
    throw std::runtime_error("输入的参数数量与模型所需数量不匹配");
  mod_arg_ = arg;
}
void Kf::Updateh(const Ins& iner) {
  for (int i = 0; i < 6; i++) {
    h_[i + i * 6] = 1;
  }
  h_[0 + 0 * 6] = iner.rmh_;
  h_[1 + 1 * 6] = iner.rnh_ * iner.cos_lat_;
}
void Kf::Updateg(const Ins& iner) {
  cblas_dscal(dim_ * 12, 0, g_.data(), 1);  // TODO: 考虑删除
  for (int i = 3; i < 6; i++) {
    g_[i + dim_ * 0] = iner.cnb_[(i - 3) + 3 * 0];
    g_[i + dim_ * 1] = iner.cnb_[(i - 3) + 3 * 1];
    g_[i + dim_ * 2] = iner.cnb_[(i - 3) + 3 * 2];
  }
  for (int i = 6; i < 9; i++) {
    g_[i + dim_ * 3] = -iner.cnb_[(i - 6) + 3 * 0];
    g_[i + dim_ * 4] = -iner.cnb_[(i - 6) + 3 * 1];
    g_[i + dim_ * 5] = -iner.cnb_[(i - 6) + 3 * 2];
  }
  // bias
  for (int i = 9; i < 15; i++) {
    g_[i + dim_ * (i - 3)] = 1;
  }
}
void Kf::Update(Ins& iner, Gnss& sat, ImuFile& imu) {
  shared_ptr<vector<ImuDat>> pimudat = imu.Data();
  shared_ptr<vector<GnssDat>> pgnssdat = sat.Data();
  shared_ptr<vector<Pva>> presult = imu.Result();
  // 根据imu数据的时间范围，确定GNSS数据在IMU时间范围内的区间
  auto itr_gnssright = std::find_if(
      pgnssdat->begin(), pgnssdat->end(),
      [pimudat](GnssDat i) { return i.time >= pimudat->front().time; });
  auto itr_gnssend = std::find_if(pgnssdat->rbegin(), pgnssdat->rend(),
                                  [pimudat](GnssDat i) {
                                    return i.time <= pimudat->back().time;
                                  })
                         .base();
  // 率先更新一个imu历元，便于接下来的时间内插
  iner.SingleUpdate(pimudat->front());
  presult->push_back(iner.nav_);
  for (auto itr = pimudat->begin() + 1; itr != pimudat->end(); ++itr) {
    iner.SingleUpdate(*itr);
    if (SingleUpdate(iner, *itr_gnssright, presult->back()) &&
        itr_gnssright != itr_gnssend - 1) {
#ifndef NDEBUG
      // 保存组合导航观测更新的状态量
      x_total_.push_back(itr_gnssright->time);
      std::copy(x_.begin(), x_.begin() + 9, std::back_inserter(x_total_));
#endif                               // 调试
      iner.nav_.time = itr_gnssright->time;  // 重设当前更新的导航时间
      ++itr_gnssright;
    }
    else if (GnssBreak(itr_gnssright->time)) {
      ++itr_gnssright;
    }  // 防止人为屏蔽部分GNSS历元后，未执行观测更新进而GNSS历元永远停止在该历元，导致接下来的所有IMU历元不进行观测更新
    Feedback(iner);
    presult->push_back(iner.nav_);
  }
}
bool Kf::SingleUpdate(Ins& iner, const GnssDat& satdat, const Pva& previous) {
  Updateg(iner);
  Updateqk(iner.ts_);
  Updateh(iner);
  Updatestm(iner);
  KfTimeUpdate(stm_.data(), qk_.data(), dim_, x_.data(), var_.data());
  if (!GnssBreak(satdat.time) && previous.time <= satdat.time &&
      iner.nav_.time > satdat.time) {
    double ob[6]{0};                              // ob = ins - gnss
    Inter(previous, iner.nav_, satdat.time, ob);  //线性内插
    cblas_daxpy(3, -1, satdat.pos, 1, ob, 1);
    cblas_daxpy(3, -1, satdat.vel, 1, &ob[3], 1);
    ob[0] *= iner.rmh_;
    ob[1] *= iner.rnh_ * iner.cos_lat_;
    KfObUpdate(ob, r_, 6, h_.data(), 6, dim_, x_.data(), var_.data(), dim_);
    return true;
  }
  return false;
}
void Kf::Updateqk(double ts) {
  cblas_dscal(dim_ * dim_, 0, qk_.data(), 1);  // set 0
  // tmp1 = stm_ * G, [dim_ * 12]
  vector<double> tmp1(dim_ * 12);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, dim_, 12, dim_, 1,
              stm_.data(), dim_, g_.data(), dim_, 0, tmp1.data(), dim_);
  // tmp2 = tmp1 * Q, [dim_ * 12]
  vector<double> tmp2(dim_ * 12);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, dim_, 12, 12, 1,
              tmp1.data(), dim_, q_, 12, 0, tmp2.data(), dim_);
  // tmp3 = tmp2 * G', [dim_ * dim_]
  vector<double> tmp3(dim_ * dim_);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, dim_, dim_, 12, 1,
              tmp2.data(), dim_, g_.data(), dim_, 0, tmp3.data(), dim_);
  // Qk = tmp3 * stm_' * t
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, dim_, dim_, dim_, 1,
              tmp3.data(), dim_, stm_.data(), dim_, 0, qk_.data(), dim_);
  cblas_dscal(dim_ * dim_, ts, qk_.data(), 1);
  switch (mod_) {
    case kGaussMarkov:
      qk_[9 + dim_ * 9] = 2 * 1 / mod_arg_[0] * pow(q_[9 * 12 + 9], 2);
      qk_[10 + dim_ * 10] = 2 * 1 / mod_arg_[1] * pow(q_[10 * 12 + 10], 2);
      qk_[11 + dim_ * 11] = 2 * 1 / mod_arg_[2] * pow(q_[11 * 12 + 11], 2);
      qk_[12 + dim_ * 12] = 2 * 1 / mod_arg_[3] * pow(q_[6 * 12 + 6], 2);
      qk_[13 + dim_ * 13] = 2 * 1 / mod_arg_[4] * pow(q_[7 * 12 + 7], 2);
      qk_[14 + dim_ * 14] = 2 * 1 / mod_arg_[5] * pow(q_[8 * 12 + 8], 2);
      break;
    case kArma11:
      qk_[9 + dim_ * 9] = (1 + pow(mod_arg_[0 * 2 + 1], 2)) * qk_[9 + dim_ * 9];
      qk_[10 + dim_ * 10] =
          (1 + pow(mod_arg_[1 * 2 + 1], 2)) * qk_[10 + dim_ * 10];
      qk_[11 + dim_ * 11] =
          (1 + pow(mod_arg_[2 * 2 + 1], 2)) * qk_[11 + dim_ * 11];
      qk_[12 + dim_ * 12] =
          (1 + pow(mod_arg_[3 * 2 + 1], 2)) * qk_[12 + dim_ * 12];
      qk_[13 + dim_ * 13] =
          (1 + pow(mod_arg_[4 * 2 + 1], 2)) * qk_[13 + dim_ * 13];
      qk_[14 + dim_ * 14] =
          (1 + pow(mod_arg_[5 * 2 + 1], 2)) * qk_[14 + dim_ * 14];
      break;
    case kArma22:
      qk_[9 + dim_ * 9] =
          (1 + pow(mod_arg_[2 + 0 * 4], 2) + pow(mod_arg_[3 + 0 * 4], 2)) *
          qk_[9 + dim_ * 9];
      qk_[10 + dim_ * 10] =
          (1 + pow(mod_arg_[2 + 1 * 4], 2) + pow(mod_arg_[3 + 1 * 4], 2)) *
          qk_[10 + dim_ * 10];
      qk_[11 + dim_ * 11] =
          (1 + pow(mod_arg_[2 + 2 * 4], 2) + pow(mod_arg_[3 + 2 * 4], 2)) *
          qk_[11 + dim_ * 11];
      qk_[12 + dim_ * 12] =
          (1 + pow(mod_arg_[2 + 3 * 4], 2) + pow(mod_arg_[3 + 3 * 4], 2)) *
          qk_[12 + dim_ * 12];
      qk_[13 + dim_ * 13] =
          (1 + pow(mod_arg_[2 + 4 * 4], 2) + pow(mod_arg_[3 + 4 * 4], 2)) *
          qk_[13 + dim_ * 13];
      qk_[14 + dim_ * 14] =
          (1 + pow(mod_arg_[2 + 5 * 4], 2) + pow(mod_arg_[3 + 5 * 4], 2)) *
          qk_[14 + dim_ * 14];
      break;
    case kArma33:
      qk_[9 + 9 * dim_] =
          (1 + pow(mod_arg_[3 + 0 * 6], 2) + pow(mod_arg_[4 + 0 * 6], 2) +
           pow(mod_arg_[5 + 0 * 6], 2)) *
          qk_[9 + 9 * dim_];
      qk_[10 + 10 * dim_] =
          (1 + pow(mod_arg_[3 + 1 * 6], 2) + pow(mod_arg_[4 + 1 * 6], 2) +
           pow(mod_arg_[5 + 1 * 6], 2)) *
          qk_[10 * 10 * dim_];
      qk_[11 + 11 * dim_] =
          (1 + pow(mod_arg_[3 + 2 * 6], 2) + pow(mod_arg_[4 + 2 * 6], 2) +
           pow(mod_arg_[5 + 2 * 6], 2)) *
          qk_[11 * 11 * dim_];
      qk_[12 + 12 * dim_] =
          (1 + pow(mod_arg_[3 + 3 * 6], 2) + pow(mod_arg_[4 + 3 * 6], 2) +
           pow(mod_arg_[5 + 3 * 6], 2)) *
          qk_[12 * 12 * dim_];
      qk_[13 + 13 * dim_] =
          (1 + pow(mod_arg_[3 + 4 * 6], 2) + pow(mod_arg_[4 + 4 * 6], 2) +
           pow(mod_arg_[5 + 4 * 6], 2)) *
          qk_[13 * 13 * dim_];
      qk_[14 + 14 * dim_] =
          (1 + pow(mod_arg_[3 + 5 * 6], 2) + pow(mod_arg_[4 + 5 * 6], 2) +
           pow(mod_arg_[5 + 5 * 6], 2)) *
          qk_[14 * 14 * dim_];
      break;
    default:
      // nothing to do
      break;
  }
}
void Kf::Updatestm(const Ins& iner) {
  cblas_dscal(dim_ * dim_, 0, stm_.data(), 1);  // set 0
  double frr[3 * 3]{0};
  frr[1] = iner.nav_.vel[1] * iner.sin_lat_ / (iner.rnh_ * iner.cos_lat2_);
  frr[6] = -iner.nav_.vel[0] / pow(iner.rmh_, 2);
  frr[7] = -iner.nav_.vel[1] / (pow(iner.rnh_, 2) * iner.cos_lat_);
  double frv[3 * 3]{0};
  frv[0] = 1.0 / iner.rmh_;
  frv[4] = 1.0 / (iner.cos_lat_ * iner.rnh_);
  frv[8] = -1.0;
  double fvr[3 * 3]{0};
  fvr[0] = -2 * iner.nav_.vel[1] * WIE * iner.cos_lat_ -
           pow(iner.nav_.vel[1], 2) / (iner.rnh_ * iner.cos_lat2_);
  fvr[1] = 2 * WIE *
               (iner.nav_.vel[0] * iner.cos_lat_ -
                iner.nav_.vel[2] * iner.sin_lat_) +
           (iner.nav_.vel[1] * iner.nav_.vel[0]) / (iner.rnh_ * iner.cos_lat2_);
  fvr[2] = 2 * iner.nav_.vel[1] * WIE * iner.sin_lat_;
  fvr[6] = -(iner.nav_.vel[0] * iner.nav_.vel[2] / pow(iner.rmh_, 2)) +
           (pow(iner.nav_.vel[1], 2) * iner.tan_lat_ / pow(iner.rnh_, 2));
  fvr[7] =
      -(iner.nav_.vel[1] * iner.nav_.vel[2] / pow(iner.rnh_, 2)) -
      (iner.nav_.vel[0] * iner.nav_.vel[1] * iner.tan_lat_ / pow(iner.rnh_, 2));
  double rmn = sqrt(iner.rm_ * iner.rn_);
  double gama = GRAVITY * pow(rmn / (rmn + iner.nav_.pos[2]), 2);
  fvr[8] = pow(iner.nav_.vel[1] / iner.rnh_, 2) +
           pow(iner.nav_.vel[0] / iner.rmh_, 2) -
           2 * gama / (rmn + iner.nav_.pos[2]);
  double fvv[3 * 3]{0};
  fvv[0] = iner.nav_.vel[2] / iner.rmh_;
  fvv[1] =
      2 * WIE * iner.sin_lat_ + (iner.nav_.vel[1] * iner.tan_lat_ / iner.rnh_);
  fvv[2] = -2 * iner.nav_.vel[0] / iner.rmh_;
  fvv[3] = -2 * WIE * iner.sin_lat_ -
           (2 * iner.nav_.vel[1] * iner.tan_lat_ / iner.rnh_);
  fvv[4] = (iner.nav_.vel[2] + iner.nav_.vel[0] * iner.tan_lat_) / iner.rnh_;
  fvv[5] = -2 * WIE * iner.cos_lat_ - 2 * iner.nav_.vel[1] / iner.rnh_;
  fvv[6] = iner.nav_.vel[0] / iner.rmh_;
  fvv[7] = 2 * WIE * iner.cos_lat_ + iner.nav_.vel[1] / iner.rnh_;
  double fer[3 * 3]{0};
  fer[0] = -WIE * iner.sin_lat_;
  fer[2] =
      -WIE * iner.cos_lat_ - iner.nav_.vel[1] / (iner.rnh_ * iner.cos_lat2_);
  fer[6] = -iner.nav_.vel[1] / pow(iner.rnh_, 2);
  fer[7] = iner.nav_.vel[0] / pow(iner.rmh_, 2);
  fer[8] = iner.nav_.vel[1] * iner.tan_lat_ / pow(iner.rnh_, 2);
  double fev[3 * 3]{0};
  fev[1] = -1.0 / iner.rmh_;
  fev[3] = 1.0 / iner.rnh_;
  fev[5] = -iner.tan_lat_ / iner.rnh_;
  double fn_anti_mat[3 * 3]{0};  // the antiskew symmetric matrix of fn
  matrix::VtrToAntisymMat(iner.fn_, fn_anti_mat);
  double wnin_anti_mat[3 * 3]{0};  // the antiskew symmetric matrix of wnin
  matrix::VtrToAntisymMat(iner.wnin_, wnin_anti_mat);
  // row 1-9
  //   col 1-3
  for (size_t i = 0; i < 3; ++i) {
    stm_[0 + i * dim_] = frr[0 + 3 * i];
    stm_[1 + i * dim_] = frr[1 + 3 * i];
    stm_[2 + i * dim_] = frr[2 + 3 * i];
    stm_[3 + i * dim_] = fvr[0 + 3 * i];
    stm_[4 + i * dim_] = fvr[1 + 3 * i];
    stm_[5 + i * dim_] = fvr[2 + 3 * i];
    stm_[6 + i * dim_] = fer[0 + 3 * i];
    stm_[7 + i * dim_] = fer[1 + 3 * i];
    stm_[8 + i * dim_] = fer[2 + 3 * i];
  }
  //   col 4-6
  for (size_t i = 3; i < 6; ++i) {
    stm_[0 + i * dim_] = frv[0 + 3 * (i - 3)];
    stm_[1 + i * dim_] = frv[1 + 3 * (i - 3)];
    stm_[2 + i * dim_] = frv[2 + 3 * (i - 3)];
    stm_[3 + i * dim_] = fvv[0 + 3 * (i - 3)];
    stm_[4 + i * dim_] = fvv[1 + 3 * (i - 3)];
    stm_[5 + i * dim_] = fvv[2 + 3 * (i - 3)];
    stm_[6 + i * dim_] = fev[0 + 3 * (i - 3)];
    stm_[7 + i * dim_] = fev[1 + 3 * (i - 3)];
    stm_[8 + i * dim_] = fev[2 + 3 * (i - 3)];
  }
  //   col 7-9
  for (size_t i = 6; i < 9; ++i) {
    stm_[3 + i * dim_] = fn_anti_mat[0 + 3 * (i - 6)];
    stm_[4 + i * dim_] = fn_anti_mat[1 + 3 * (i - 6)];
    stm_[5 + i * dim_] = fn_anti_mat[2 + 3 * (i - 6)];
    stm_[6 + i * dim_] = -wnin_anti_mat[0 + 3 * (i - 6)];
    stm_[7 + i * dim_] = -wnin_anti_mat[1 + 3 * (i - 6)];
    stm_[8 + i * dim_] = -wnin_anti_mat[2 + 3 * (i - 6)];
  }
  //   if dim_ > 9, so state contains bias
  //   col 10-12
  if (dim_ > 9) {
    for (size_t i = 12; i < 15; ++i) {
      stm_[3 + i * dim_] = iner.cnb_[0 + 3 * (i - 12)];
      stm_[4 + i * dim_] = iner.cnb_[1 + 3 * (i - 12)];
      stm_[5 + i * dim_] = iner.cnb_[2 + 3 * (i - 12)];
    }
    for (size_t i = 9; i < 12; ++i) {
      stm_[6 + i * dim_] = -iner.cnb_[0 + 3 * (i - 9)];
      stm_[7 + i * dim_] = -iner.cnb_[1 + 3 * (i - 9)];
      stm_[8 + i * dim_] = -iner.cnb_[2 + 3 * (i - 9)];
    }
  }
  // discretization
  vector<double> I(dim_ * dim_, 0);  // unit matrix
  for (size_t i = 0; i < dim_; i++) {
    I[i * dim_ + i] = 1;
  }
  cblas_dscal(dim_ * dim_, iner.ts_, stm_.data(), 1);
  cblas_daxpy(9 * dim_, 1, I.data(), 1, stm_.data(), 1);
  switch (mod_) {
    case kGaussMarkov:
      stm_[9 + dim_ * 9] = 1 - 1 / mod_arg_[0] * iner.ts_;
      stm_[10 + dim_ * 10] = 1 - 1 / mod_arg_[1] * iner.ts_;
      stm_[11 + dim_ * 11] = 1 - 1 / mod_arg_[2] * iner.ts_;
      stm_[12 + dim_ * 12] = 1 - 1 / mod_arg_[3] * iner.ts_;
      stm_[13 + dim_ * 13] = 1 - 1 / mod_arg_[4] * iner.ts_;
      stm_[14 + dim_ * 14] = 1 - 1 / mod_arg_[5] * iner.ts_;
      break;
    case kRandomWalk:
      stm_[9 + dim_ * 9] = 1;
      stm_[10 + dim_ * 10] = 1;
      stm_[11 + dim_ * 11] = 1;
      stm_[12 + dim_ * 12] = 1;
      stm_[13 + dim_ * 13] = 1;
      stm_[14 + dim_ * 14] = 1;
      break;
    case kAr1:
      stm_[9 + 9 * dim_] = mod_arg_[0];
      stm_[10 + 10 * dim_] = mod_arg_[1];
      stm_[11 + 11 * dim_] = mod_arg_[2];
      stm_[12 + 12 * dim_] = mod_arg_[3];
      stm_[13 + 13 * dim_] = mod_arg_[4];
      stm_[14 + 14 * dim_] = mod_arg_[5];
      break;
    case kAr2:
      stm_[9 + 9 * dim_] = mod_arg_[0 + 2 * 0];
      stm_[9 + 15 * dim_] = mod_arg_[1 + 2 * 0];
      stm_[10 + 10 * dim_] = mod_arg_[0 + 2 * 1];
      stm_[10 + 16 * dim_] = mod_arg_[1 + 2 * 1];
      stm_[11 + 11 * dim_] = mod_arg_[0 + 2 * 2];
      stm_[11 + 17 * dim_] = mod_arg_[1 + 2 * 2];
      stm_[12 + 12 * dim_] = mod_arg_[0 + 2 * 3];
      stm_[12 + 17 * dim_] = mod_arg_[1 + 2 * 3];
      stm_[13 + 13 * dim_] = mod_arg_[0 + 2 * 4];
      stm_[13 + 19 * dim_] = mod_arg_[1 + 2 * 4];
      stm_[14 + 14 * dim_] = mod_arg_[0 + 2 * 5];
      stm_[14 + 20 * dim_] = mod_arg_[1 + 2 * 5];
      // t-1
      stm_[15 + 9 * dim_] = 1;
      stm_[16 + 10 * dim_] = 1;
      stm_[17 + 11 * dim_] = 1;
      stm_[18 + 12 * dim_] = 1;
      stm_[19 + 13 * dim_] = 1;
      stm_[20 + 14 * dim_] = 1;
      break;
    case kAr3:
      stm_[9 + 9 * dim_] = mod_arg_[0 + 3 * 0];
      stm_[9 + 15 * dim_] = mod_arg_[1 + 3 * 0];
      stm_[9 + 21 * dim_] = mod_arg_[2 + 3 * 0];
      stm_[10 + 10 * dim_] = mod_arg_[0 + 3 * 1];
      stm_[10 + 16 * dim_] = mod_arg_[1 + 3 * 1];
      stm_[10 + 22 * dim_] = mod_arg_[2 + 3 * 1];
      stm_[11 + 11 * dim_] = mod_arg_[0 + 3 * 2];
      stm_[11 + 17 * dim_] = mod_arg_[1 + 3 * 2];
      stm_[11 + 23 * dim_] = mod_arg_[2 + 3 * 2];
      stm_[12 + 12 * dim_] = mod_arg_[0 + 3 * 3];
      stm_[12 + 18 * dim_] = mod_arg_[1 + 3 * 3];
      stm_[12 + 24 * dim_] = mod_arg_[2 + 3 * 3];
      stm_[13 + 13 * dim_] = mod_arg_[0 + 3 * 4];
      stm_[13 + 19 * dim_] = mod_arg_[1 + 3 * 4];
      stm_[13 + 25 * dim_] = mod_arg_[2 + 3 * 4];
      stm_[14 + 14 * dim_] = mod_arg_[0 + 3 * 5];
      stm_[14 + 20 * dim_] = mod_arg_[1 + 3 * 5];
      stm_[14 + 26 * dim_] = mod_arg_[2 + 3 * 5];
      // t-1
      stm_[15 + 9 * dim_] = 1;
      stm_[16 + 10 * dim_] = 1;
      stm_[17 + 11 * dim_] = 1;
      stm_[18 + 12 * dim_] = 1;
      stm_[19 + 13 * dim_] = 1;
      stm_[20 + 14 * dim_] = 1;
      // t-2
      stm_[21 + 15 * dim_] = 1;
      stm_[22 + 16 * dim_] = 1;
      stm_[23 + 17 * dim_] = 1;
      stm_[24 + 18 * dim_] = 1;
      stm_[25 + 19 * dim_] = 1;
      stm_[26 + 20 * dim_] = 1;
      break;
    case kArma11:
      stm_[9 + 9 * dim_] = mod_arg_[0 + 2 * 0];
      stm_[10 + 10 * dim_] = mod_arg_[0 + 2 * 1];
      stm_[11 + 11 * dim_] = mod_arg_[0 + 2 * 2];
      stm_[12 + 12 * dim_] = mod_arg_[0 + 2 * 3];
      stm_[13 + 13 * dim_] = mod_arg_[0 + 2 * 4];
      stm_[14 + 14 * dim_] = mod_arg_[0 + 2 * 5];
      break;
    case kArma22:
      stm_[9 + 9 * dim_] = mod_arg_[0 + 4 * 0];
      stm_[9 + 15 * dim_] = mod_arg_[1 + 4 * 0];
      stm_[10 + 10 * dim_] = mod_arg_[0 + 4 * 1];
      stm_[10 + 16 * dim_] = mod_arg_[1 + 4 * 1];
      stm_[11 + 11 * dim_] = mod_arg_[0 + 4 * 2];
      stm_[11 + 17 * dim_] = mod_arg_[1 + 4 * 2];
      stm_[12 + 12 * dim_] = mod_arg_[0 + 4 * 3];
      stm_[12 + 17 * dim_] = mod_arg_[1 + 4 * 3];
      stm_[13 + 13 * dim_] = mod_arg_[0 + 4 * 4];
      stm_[13 + 19 * dim_] = mod_arg_[1 + 4 * 4];
      stm_[14 + 14 * dim_] = mod_arg_[0 + 4 * 5];
      stm_[14 + 20 * dim_] = mod_arg_[1 + 4 * 5];
      // t-1
      stm_[15 + 9 * dim_] = 1;
      stm_[16 + 10 * dim_] = 1;
      stm_[17 + 11 * dim_] = 1;
      stm_[18 + 12 * dim_] = 1;
      stm_[19 + 13 * dim_] = 1;
      stm_[20 + 14 * dim_] = 1;
      break;
    case kArma33:
      stm_[9 + 9 * dim_] = mod_arg_[0 + 6 * 0];
      stm_[9 + 15 * dim_] = mod_arg_[1 + 6 * 0];
      stm_[9 + 21 * dim_] = mod_arg_[2 + 6 * 0];
      stm_[10 + 10 * dim_] = mod_arg_[0 + 6 * 1];
      stm_[10 + 16 * dim_] = mod_arg_[1 + 6 * 1];
      stm_[10 + 22 * dim_] = mod_arg_[2 + 6 * 1];
      stm_[11 + 11 * dim_] = mod_arg_[0 + 6 * 2];
      stm_[11 + 17 * dim_] = mod_arg_[1 + 6 * 2];
      stm_[11 + 23 * dim_] = mod_arg_[2 + 6 * 2];
      stm_[12 + 12 * dim_] = mod_arg_[0 + 6 * 3];
      stm_[12 + 18 * dim_] = mod_arg_[1 + 6 * 3];
      stm_[12 + 24 * dim_] = mod_arg_[2 + 6 * 3];
      stm_[13 + 13 * dim_] = mod_arg_[0 + 6 * 4];
      stm_[13 + 19 * dim_] = mod_arg_[1 + 6 * 4];
      stm_[13 + 25 * dim_] = mod_arg_[2 + 6 * 4];
      stm_[14 + 14 * dim_] = mod_arg_[0 + 6 * 5];
      stm_[14 + 20 * dim_] = mod_arg_[1 + 6 * 5];
      stm_[14 + 26 * dim_] = mod_arg_[2 + 6 * 5];
      // t-1
      stm_[15 + 9 * dim_] = 1;
      stm_[16 + 10 * dim_] = 1;
      stm_[17 + 11 * dim_] = 1;
      stm_[18 + 12 * dim_] = 1;
      stm_[19 + 13 * dim_] = 1;
      stm_[20 + 14 * dim_] = 1;
      // t-2
      stm_[21 + 15 * dim_] = 1;
      stm_[22 + 16 * dim_] = 1;
      stm_[23 + 17 * dim_] = 1;
      stm_[24 + 18 * dim_] = 1;
      stm_[25 + 19 * dim_] = 1;
      stm_[26 + 20 * dim_] = 1;
      break;
    default:
      break;
  }
}
bool Kf::GnssBreak(double gnss_time) {
  for (size_t i = 0; i < gnss_break_time_.size() / 2; i++) {
    if (gnss_time >= gnss_break_time_[2 * i] &&
        gnss_time <= gnss_break_time_[2 * i + 1]) {
      return true;
    }
  }
  return false;
}
void Kf::SetGnssBreakTime(const vector<double>& break_time) {
  if (break_time.size() % 2 != 0) {
    throw std::runtime_error("ERROR! Failed to set GNSS break time!");
  }
  for (auto itr = break_time.begin() + 1; itr != break_time.end() - 1; itr++) {
    if (*itr < *(itr - 1) || *itr > *(itr + 1)) {
      throw std::runtime_error("ERROR! Failed to set GNSS break time!");
    }
  }
  gnss_break_time_ = break_time;
}
void Kf::Feedback(Ins& inertial) {
  // velocity and position compensation
  for (int i = 0; i < 3; i++) {
    inertial.nav_.pos[i] -= x_[i];
    inertial.nav_.vel[i] -= x_[i + 3];
  }
  // attitude compensation
  double sita_error[3]{x_[6], x_[7], x_[8]};
  double sita_error_Antisymmat[9]{0};
  matrix::VtrToAntisymMat(sita_error, sita_error_Antisymmat);
  double I[9]{0};
  for (size_t i = 0; i < 3; i++) {
    I[i + i * 3] = 1;
  }
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1, I, 3,
              sita_error_Antisymmat, 3, 1, I, 3);
  double cnb[9]{0};
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1, I, 3,
              inertial.cnb_, 3, 0, cnb, 3);
  cblas_dcopy(9, cnb, 1, inertial.cnb_, 1);
  DCMToAttitude(inertial.cnb_, inertial.nav_.att);
  DCMToQuaternion(inertial.cnb_, inertial.qnb_);

  // bias compensation
  for (int i = 0; i < 3; i++) {
    inertial.err_.gyro_bais[i] += x_[9 + i];
    inertial.err_.acce_bais[i] += x_[12 + i];
  }
  cblas_dscal(dim_, 0, x_.data(), 1);
}
void KfTimeUpdate(const double* stm, const double* qk, const size_t& dim,
                  double* state, double* var) {
  // state = stm * state
  double* tmp_state = new double[dim]{0};
  cblas_dgemv(CblasColMajor, CblasNoTrans, dim, dim, 1, stm, dim, state, 1, 0,
              tmp_state, 1);
  memcpy(state, tmp_state, dim * sizeof(double));
  delete[] tmp_state;
  // var = stm* var * stm'+ qk
  double* tmp = new double[dim * dim]{0};  // tmp = stm * var
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, 1, stm,
              dim, var, dim, 0, tmp, dim);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, dim, dim, dim, 1, tmp,
              dim, stm, dim, 0, var, dim);
  delete[] tmp;
  cblas_daxpy(dim * dim, 1, qk, 1, var, 1);
}
void KfObUpdate(const double* obval, const double* r, const size_t& row_obval,
                const double* h, const size_t& row_h, const size_t& col_h,
                double* state, double* var, const size_t& row_state) {
  // row_state == col_h
  // row_obval == row_h
  // k = var * h' * (h * var * h' + r)^-1
  double* tmp1 = new double[row_state * row_h]{0};  // tmp1 = var * h'
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, row_state, row_h,
              row_state, 1, var, row_state, h, row_h, 0, tmp1, row_state);
  double* tmp2 = new double[row_h * row_h]{0};  // tmp2 = h * tmp1
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, row_h, row_h,
              row_state, 1, h, row_h, tmp1, row_state, 0, tmp2, row_h);
  cblas_daxpy(row_h * row_h, 1, r, 1, tmp2, 1);  // tmp2 += r
  matrix::Inv(tmp2, row_h);
  double* k = new double[row_state * row_h]{0};  // k = tmp1*tmp2
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, row_state, row_h,
              row_h, 1, tmp1, row_state, tmp2, row_h, 0, k, row_state);
  // v = (obval - h * state)
  double* v = new double[row_obval]{0};
  cblas_dcopy(row_obval, obval, 1, v, 1);
  cblas_dgemv(CblasColMajor, CblasNoTrans, row_h, col_h, -1, h, row_h, state, 1,
              1, v, 1);
  // state = state + k * v
  cblas_dgemv(CblasColMajor, CblasNoTrans, row_state, row_h, 1, k, row_state, v,
              1, 1, state, 1);
  // var = (i - k * h) * var
  double* unit_mat = new double[row_state * row_state]{0};
  for (size_t i = 0; i < row_state; i++) {
    unit_mat[i + i * row_state] = 1;
  }
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, row_state, col_h,
              row_h, -1, k, row_state, h, row_h, 1, unit_mat, row_state);
  double* tmp_var = new double[row_state * row_state]{0};
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, row_state, row_state,
              row_state, 1, unit_mat, row_state, var, row_state, 0, tmp_var,
              row_state);
  cblas_dcopy(row_state * row_state, tmp_var, 1, var, 1);
  delete[] tmp_var;
  delete[] unit_mat;
  delete[] v;
  delete[] k;
  delete[] tmp2;
  delete[] tmp1;
}
}  // namespace navigation
}  // namespace nstd