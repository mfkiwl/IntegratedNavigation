// Filename:      ins.cc
// Description:   This file implements ins.h
// Author:        Long Xingyu
// date:          July 1,2020

#include "ins.h"

namespace nstd {
namespace navigation {
using std::vector;
namespace {
void GeneralAllan(const vector<double> &imudat, double ts,
                  vector<double> *avar_curve) {
  size_t allan_num = 100;       // 计算多少个allan方差
  avar_curve->resize(7 * 100);  // reset the vector size
  // n次allan方差的计算需要n组区间长度
  // 为了保证在双对数图像中点数均匀，所以区间长度需要以10为幂的指数增长
  std::vector<size_t> len_per_group(allan_num);
  // 以10为幂的不超过数据长度一半的最大指数
  size_t row_imudat = imudat.size() / 6;
  double max = log10(row_imudat / 2);
  // 指数取exponent_mean_length为间隔的等差数列
  double exponent_mean_length = (max - 1) * 1.0 / allan_num;
  for (size_t i = 0; i < allan_num; ++i) {
    len_per_group[i] = pow(10, (1 + i * exponent_mean_length));
    (*avar_curve)[i] = len_per_group[i] * ts;
  }
  std::vector<double> imudat_csum(row_imudat * 6);
  matrix::Csum(imudat.data(), imudat_csum.data(), row_imudat, 6);
  // computing Allan Variance
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j < allan_num; j++) {
      for (size_t k = 0; k < row_imudat - 2 * len_per_group[j]; ++k) {
        size_t sub1 = k + 2 * len_per_group[j] + i * row_imudat;
        size_t sub2 = k + len_per_group[j] + i * row_imudat;
        size_t sub3 = k + i * row_imudat;
        (*avar_curve)[j + (i + 1) * allan_num] += pow(
            imudat_csum[sub1] - 2 * imudat_csum[sub2] + imudat_csum[sub3], 2);
      }  // compute allan variance
      (*avar_curve)[j + (i + 1) * allan_num] /=
          (2 * (row_imudat - 2 * len_per_group[j]) * pow((*avar_curve)[j], 2));
      (*avar_curve)[j + (i + 1) * allan_num] =
          sqrt((*avar_curve)[j + (i + 1) * allan_num]);
    }  // there are allan_num avar points
  }    // each axis is computed allan variance. There are six axes
}
// 向量转四元数
void VectorToQuaternion(const double *vector, double *quaternion) {
  double norm_square =
      pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2);
  double s = 0;
  if (norm_square < 1.0e-10) {
    quaternion[0] = 1 - norm_square * (1.0 / 8.0 - norm_square / 384.0);
    s = 0.5 - norm_square * (1.0 / 48.0 - norm_square / 3840);
  }
  else {
    double n = sqrt(norm_square);
    quaternion[0] = cos(n / 2.0);
    s = sin(n / 2.0) / n;
  }
  quaternion[1] = vector[0] * s;
  quaternion[2] = vector[1] * s;
  quaternion[3] = vector[2] * s;
  norm_square = pow(quaternion[0], 2) + pow(quaternion[1], 2) +
                pow(quaternion[2], 2) + pow(quaternion[3], 2);  // norm of Q3
  if ((norm_square > 1.000001) || (norm_square < 0.999999)) {
    cblas_dscal(4, 1.0 / sqrt(norm_square), quaternion, 1);
  }
}
// 四元数乘法
void QuaternionMulti(const double *Q1, const double *Q2, double *Q3) {
  Q3[0] = Q1[0] * Q2[0] - Q1[1] * Q2[1] - Q1[2] * Q2[2] - Q1[3] * Q2[3];
  Q3[1] = Q1[1] * Q2[0] + Q1[0] * Q2[1] - Q1[3] * Q2[2] + Q1[2] * Q2[3];
  Q3[2] = Q1[2] * Q2[0] + Q1[3] * Q2[1] + Q1[0] * Q2[2] - Q1[1] * Q2[3];
  Q3[3] = Q1[3] * Q2[0] - Q1[2] * Q2[1] + Q1[1] * Q2[2] + Q1[0] * Q2[3];
  /*quaternion normalization*/
  double norm_square = pow(Q3[0], 2) + pow(Q3[1], 2) + pow(Q3[2], 2) +
                       pow(Q3[3], 2); /*norm of Q3*/
  if ((norm_square > 1.000001) || (norm_square < 0.999999)) {
    cblas_dscal(4, 1.0 / sqrt(norm_square), Q3, 1);
  }
}
// 四元数转旋转余弦矩阵
void QuaternionToDCM(double *Quaternion, double *DCM) {
  DCM[0] = pow(Quaternion[0], 2) + pow(Quaternion[1], 2) -
           pow(Quaternion[2], 2) - pow(Quaternion[3], 2);
  DCM[1] = 2 * (Quaternion[1] * Quaternion[2] + Quaternion[0] * Quaternion[3]);
  DCM[2] = 2 * (Quaternion[1] * Quaternion[3] - Quaternion[0] * Quaternion[2]);
  DCM[3] = 2 * (Quaternion[1] * Quaternion[2] - Quaternion[0] * Quaternion[3]);
  DCM[4] = pow(Quaternion[0], 2) - pow(Quaternion[1], 2) +
           pow(Quaternion[2], 2) - pow(Quaternion[3], 2);
  DCM[5] = 2 * (Quaternion[2] * Quaternion[3] + Quaternion[0] * Quaternion[1]);
  DCM[6] = 2 * (Quaternion[1] * Quaternion[3] + Quaternion[0] * Quaternion[2]);
  DCM[7] = 2 * (Quaternion[2] * Quaternion[3] - Quaternion[0] * Quaternion[1]);
  DCM[8] = pow(Quaternion[0], 2) - pow(Quaternion[1], 2) -
           pow(Quaternion[2], 2) + pow(Quaternion[3], 2);
}
// 姿态角转旋转余弦矩阵
// void AttitudeToQuaternion(const double *att, double *quaternion) {
//   double c[3]{cos(att[0] * 0.5), cos(att[1] * 0.5), cos(att[2] * 0.5)};
//   double s[3]{sin(att[0] * 0.5), sin(att[1] * 0.5), sin(att[2] * 0.5)};
//   quaternion[0] = c[0] * c[1] * c[2] + s[0] * s[1] * c[2];
//   quaternion[1] = s[0] * c[1] * c[2] - c[0] * s[1] * s[2];
//   quaternion[2] = c[0] * s[1] * c[2] + s[0] * c[1] * s[2];
//   quaternion[3] = c[0] * c[1] * s[2] - s[0] * s[1] * c[2];
//   // normalization
//   double norm_square = pow(quaternion[0], 2) + pow(quaternion[1], 2) +
//                        pow(quaternion[2], 2) + pow(quaternion[3], 2);
//   if ((norm_square > 1.0001) || (norm_square < 0.9999)) {
//     cblas_dscal(4, 1.0 / sqrt(norm_square), quaternion, 1);
//   }
// }
// attitude convert to DCM, [roll, pitch, head]=>Cnb,
//   n:NED, b:front-right-down
void AttitudeToDCM(const double *att, double (*cnb)[9]) {
  double c[3]{cos(att[0]), cos(att[1]), cos(att[2])};
  double s[3]{sin(att[0]), sin(att[1]), sin(att[2])};
  (*cnb)[0] = c[1] * c[2];
  (*cnb)[1] = c[1] * s[2];
  (*cnb)[2] = -s[1];
  (*cnb)[3] = -c[0] * s[2] + s[0] * s[1] * c[2];
  (*cnb)[4] = c[0] * c[2] + s[0] * s[1] * s[2];
  (*cnb)[5] = s[0] * c[1];
  (*cnb)[6] = s[0] * s[2] + c[0] * s[1] * c[2];
  (*cnb)[7] = -s[0] * c[2] + c[0] * s[1] * s[2];
  (*cnb)[8] = c[0] * c[1];
}
}  // namespace
Ins::Ins(const Pva &nav, const double &fs, const ImuErr &err)
    : nav_(nav), ts_(1.0 / fs), err_(err) {
  AttitudeToDCM(nav_.att, &cnb_);
  DCMToQuaternion(cnb_, qnb_);
  EarthArgCompute(nav_.pos, nav_.vel);  // 初始化地球相关参数
}
void Ins::Updata(ImuFile dat) {
  auto pdat = dat.Data();    // 原始数据智能指针
  auto pnav = dat.Result();  // 导航结果智能指针
  pnav->reserve(pdat->size());  // 为导航结果开辟不小于原数据大小的内存
  for (auto itr = pdat->begin(); itr != pdat->end(); ++itr) {
    SingleUpdate(*itr);
    pnav->push_back(nav_);
  }
}
void Ins::SingleUpdate(const ImuDat &dat) {
  nav_.time = dat.time;
  memcpy(&olddat_, &dat_, sizeof(ImuDat));
  memcpy(&dat_, &dat, sizeof(ImuDat));
  memcpy(&oldnav_, &nav_, sizeof(Pva));
  // bias compensation
  cblas_daxpy(3, -ts_, err_.gyro_bais, 1, dat_.gyro, 1);
  cblas_daxpy(3, -ts_, err_.acce_bais, 1, dat_.acce, 1);
  EarthArgCompute(nav_.pos, nav_.vel);
  VelocityUpdate();
  PositionUpdate();
  AttitudeUpdate();
}
void Ins::VelocityUpdate() {
  // vel_ex = vt-1 + 0.5(vt-1 - vt-2) = vt-1 + 0.5*dvel_
  // dvel_: 在dvel_未完成更新前, 为vt-1与vt-2间的速度增量
  //   dvel_更新后为上一时刻至当前时刻的速度增量
  // 根据外推速度，完成外推位置
  // 根据外推速度与外推位置完成地球相关参数外推，主要wnin与有害加速度的速度增量
  double vel_ex[3]{0};  // velocity extrapolating
  double pos_ex[3]{0};  // position extrapolating
  cblas_dcopy(3, nav_.vel, 1, vel_ex, 1);
  cblas_daxpy(3, 0.5, dvel_, 1, vel_ex, 1);
  cblas_dcopy(3, nav_.pos, 1, pos_ex, 1);
  cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 0.5 * ts_, mpv_, 3, vel_ex, 1,
              1, pos_ex, 1);
  EarthArgCompute(pos_ex, vel_ex);
  // vt = vt-1 + vsf + vcorg
  // vt-1: 上时刻的速度
  // vsf： 加速度计输出经过改正、经过坐标系投影变换后再经过积分变换的结果
  //   vsf = (I-T/2wninX)DCM * (dat_.acce + conic + sulling)
  //     dat_.acce + conic + sulling: 经过圆锥误差和划桨误差改正后的加速度计输出
  // vcorg: 有害加速度的速度增量
  // 即: dvel = vsf + vcorg为上一时刻至当前时刻的速度增量
  // computing conic error
  double conic[3]{0};  // conic error
  matrix::VtrCross(dat_.gyro, dat_.acce, conic);
  cblas_dscal(3, 0.5, conic, 1);
  // computing sulling error
  double sulling[3]{0};  // sculling error
  double sulling_part1[3]{0}, sulling_part2[3]{0};
  matrix::VtrCross(olddat_.gyro, dat_.acce, sulling_part1);
  matrix::VtrCross(olddat_.acce, dat_.gyro, sulling_part2);
  cblas_dcopy(3, sulling_part1, 1, sulling, 1);
  cblas_daxpy(3, 1, sulling_part2, 1, sulling, 1);
  cblas_dscal(3, 1.0 / 12, sulling, 1);
  // 加速度计输出改正 dat_.acce + conic + sulling
  double vsf_part1[3]{0};
  cblas_dcopy(3, dat_.acce, 1, vsf_part1, 1);
  cblas_daxpy(3, 1, sulling, 1, vsf_part1, 1);
  cblas_daxpy(3, 1, conic, 1, vsf_part1, 1);
  // I = (I-T/2wninX)
  double I[9]{1, 0, 0, 0, 1, 0, 0, 0, 1};  // eye(3)
  double tmp[3]{0};
  cblas_daxpy(3, -0.5 * ts_, wnin_, 1, tmp, 1);
  double tmp_anti_mat[9]{0};  // 反斜对称矩阵
  matrix::VtrToAntisymMat(tmp, tmp_anti_mat);
  cblas_daxpy(9, 1, tmp_anti_mat, 1, I, 1);
  // part2 = I * DCM
  double vsf_part2[9]{0};
  double vsf[3]{0};
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1, I, 3, cnb_,
              3, 0, vsf_part2, 3);
  // vsf = part2 * part1
  cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1, vsf_part2, 3, vsf_part1, 1,
              0, vsf, 1);
  // dvel_ = dvel_sf + vcorg
  for (int i = 0; i < 3; i++) {
    dvel_[i] = vsf[i] + vcorg_[i];
  }
  cblas_daxpy(3, 1, dvel_, 1, nav_.vel, 1);  // vt = vt-1 + dvel
  // fn_： 比力， 在组合导航算法中需要使用
  cblas_dscal(3, 1.0 / ts_, vsf, 1);
  cblas_dcopy(3, vsf, 1, fn_, 1);
}
void Ins::PositionUpdate() {
  // pt = pt-1 + mpv_t-1/2 * (vt-1 + vt) * 0.5
  // mpv_t-1/2: 外推的位置变化率系数矩阵
  // mpv_ = [1 / rmh_, 0, 0; 0, 1 / (rnh_ * cos(latitude)), 0; 0, 0, -1]
  // mpv_仅与cos(latitude)， rmh_, rnh_相关
  // rmh_ = rm + h, rnh_ = rm + h
  // rm_, rn_, cos(latitude)仅与维度相关
  // 由于在速度更新中已经更新了外推维度相关的参数rm_, rn_, cos(latitude)
  // 且本处仅重新计算外推高度，所以rm_, rn, cos(latitude)未变化
  // 所以本处仅计算外推的rmh_, rnh_即可获得mpv_，
  // 无需为了计算mpv_再次调用EarthArgCompute(), 导致增加额外开销
  double vel_2[3]{0};  // (vt-1+vt)/2
  cblas_daxpy(3, 0.5, nav_.vel, 1, vel_2, 1);
  cblas_daxpy(3, 0.5, oldnav_.vel, 1, vel_2, 1);
  double h_ex = nav_.pos[2] - vel_2[2] * ts_ * 0.5;
  rmh_ = rm_ + h_ex;
  rnh_ = rn_ + h_ex;
  mpv_[0] = 1.0 / rmh_;
  mpv_[4] = 1.0 / (rnh_ * cos_lat_);
  cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, ts_, mpv_, 3, vel_2, 1, 1,
              nav_.pos, 1);
}
void Ins::AttitudeUpdate() {
  // 二子样计算
  double sita_old_cross_sita[3]{0};  // error
  matrix::VtrCross(olddat_.gyro, dat_.gyro, sita_old_cross_sita);
  cblas_dscal(3, 1.0 / 12.0, sita_old_cross_sita, 1);
  double sita_correct[3]{0};  // equivalent rotation vector
  cblas_dcopy(3, dat_.gyro, 1, sita_correct, 1);
  cblas_daxpy(3, 1, sita_old_cross_sita, 1, sita_correct, 1);
  double Qbb[4]{0};  // the quaternion of body coordinate from b(t-1) to b(t)
  VectorToQuaternion(sita_correct, Qbb);
  double vector_sita_win[3]{0};  // sita_navigation = win*T
  //从t-1至t时刻导航坐标系旋转的角度sita_navigation所表达的四元数为Qn(t-1)n(t),我们需要的是Qn(t)n(t-1),所以是转过角度的负数
  cblas_daxpy(3, -ts_, wnin_, 1, vector_sita_win, 1);
  // the quaternion of navigation coordinate from n(t) to n(t-1)
  double Qnn[4]{0};
  VectorToQuaternion(vector_sita_win, Qnn);
  double Q_nt_bt_1[4]{0};
  QuaternionMulti(Qnn, qnb_, Q_nt_bt_1);
  QuaternionMulti(Q_nt_bt_1, Qbb, qnb_);
  QuaternionToDCM(qnb_, cnb_);
  DCMToAttitude(cnb_, nav_.att);
}
void Ins::EarthArgCompute(const double *p, const double *v) {
  sin_lat_ = sin(p[0]);
  cos_lat_ = cos(p[0]);
  tan_lat_ = sin_lat_ / cos_lat_;
  sin_lat2_ = pow(sin_lat_, 2);
  cos_lat2_ = pow(cos_lat_, 2);
  rm_ = RE * (1 - ECCENTRICITY2) / pow((1 - ECCENTRICITY2 * sin_lat2_), 1.5);
  rn_ = RE / sqrt(1 - ECCENTRICITY2 * sin_lat2_);
  rmh_ = rm_ + p[2];
  rnh_ = rn_ + p[2];
  wnie_[0] = WIE * cos_lat_;
  wnie_[1] = 0;
  wnie_[2] = -WIE * sin_lat_;
  wnen_[0] = v[1] / rnh_;
  wnen_[1] = -v[0] / rmh_;
  wnen_[2] = -(v[1] * tan_lat_) / rnh_;
  g_ = GRAVITY * (1 + 0.0052790414 * sin_lat2_ +
                  (-0.000003087691089) * pow(sin_lat2_, 2)) +
       (0.000003087691089 + 0.000000004397731 * sin_lat2_) * p[2] +
       0.000000000000721 * pow(p[2], 2);
  gn_[2] = g_;
  double w2ieen[3]{0};
  for (int i = 0; i < 3; i++) {
    wnin_[i] = wnie_[i] + wnen_[i];
    w2ieen[i] = 2 * wnie_[i] + wnen_[i];
  }
  // vcorg_
  vcorg_[0] = w2ieen[2] * v[1] - w2ieen[1] * v[2];
  vcorg_[1] = w2ieen[0] * v[2] - w2ieen[2] * v[0];
  vcorg_[2] = w2ieen[1] * v[0] - w2ieen[0] * v[1] + gn_[2];
  cblas_dscal(3, ts_, vcorg_, 1);
  mpv_[0] = 1.0 / rmh_;
  mpv_[4] = 1.0 / (rnh_ * cos_lat_);
  mpv_[8] = -1.0;
}
void DCMToAttitude(const double *DCM, double *att) {
  att[1] = atan(-DCM[2] / sqrt(pow(DCM[5], 2) + pow(DCM[8], 2)));  // pitch
  if (abs(DCM[2]) < 0.999) {
    att[0] = atan2(DCM[5], DCM[8]);  // roll
    att[2] = atan2(DCM[1], DCM[0]);  // head
  }
  else if (DCM[2] <= -0.999) {
    att[0] = 0;
    att[2] = atan2((DCM[7] - DCM[3]), (DCM[6] + DCM[4]));
  }
  else {
    att[0] = 0;
    att[2] = atan2((DCM[7] + DCM[3]), (DCM[6] - DCM[4])) + PI;
  }  // DCM[2] >= 0.999
  if (att[2] < 0) {
    att[2] += 2 * PI;
  }
}
void DCMToQuaternion(const double *DCM, double *quaternion) {
  //  naser
  double DCM_trace = DCM[0] + DCM[4] + DCM[8];  // the trace of DCM
  double P1 = 1 + DCM_trace;
  double P2 = 1 + 2 * DCM[0] - DCM_trace;
  double P3 = 1 + 2 * DCM[4] - DCM_trace;
  double P4 = 1 + 2 * DCM[8] - DCM_trace;
  if ((P1 > P2) && (P1 > P3) && (P1 > P4)) {
    quaternion[0] = 0.5 * sqrt(P1);
    quaternion[1] = 0.25 * (DCM[5] - DCM[7]) / quaternion[0];
    quaternion[2] = 0.25 * (DCM[6] - DCM[2]) / quaternion[0];
    quaternion[3] = 0.25 * (DCM[1] - DCM[3]) / quaternion[0];
  }
  else if ((P2 > P1) && (P2 > P3) && (P2 > P4)) {
    quaternion[1] = 0.5 * sqrt(P2);
    quaternion[2] = 0.25 * (DCM[1] + DCM[3]) / quaternion[1];
    quaternion[3] = 0.25 * (DCM[6] + DCM[2]) / quaternion[1];
    quaternion[0] = 0.25 * (DCM[5] - DCM[7]) / quaternion[1];
  }
  else if ((P3 > P1) && (P3 > P2) && (P3 > P4)) {
    quaternion[2] = 0.5 * sqrt(P3);
    quaternion[3] = 0.25 * (DCM[5] + DCM[7]) / quaternion[2];
    quaternion[0] = 0.25 * (DCM[6] - DCM[2]) / quaternion[2];
    quaternion[1] = 0.25 * (DCM[1] + DCM[3]) / quaternion[2];
  }
  else if ((P4 > P1) && (P4 > P2) && (P4 > P3)) {
    quaternion[3] = 0.5 * sqrt(P4);
    quaternion[0] = 0.25 * (DCM[1] - DCM[3]) / quaternion[3];
    quaternion[1] = 0.25 * (DCM[6] + DCM[2]) / quaternion[3];
    quaternion[2] = 0.25 * (DCM[5] + DCM[7]) / quaternion[3];
  }
  // normalization
  double norm = sqrt(pow(quaternion[0], 2) + pow(quaternion[1], 2) +
                     pow(quaternion[2], 2) + pow(quaternion[3], 2));
  if (abs(norm - 1) > 0.000000001) {
    cblas_dscal(4, 1.0 / norm, quaternion, 1);
  }
}
void Alignment(const ImuFile &dat, const double *pos, double (*attitude)[3]) {
  double gyro_mean[3]{0};
  double acce_mean[3]{0};
  double ts = (dat.Data()->back().time - dat.Data()->front().time) * 1.0 /
              (dat.Data()->size() - 1);
  double fs = static_cast<size_t>(1.0 / ts);
  // mean
  for (auto itr = dat.Data()->begin(); itr != dat.Data()->begin() + 300 * fs;
       itr++) {
    cblas_daxpy(3, 1, itr->gyro, 1, gyro_mean, 1);  // 陀螺仪求和
    cblas_daxpy(3, 1, itr->acce, 1, acce_mean, 1);  // 加表求和
  }
  cblas_dscal(3, 1.0 / (300 * fs), gyro_mean, 1);
  cblas_dscal(3, 1.0 / (300 * fs), acce_mean, 1);
  // increment convert to intanous
  cblas_dscal(3, 1 / ts, gyro_mean, 1);
  cblas_dscal(3, 1 / ts, acce_mean, 1);

  double vtr_corss_imu[3];       // data outputed by IMU corss
  double vector_corss_earth[3];  // data outputed by earth corss
  matrix::VtrCross(acce_mean, gyro_mean, vtr_corss_imu);  // imu
  // the earth gravity
  double g =
      GRAVITY * (1 + 0.0052790414 * pow(sin(pos[0]), 2) +
                 (-0.000003087691089) * pow(sin(pos[0]), 4)) +
      (0.000003087691089 + 0.000000004397731 * pow(sin(pos[0]), 2)) * pos[2] +
      0.000000000000721 * pow(pos[2], 2);
  double gb[3]{0, 0, -g};  // 比力
  double wnie[3]{WIE * cos(pos[0]), 0, -WIE * sin(pos[0])};
  matrix::VtrCross(gb, wnie, vector_corss_earth);
  double rb[9] = {acce_mean[0], gyro_mean[0], vtr_corss_imu[0],
                  acce_mean[1], gyro_mean[1], vtr_corss_imu[1],
                  acce_mean[2], gyro_mean[2], vtr_corss_imu[2]};  // col major
  double rn[9] = {gb[0], wnie[0], vector_corss_earth[0],
                  gb[1], wnie[1], vector_corss_earth[1],
                  gb[2], wnie[2], vector_corss_earth[2]};  // col major
  matrix::Inv(rn, 3);
  double DCM[9]{0};
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1, rn, 3, rb,
              3, 0, DCM, 3);
  DCMToAttitude(DCM, *attitude);
  cblas_dscal(3, RAD_TO_DEG, *attitude, 1);
}
void Allan(const ImuFile &dat, double fs, vector<double> *avar_curve) {
  // copy vector<ImuDat>  to vector<double>([line][6]) without time,
  // and change the unit of gyroscope from rad to deg/h
  // matrix layout by line priority
  auto pdat = dat.Data();
  size_t col = pdat->size();
  vector<double> imudat(col * 6);
  for (size_t i = 0; i < col; ++i) {
    imudat[i + 0 * col] = (*pdat)[i].gyro[0] * RAD_TO_DEG * 3600;
    imudat[i + 1 * col] = (*pdat)[i].gyro[1] * RAD_TO_DEG * 3600;
    imudat[i + 2 * col] = (*pdat)[i].gyro[2] * RAD_TO_DEG * 3600;
    imudat[i + 3 * col] = (*pdat)[i].acce[0];
    imudat[i + 4 * col] = (*pdat)[i].acce[1];
    imudat[i + 5 * col] = (*pdat)[i].acce[2];
  }
  GeneralAllan(imudat, 1.0 / fs, avar_curve);
  return;
}
}  // namespace navigation
}  // namespace nstd