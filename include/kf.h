// Filename:      kalmanfilter.h
// Description:   the file declares Kf class
//                about INS/GNSS based loosely-coupled method.
// Author:        Long Xingyu
// date:          July 8, 2020

#ifndef _KALMANFILTER_H_
#define _KALMANFILTER_H_

#include <algorithm>
#include <iterator>
#include <memory>
#include <vector>

#include "gnss.h"
#include "ins.h"

namespace nstd {
namespace navigation {
// 松组合零偏补偿使用的零偏误差函数模型
enum SensorErrMod {
  kNoCompensation = 9,  // 无传感器误差（零偏）补偿
  kGaussMarkov,         // Gauss-Markov
  kRandomWalk,          // RandomWalk
  kAr1,                 // AR(1)
  kAr2,                 // AR(2)
  kAr3,                 // AR(3)
  kArma11,              // ARMA(1,1)
  kArma22,              // ARMA(2,2)
  kArma33,              // ARMA(3,3)
};

class Kf {
 public:
  Kf(const ImuErr &, const double *observeNoise,
     const SensorErrMod = kRandomWalk,
     const std::vector<double> mod_arg = {1, 1, 1, 1, 1, 1});
  ~Kf() = default;
  void Update(Ins &, Gnss &, ImuFile &);
  void SetGnssBreakTime(const std::vector<double> &);
#ifndef NDEBUG
  std::vector<double> x_total_;
#endif

 private:
  void Resize();  // 重设与dim相关的容器成员大小
  void Iniq(const ImuErr &);
  void Inir(const double *noise);
  void IniSensorModArg(const std::vector<double> &arg);
  void Updateh(const Ins &);
  void Updateg(const Ins &);
  void Updateqk(double ts);
  void Updatestm(const Ins &);
  bool SingleUpdate(Ins &, const GnssDat &, const Pva &pre);
  bool GnssBreak(double gnss_time);
  void Feedback(Ins &inertial);
  const SensorErrMod mod_ = kRandomWalk;  // sensor compensition mode
  const size_t dim_ = 15;                 // dim of state
  std::vector<double> mod_arg_;  // argument of sensor compensition mode
  // state_est_ = [pos,vel,att]+[gyro bias,acce bias]
  // state vector estimation,[dim_*1]
  std::vector<double> x_;    // state vector, [dim_ * 1]
  std::vector<double> stm_;  // state transition matrix, [dim_* dim_]
  std::vector<double> g_;    // design matrix in system equation, [dim_ * 12]
  double q_[12 * 12]{0};     // spectral density matrix
  std::vector<double> qk_;   // driving noise matrix, [dim_*dim_]
  std::vector<double> k_;    // [dim_ * 6]
  std::vector<double> var_;  // eatimation weight [dim_ * dim_]
  std::vector<double> h_;    // design matrix in observation equation,[6 * dim_]
  double v_[6 * 1]{0};       // 新息
  double r_[6 * 6]{0};       // observation variance
  std::vector<double> gnss_break_time_;  // GNSS中断的时候

};  // class Kf
void KfTimeUpdate(const double *stm, const double *qk, const size_t &dim,
                  double *state, double *var);
void KfObUpdate(const double *obval, const double *r, const size_t &row_obval,
                const double *h, const size_t &row_h, const size_t &rol_h,
                double *state, double *var, const size_t &row_state);
}  // namespace navigation
}  // namespace nstd
#endif  // _KALMANFILTER_H_