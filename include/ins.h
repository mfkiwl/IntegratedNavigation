// Filename:      ins.h
// Description:   the file declares class Ins
//                about Inertial navigation system.
//                body coordinate: front-right-under;
//                navigation coordinate: north-east-down
// Author:        Long Xingyu
// date:          July 1, 2020

#ifndef _INS_H_
#define _INS_H_

#include <cmath>
#include <cstring>
#include <memory>
#include <vector>

#include "arg.h"
#include "imufile.h"
#include "matrix.h"
#include "openblas/cblas.h"

namespace nstd {
namespace navigation {

// IMU error
// gyroscope_bais: bais of gyroscope, unit: deg/h
// accelerometer_bais: bais of accelerometer, unit: m/s
// a: angular random walk, TODO: 单位说明
// v: velocity random walk TODO: 单位说明
struct ImuErr {
  ImuErr() = default;
  ImuErr(const double *gyroscope_bais, const double *accelerometer_bais,
         const double *a, const double *v) {
    for (size_t i = 0; i < 3; i++) {
      gyro_bais[i] = gyroscope_bais[i] * DEG_TO_RAD / 3600;
      acce_bais[i] = accelerometer_bais[i];
      arw[i] = a[i];
      vrw[i] = v[i];
    }
  }
  double gyro_bais[3]{0};  // bais of gyroscope
  double acce_bais[3]{0};  // bais of accelerometer
  double arw[3]{0};        // angular random walk
  double vrw[3]{0};        // velocity random walk
};

// Inertial navigation system
// body coordinate: front-right-down
// navigation coordinate: north-east-down
// function:
//   ins update, alignment, compute Allan Variance
// example:
//   double pos[3]{39.9886389, 116.34608889, 54};
//   double vel[3]{0};
//   double att[3]{-0.492417, -0.244984, 174.149};
//   Pva initial(pos, vel, att);
//   double fs = 200;
//   Ins inertial(initial, fs, (ImuErr: 可选参数));
//   ImuFile data("readpath");
//   inertial.Updata(data);
//   data.Save("writepath");
//   inertial.Update()
//     parameter description:
//       pos:[lattitude, longitude, elevation]
//       vel:[velociy-north, velociy-east, velocity-down]
//       att:[roll, pitch, heading]
//       fs: sample frequency
//       readpath, writepath: the path of file
class Ins {
  friend class Kf;

 public:
  Ins(const Pva &nav, const double &fs) : Ins(nav, fs, ImuErr()) {}
  Ins(const Pva &nav, const double &fs, const ImuErr &err);
  ~Ins() = default;
  void Updata(ImuFile dat);

 private:
  // only one epoch is updated
  void SingleUpdate(const ImuDat &dat);
  // attitude update by two-sample
  void AttitudeUpdate();
  // velocity update by two-sample
  void VelocityUpdate();
  // position update
  void PositionUpdate();
  // 计算地球相关参数
  // p: 位置
  // v: 速度
  void EarthArgCompute(const double *p, const double *v);
  // 数据成员
  Pva nav_;           // 当前时刻的导航结果
  Pva oldnav_;        // 上个时刻的导航结果
  double ts_ = 0.01;  // 采样周期,unit:s
  double cnb_[9]{0};  // 载体坐标系到导航坐标系的旋转余弦矩阵(DCM)
  double qnb_[4]{0};  // 载体坐标系到导航坐标系的四元数
  ImuErr err_;        // 传感器误差
  double mpv_[9]{0};  // 位置变化率系数矩阵, mpv_* v为位置的变化角速度
  ImuDat olddat_;      // 上个时刻的IMU数据
  ImuDat dat_;         // 当前时刻的IMU数据
  double dvel_[3]{0};  // 两时刻速度的增量(vt- vt-1)
  double fn_[3]{0};    // 比力， 在组合导航算法中需要使用
  // 地球相关参数
  double wnen_[3]{
      0};  // 导航坐标系相对地心地固坐标系(wen)的角速度在导航系的表示
  double weie_[3]{
      0, 0,
      WIE};  // 地心地固坐标系相对惯性系的角速度(wie)在地心地固坐标系下表示
  double wnie_[3]{0};  // wie在导航系的表示
  double wnin_[3]{0};  // 导航坐标系相对惯性系的角速度(win)在导航系下表示
  double sin_lat_{0};   // sin(Latitude)
  double cos_lat_{0};   // sin^2(latitude)
  double tan_lat_{0};   // tan(latitude)
  double sin_lat2_{0};  // sin^2(Latitude)
  double cos_lat2_{0};  // cos^2(Latitude)
  double rm_{0};        // 子午圈曲率半径
  double rn_{0};        // 卯酉圈曲率半径
  double rmh_{0};       // Rm+h
  double rnh_{0};       // Rn+h
  double g_ = GRAVITY;  // 当前位置的重力(不同于平均重力)
  double gn_[3]{0};     // 重力投影到导航系
  double vcorg_[3]{0};  // 有害加速度的速度增量

};  // class Ins
// 旋转余弦矩阵转姿态角
// Cnb => [roll, pitch, head]
//   n:NED, b:front-right-down
void DCMToAttitude(const double *DCM, double *att);
// 旋转余弦矩阵转四元数, 导航坐标系: NED
void DCMToQuaternion(const double *DCM, double *quaternion);
// 粗对准, 须确保IMU数据前5分钟静止
void Alignment(const ImuFile &dat, const double *pos, double (*attitude)[3]);
void Allan(const ImuFile &dat, double fs, std::vector<double> *avar_curve);
}  // namespace navigation
}  // namespace nstd

#endif  // _INS_H_