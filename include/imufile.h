// Filename:      imufile.h
// Description:   the file declares class imufile about reading imu file and
//                saving imu data in memory.
// Author:        Long Xingyu
// date:          July 1,2020
#ifndef _IMUFILE_H_
#define _IMUFILE_H_

#include <openblas/cblas.h>
#include <algorithm>
#include <fstream>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "arg.h"
#include "function.h"

namespace nstd {
namespace navigation {
// navigation result per epoch
// tmptime: time
// p: position, [latitude, longitude, altitude] unit: deg, deg, m
// v: velocity, [north, east, down] unit: m/s
// a: attitude, [roll, pitch, heading] unit: deg
struct Pva {
  Pva() = default;
  Pva(const double *p, const double *v, const double *a) : Pva(0, p, v, a) {}
  Pva(const double tmptime, const double *p, const double *v, const double *a)
      : time(tmptime) {
    cblas_dcopy(3, p, 1, pos, 1);
    cblas_dcopy(3, v, 1, vel, 1);
    cblas_dcopy(3, a, 1, att, 1);
    cblas_dscal(2, DEG_TO_RAD, pos, 1);
    cblas_dscal(3, DEG_TO_RAD, att, 1);
  }
  double time = 0;
  double pos[3]{0};  // vehicle postion, [latitude,longitude,altitude]
  double vel[3]{0};  // vehicle velocity,[N,E,D]
  double att[3]{0};  // vehicle attitude, [roll, pitch, course]
};

// format of imu data
enum DataFormat {
  // unkown format
  kNullImuDataFormat,
  // front-right-down, increment, gyroscope unit:rad accelerometer unit:m/s
  kFrdIncre,
  // front-right-down, rate, gyroscope unit:rad/s accelerometer unit:m/s^2
  kFrdRate,
  // ront-right-down, rate, gyroscope unit:rad/s, accelerometer unit:g
  kFrdRateG,
  // right-front-up, increment, gyroscope unit:rad accelerometer unit:m/s
  kRfuIncre,
  // right-front-up, rate, gyroscope unit:rad/s accelerometer unit:m/s^2
  kRfuRate,
  // right-front-up, rate, gyroscope unit:rad/s, accelerometer unit:g
  kRfuRateG,
};

// IMU data
// 此数据结构不强制要求陀螺仪与加速度计的单位，可根据自己的需要改变
struct ImuDat {
  ImuDat() = default;
  ImuDat(const double *g, const double *a) : ImuDat(0, g, a) {}
  ImuDat(const double t, const double *g, const double *a) {
    time = t;
    cblas_dcopy(3, g, 1, gyro, 1);
    cblas_dcopy(3, a, 1, acce, 1);
  }
  // 由于后续使用了指针访问连续存储的数据，所以声明顺序不能变
  double time = 0;    // sampling time
  double gyro[3]{0};  // data of gyroscope
  double acce[3]{0};  // data of accelerometer
};

// Warning:
//   the format of imu file must be time-goryscope-acceleration, which has 7
//   columns! the coordinate of imu data in file can be right-front-up or
//   front-righr-down. the data can be increment or rate of change. coordinate
//   of imu data being read in will be converted to increment in
//   FRONT-RIGHT-DOWN.
// Description:
//   save the IMU data and navigation results
// Example:
//   ImuFile data(path);
class ImuFile {
 public:
  explicit ImuFile(std::string path);
  ~ImuFile() = default;
  bool Empty() { return data_->empty(); }
  bool Save(std::string path);
  void Clear() { nav_result_.reset(); }
  std::shared_ptr<std::vector<ImuDat>> Data() const { return data_; }
  std::shared_ptr<std::vector<Pva>> Result() { return nav_result_; }

 private:
  std::shared_ptr<std::vector<ImuDat>>
      data_;  // imu数据, 陀螺仪单位rad, 加速度计单位m
  std::shared_ptr<std::vector<Pva>> nav_result_;  // 本次导航的定位定速定姿结果
};
}  // namespace navigation
}  // namespace nstd

#endif  // _IMUFILE_H_