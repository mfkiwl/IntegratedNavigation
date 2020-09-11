// Filename:      imufile.cc
// Description:   This file implements the class imufile
// Author:        Long Xingyu
// date:          July 1,2020

#include <imufile.h>

namespace nstd {
namespace navigation {

using nstd::navigation::ImuDat;
using std::shared_ptr;
using std::string;
using std::vector;

namespace {
bool ReadBinaryFile(const string path, shared_ptr<vector<ImuDat>> *pdata) {
  std::ifstream file_stream(path, std::ios::binary);
  if (!file_stream.is_open()) return false;
  size_t file_size = function::GetFileSize(file_stream);
  if (file_size % sizeof(ImuDat) != 0) throw string("数据格式不正确");
  (*pdata)->resize(file_size / sizeof(ImuDat));
  for (vector<ImuDat>::size_type i = 0; !file_stream.eof(); ++i) {
    file_stream.read(reinterpret_cast<char *>(&(*(*pdata))[i]), sizeof(ImuDat));
  }
  file_stream.close();
  return true;
}
void LineToImuData(string *line, ImuDat *p) {
  std::istringstream istring(*line);
  auto ptr = &(p->time);
  int i = 0;  // 控制循环
  for (string word; istring >> word && i < 7; ++ptr, ++i) {
    *ptr = std::stod(word);
  }
}
bool ReadTextFile(const string path, shared_ptr<vector<ImuDat>> *pdata) {
  std::ifstream file_stream(path, std::ios::binary);
  if (!file_stream.is_open()) return false;
  ImuDat tmppva;
  for (string line, word; getline(file_stream, line);) {
    LineToImuData(&line, &tmppva);
    (*pdata)->push_back(tmppva);
  }
  file_stream.close();
  return true;
}
DataFormat CheckDataFormat(shared_ptr<vector<ImuDat>> pdata) {
  double acce_z = pdata->front().acce[2];
  // 计算近似采样时间
  double t = (pdata->back().time - pdata->front().time) / (pdata->size() - 1);
  // 通过加速度计Z轴值大小判断坐标系以及数据类型（变化率或者变化量）
  if (acce_z < (-8.0 * t) && acce_z > (-12.0 * t)) return kFrdIncre;
  if (acce_z > -12 && acce_z < -8) return kFrdRate;
  if (acce_z > -1.3 && acce_z < -0.97) return kFrdRateG;
  if (acce_z > (8.0 * t) && acce_z < (12.0 * t)) return kRfuIncre;
  if (acce_z > 8 && acce_z < 12) return kRfuRate;
  if (acce_z > 0.98 && acce_z < 1.2) return kRfuRateG;
  return kNullImuDataFormat;
}
void DataFormatConvert(DataFormat format_type,
                       shared_ptr<vector<ImuDat>> *pdata) {
  double t =
      ((*pdata)->back().time - (*pdata)->front().time) / ((*pdata)->size() - 1);
  for (auto itr = (*pdata)->begin(); itr != (*pdata)->end(); ++itr) {
    double tmp = 0;
    switch (format_type) {
      case kFrdIncre:
        // nothing to do
        break;
      case kFrdRate:
        cblas_dscal(6, t, itr->gyro, 1);
        break;
      case kFrdRateG:
        cblas_dscal(3, t, itr->gyro, 1);
        cblas_dscal(3, t * GRAVITY, itr->acce, 1);
        break;
      case kRfuIncre:
        tmp = itr->gyro[0];
        itr->gyro[0] = itr->gyro[1];
        itr->gyro[1] = tmp;
        itr->gyro[2] = -itr->gyro[2];
        tmp = itr->acce[0];
        itr->acce[0] = itr->acce[1];
        itr->acce[1] = tmp;
        itr->acce[2] = -itr->acce[2];
        break;
      case kRfuRate:
        tmp = itr->gyro[0];
        itr->gyro[0] = itr->gyro[1];
        itr->gyro[1] = tmp;
        itr->gyro[2] = -itr->gyro[2];
        tmp = itr->acce[0];
        itr->acce[0] = itr->acce[1];
        itr->acce[1] = tmp;
        itr->acce[2] = -itr->acce[2];
        cblas_dscal(6, t, itr->gyro, 1);
        break;
      case kRfuRateG:
        tmp = itr->gyro[0];
        itr->gyro[0] = itr->gyro[1];
        itr->gyro[1] = tmp;
        itr->gyro[2] = -itr->gyro[2];
        tmp = itr->acce[0];
        itr->acce[0] = itr->acce[1];
        itr->acce[1] = tmp;
        itr->acce[2] = -itr->acce[2];
        cblas_dscal(3, t, itr->gyro, 1);
        cblas_dscal(3, t * GRAVITY, itr->acce, 1);
        break;
      default:
        break;
    }
  }
}
}  // namespace
ImuFile::ImuFile(const string path)
    : data_(new vector<ImuDat>), nav_result_(new vector<Pva>) {
  auto file_type = function::GetFileType(path);
  // 将数据保存至data_中,未检查数据单位及坐标系
  bool idx = false;
  if (file_type == function::kText)
    idx = ReadTextFile(path, &data_);
  else if (file_type == function::kBinary)
    idx = ReadBinaryFile(path, &data_);
  else
    throw string("未知文件类型");
  // 检查并获得数据格式，根据数据格式将原数据转化为前右下坐标系的变化量形式
  if (!idx) return;
  auto data_format = CheckDataFormat(data_);
  DataFormatConvert(data_format, &data_);
}
bool ImuFile::Save(string path) {
  std::ofstream os(path);
  if (!os.is_open()) return false;
  for (auto itr = nav_result_->begin(); itr != nav_result_->end(); ++itr) {
    os.write(reinterpret_cast<char *>(&*itr), sizeof(Pva));
  }
  return true;
}
}  // namespace navigation
}  // namespace nstd