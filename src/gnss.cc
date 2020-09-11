#include <cstring>

#include <fstream>
#include <iostream>

#include "function.h"
#include "gnss.h"

namespace nstd {
namespace navigation {
using std::string;
using std::vector;
// read Gnss file
bool Gnss::ReadFile(std::string file_path) {
  std::ifstream file_stream(file_path, std::ios::binary);
  if (!file_stream.is_open()) return false;
  string tmp_str;
  vector<string> tmp_str_devided;
  GnssDat tmp_gnss_vtr;
  while (getline(file_stream, tmp_str)) {
    nstd::function::SplitString(tmp_str, tmp_str_devided);
    tmp_gnss_vtr.time = stod(tmp_str_devided[0]);
    tmp_gnss_vtr.pos[0] = stod(tmp_str_devided[1]) * DEG_TO_RAD;
    tmp_gnss_vtr.pos[1] = stod(tmp_str_devided[2]) * DEG_TO_RAD;
    tmp_gnss_vtr.pos[2] = stod(tmp_str_devided[3]);
    tmp_gnss_vtr.vel[0] = stod(tmp_str_devided[4]);
    tmp_gnss_vtr.vel[1] = stod(tmp_str_devided[5]);
    tmp_gnss_vtr.vel[2] = stod(tmp_str_devided[6]);
    tmp_str_devided.clear();
    data_->push_back(tmp_gnss_vtr);
  }
  return true;
}
}  // namespace navigation
}  // namespace nstd