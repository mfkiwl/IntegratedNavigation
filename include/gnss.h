#ifndef _GPS_H_
#define _GPS_H_

#include <memory>
#include <vector>

#include "arg.h"
#include "function.h"
#include "ins.h"

namespace nstd {
namespace navigation {

// Gnss data
struct GnssDat {
  double time = 0;
  double precision{0};
  double vel[3]{0};
  double pos[3]{0};
};

class Gnss {
 public:
  Gnss() : data_(new std::vector<GnssDat>) {}
  bool ReadFile(std::string GPSFileName);
  // input: epoch is symbol of subscript, 0~end
  // output: a copy of specified epoch
  GnssDat GetData(const size_t epoch) { return (*data_)[epoch]; }
  std::shared_ptr<std::vector<GnssDat>> Data() { return data_; }
  std::vector<GnssDat>::iterator begin() { return data_->begin(); }
  std::vector<GnssDat>::iterator end() { return data_->end(); }

 private:
  std::shared_ptr<std::vector<GnssDat>> data_;
};
}  // namespace navigation
}  // namespace nstd
#endif  //_GPS_H_
