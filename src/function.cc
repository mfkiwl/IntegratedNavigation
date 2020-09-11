#include "function.h"
namespace nstd {
namespace function {
using std::complex;
using std::vector;

FileType GetFileType(const std::string file_path) {
  std::string::size_type pos = file_path.rfind(".");
  std::string file_xetension_name = file_path.substr(pos + 1, 3);
  if (file_xetension_name == "bin") {
    return kBinary;
  } else {
    return kText;
  }
}
void SplitString(const std::string &s, std::vector<std::string> &v,
                 const std::string &c) {
  std::string::size_type pos1, pos2;
  pos2 = s.find(c);
  pos1 = 0;
  while (std::string::npos != pos2) {
    if (pos2 != pos1) {
      v.push_back(s.substr(pos1, pos2 - pos1));
    }
    pos1 = pos2 + c.size();
    pos2 = s.find(c, pos1);
  }
  if (pos1 != s.length()) v.push_back(s.substr(pos1));
}
Fraction DToFra(double x, const int lim) {
  std::string s = DToS(x);
  int sym = 1;
  size_t len = s.size();
  if (s[0] == '-') {
    s = s.substr(1, len - 1);
    sym = -1;
    --len;
  }
  std::string::size_type pos = s.find(".");
  std::string snum1 = s.substr(0, pos);
  std::string snum2 = s.substr(pos + 1, lim);
  size_t len_snum2 = snum2.size();
  long long den = pow(10, len_snum2);
  long long num = std::stol(snum1) * den + std::stol(snum2);
  long long gcd = GetGcd(num, den);
  num /= gcd;
  den /= gcd;
  return Fraction(sym * num, den);
}
void sort(double *y, const size_t &ny, const int index) {
  if (index > 0) {
    for (size_t i = 0; i < ny - 1; ++i) {
      for (size_t j = 1; j < ny - i; ++j) {
        if (y[j] < y[j - 1]) {
          double tmp = y[j];
          y[j] = y[j - 1];
          y[j - 1] = tmp;
        }
      }
    }
  }  // ascending order
  else {
    for (size_t i = 0; i < ny - 1; ++i) {
      for (size_t j = 1; j < ny - i; ++j) {
        if (y[j] > y[j - 1]) {
          double tmp = y[j];
          y[j] = y[j - 1];
          y[j - 1] = tmp;
        }
      }
    }
  }  // descending order
}
size_t MaxScript(const double *x, const size_t &n) {
  size_t max_script = 0;
  for (size_t i = 1; i < n; ++i) {
    if (x[i] > x[max_script]) max_script = i;
  }
  return max_script;
}
size_t MinScript(const double *x, const size_t &n) {
  size_t min_script = 0;
  for (size_t i = 1; i < n; ++i) {
    if (x[i] < x[min_script]) min_script = i;
  }
  return min_script;
}
double Min(const double *x, const size_t &nx) {
  double min = x[0];
  for (size_t i = 1; i < nx; ++i) {
    if (x[i] < min) min = x[i];
  }
  return min;
}
std::vector<double> Min(const double *x, const size_t &nx, const size_t &nmin) {
  double *x_cp = new double[nx]{0};
  memcpy(x_cp, x, nx * sizeof(double));
  sort(x_cp, nx, 1);
  std::vector<double> x_min;
  for (size_t i = 0; i < nmin; ++i) {
    x_min.push_back(x_cp[i]);
  }
  delete[] x_cp;
  return x_min;
}
}  // namespace function
}  // namespace nstd