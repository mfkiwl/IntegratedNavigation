#ifndef _FUNCTION_H_
#define _FUNCTION_H_

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace nstd {
namespace function {
const double kPi = acos(-1);
enum FileType { kText, kBinary };
struct Fraction {
  Fraction() = default;
  Fraction(long long x, long long y) : num_(x), den_(y) {}
  long long num_ = 0;
  long long den_ = 0;
};  // struct Fraction
void SplitString(const std::string &s, std::vector<std::string> &v,
                 const std::string &c = " ");
FileType GetFileType(const std::string file_path);
inline size_t GetFileSize(std::ifstream &stream) {
  stream.seekg(0, std::ios::end);
  size_t read_len = stream.tellg();
  stream.seekg(0, std::ios::beg);
  return read_len;
}
inline bool Juge2Power(unsigned long n) {
  if (n < 1) return false;
  int idx = n & n * (n - 1);
  return idx == 0;
}
inline size_t Supple2Power(const size_t &n) {
  if (Juge2Power(n)) return n;
  return pow(2, ceil(log2(n)));
}
// sort by index,
//   index > 0 ascending order(default)
//   index <= 0 descending order
void sort(double *y, const size_t &ny, const int index = 1);
// double convert to fraction
Fraction DToFra(double x, const int lim = 4);
// return the script of max number of matrix
size_t MaxScript(const double *x, const size_t &nx);
// return max of matrix
template <typename T>
T Max(const T *x, const size_t &nx) {
  T max = 0;
  for (size_t i = 0; i < nx; ++i) {
    if (x[i] > max) max = x[i];
  }
  return max;
}
template <typename T>
std::vector<T> Max(const T *x, const size_t &nx, const size_t &nmax) {
  std::vector<T> x_max;
  if (nx < nmax) {
    return x_max;
  }
  T *x_cp = new T[nx]{0};
  memcpy(x_cp, x, nx * sizeof(T));
  sort(x_cp, nx, -1);
  for (size_t i = 0; i < nmax; ++i) {
    x_max.push_back(x_cp[i]);
  }
  delete[] x_cp;
  return x_max;
}
// return min
double Min(const double *x, const size_t &nx);
std::vector<double> Min(const double *x, const size_t &nx, const size_t &nmin);
// return the script of min number of matrix
size_t MinScript(const double *, const size_t &);
// 最大公约数(欧几里得法)
template <typename T>
T GetGcd(const T &a, const T &b) {
  T m, n, r;
  m = a;
  n = b;
  r = m % n;
  // m和n的最大公约数（m大于n）等于n和r的最大公约数
  while (r != 0) {
    m = n;
    n = r;
    r = m % n;
  }
  return n;
}
//  least common multiple
template <typename T>
T GetLcm(const T &a, const T &b) {
  T c = a * b;
  T gcd = GetGcd(a, b);
  return c / gcd;
}
// 多个数的最大公约数
template <typename T>
T GetGcd(const T *a, const size_t &na) {
  if (na == 1) return a[0];
  T gcd = GetGcd(a[0], a[1]);
  for (size_t i = 2; i < na; i++) {
    gcd = GetGcd(gcd, a[i]);
  }
  return gcd;
}
// 分数的最大公约数
template <>
inline Fraction GetGcd(const Fraction &a, const Fraction &b) {
  long long lcm = GetLcm(a.den_, b.den_);
  long long gcd = GetGcd(a.num_, b.num_);
  // 化简
  long long num = gcd / GetGcd(gcd, lcm);
  long long den = lcm / GetGcd(gcd, lcm);
  return Fraction(num, den);
}
// double convert to string
inline std::string DToS(double x) {
  std::stringstream s;
  s << std::setprecision(20) << x;
  return s.str();
}
// mean
template <typename T>
double Mean(const T *y, const size_t n) {
  double u = 0;
  for (size_t i = 0; i < n; i++) {
    u += y[i];
  }
  return u / n;
}
// varience
template <typename T>
double Var(const T *y, const size_t n) {
  double u = Mean(y, n);
  double sum = 0;
  for (size_t i = 0; i < n; ++i) {
    sum += pow(y[i] - u, 2);
  }
  return sum / (n - 1);
}
// standard deviation
template <typename T>
double Std(const T *y, const size_t n) {
  double var = Var(y, n);
  return pow(var, 0.5);
}
}  // namespace function
}  // namespace nstd

#endif  // _FUNCTION_H_