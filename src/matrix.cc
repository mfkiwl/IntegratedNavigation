#include "matrix.h"

namespace nstd {
namespace matrix {
void VtrToAntisymMat(const double *vtr, double *mat) {
  mat[0] = 0;
  mat[1] = vtr[2];
  mat[2] = -vtr[1];
  mat[3] = -vtr[2];
  mat[4] = 0;
  mat[5] = vtr[0];
  mat[6] = vtr[1];
  mat[7] = -vtr[0];
  mat[8] = 0;
}
void VtrCross(const double *vtr1, const double *vtr2, double *cross) {
  cross[0] = vtr1[1] * vtr2[2] - vtr1[2] * vtr2[1];
  cross[1] = vtr1[2] * vtr2[0] - vtr1[0] * vtr2[2];
  cross[2] = vtr1[0] * vtr2[1] - vtr1[1] * vtr2[0];
}
void PrintMat(const double *mat, const size_t row, const size_t col,
              const int precision, const char *name) {
  std::cout << name << std::endl;
  std::cout << std::resetiosflags(std::ios::scientific);
  for (size_t i = 0; i < row; i++) {
    for (size_t j = 0; j < col; j++) {
      std::cout << std::right << std::setw(precision + 7)
                << std::setprecision(precision) << mat[i + j * row];
    }
    std::cout << std::endl;
  }
}
bool SaveMat(const char *path, const double *mat, const size_t &row,
             const size_t &col) {
  std::ofstream stream(path, std::ios::binary);
  if (!stream.is_open()) {
    return false;
  }
  try {
    for (size_t i = 0; i < row; ++i) {
      for (size_t j = 0; j < col; ++j) {
        stream << std::setw(20) << std::setiosflags(std::ios::right)
               << std::setprecision(13) << mat[i + j * row];
      }
      stream << std::endl;
    }
  } catch (std::exception e) {
    stream.close();
    remove(path);
    return false;
  }
  return true;
}
bool ReadMat(const char *path, std::vector<double> *mat, size_t *row,
             size_t *col) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream.is_open()) {
    return false;
  }
  try {
    std::string stmp;
    std::vector<std::string> vtmp;
    getline(stream, stmp);
    function::SplitString(stmp, vtmp);
    *col = vtmp.size();
    *row = 0;
    vtmp.clear();
    stream.seekg(std::ios::beg);
    while (getline(stream, stmp)) {
      function::SplitString(stmp, vtmp);
      if (*col != vtmp.size()) {
        std::vector<double> empty;
        mat->swap(empty);
        return false;
      }
      ++(*row);
      for (auto itr = vtmp.begin(); itr != vtmp.end(); ++itr) {
        mat->push_back(std::stod(*itr));
      }
      vtmp.clear();
    }
  } catch (std::exception) {
    std::vector<double> empty;
    mat->swap(empty);
    return false;
  }
  stream.close();
  double *mat_cp = new double[*row * *col]{0};
  // row major convert to col major
  for (size_t i = 0; i < *col; ++i) {
    for (size_t j = 0; j < *row; ++j) {
      mat_cp[j + i * *row] = (*mat)[i + j * *col];
    }
  }
  memcpy(mat->data(), mat_cp, *row * *col * sizeof(double));
  delete[] mat_cp;
  return true;
}
lapack_int Inv(double *mat, unsigned n) {
  int *ipiv = new int[n + 1];
  lapack_int ret;
  ret = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, mat, n, ipiv);
  if (ret > 0) {
    throw std::runtime_error("Matrix::Inv():singular matrix");
  }
  if (ret < 0) {
    throw std::runtime_error("Matrix::Inv()::llegal argument value");
  }  // TODO:写得不好
  if (ret != 0) return ret;
  ret = LAPACKE_dgetri(LAPACK_COL_MAJOR, n, mat, n, ipiv);
  delete[] ipiv;
  return ret;
}
void Trans(double *a, const size_t &row, const size_t &col) {
  double *cp = new double[col * row]{0};
  for (size_t i = 0; i < row; ++i) {
    for (size_t j = 0; j < col; ++j) {
      cp[j + i * col] = a[i + j * row];
    }
  }
  memcpy(a, cp, row * col * sizeof(double));
  delete[] cp;
}
void Csum(const double *res, double *sum, const size_t &row,
          const size_t &col) {
  for (size_t j = 0; j < col; j++) {
    sum[0 + j * row] = res[0 + j * row];
    for (size_t i = 1; i < row; i++) {
      sum[i + j * row] = sum[i + j * row - 1] + res[i + j * row];
    }
  }
  return;
}
}  // namespace matrix
}  // namespace nstd