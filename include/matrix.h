#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cblas.h"
#include "openblas/lapacke.h"

#include "function.h"

namespace nstd {
namespace matrix {
// 向量变为反对称矩阵
void VtrToAntisymMat(const double *vtr, double *mat);
// 3d vector cross multiplication
void VtrCross(const double *vtr1, const double *vtr2, double *cross);
// print matrix
void PrintMat(const double *mat, const size_t row, const size_t col,
              const int precision = 4, const char *name = "mat");
bool SaveMat(const char *path, const double *mat, const size_t &row,
             const size_t &col);
bool ReadMat(const char *path, std::vector<double> *mat, size_t *row,
             size_t *col);
// inplace inverse n x n matrix A.
// matrix A is Row Major (i.e. firts line, second line ... *not* C[][]
// order) returns:
//   ret = 0 on success
//   ret < 0 illegal argument value
//   ret > 0 singular matrix
lapack_int Inv(double *mat, unsigned n);
// transposition
void Trans(double *a, const size_t &row, const size_t &col);
// matrix accumulative
void Csum(const double *res, double *sum, const size_t &row, const size_t &col);
}  // namespace matrix
}  // namespace nstd

#endif  // _MATRIX_H_