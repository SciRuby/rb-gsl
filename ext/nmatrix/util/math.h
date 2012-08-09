/////////////////////////////////////////////////////////////////////
// = NMatrix
//
// A linear algebra library for scientific computation in Ruby.
// NMatrix is part of SciRuby.
//
// NMatrix was originally inspired by and derived from NArray, by
// Masahiro Tanaka: http://narray.rubyforge.org
//
// == Copyright Information
//
// SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
// NMatrix is Copyright (c) 2012, Ruby Science Foundation
//
// Please see LICENSE.txt for additional copyright notices.
//
// == Contributing
//
// By contributing source code to SciRuby, you agree to be bound by
// our Contributor Agreement:
//
// * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
//
// == math.h
//
// Header file for math functions, interfacing with BLAS, etc.

#ifndef MATH_H
#define MATH_H

/*
 * Standard Includes
 */
#include <cblas.h>

/*
 * Project Includes
 */
#include "data/complex.h"
#include "data/rational.h"
#include "data/ruby_object.h"
#include "types.h"

/*
 * Macros
 */

/*
 * Types
 */

/*
#ifndef CBLAS_ENUM_DEFINED_H
# define CBLAS_ENUM_DEFINED_H
enum CBLAS_TRANSPOSE { CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, AtlasConj=114 };
enum CBLAS_ORDER { CblasRowMajor=101, CblasColMajor=102 };
#endif
*/

/*
 * Data
 */

/*
 * Functions
 */
template <typename DType>
bool gemm(const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const DType* alpha, const DType* A, const int lda, const DType* B, const int ldb, const DType* beta, DType* C, const int ldc);

template <>
bool gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const float* alpha, const float* A, const int lda, const float* B, const int ldb, const float* beta, float* C, const int ldc);

template <>
bool gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const double* alpha, const double* A, const int lda, const double* B, const int ldb, const double* beta, double* C, const int ldc);

template <>
bool gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const Complex64* alpha, const Complex64* A, const int lda, const Complex64* B, const int ldb, const Complex64* beta, Complex64* C, const int ldc);

template <>
bool gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const Complex128* alpha, const Complex128* A, const int lda, const Complex128* B, const int ldb, const Complex128* beta, Complex128* C, const int ldc);

template <typename DType>
bool gemv(const CBLAS_TRANSPOSE Trans, const int M, const int N, const DType* alpha, const DType* A, const int lda,
          const DType* X, const int incX, const DType* beta, DType* Y, const int incY);

template <>
bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const float* alpha, const float* A, const int lda,
          const float* X, const int incX, const float* beta, float* Y, const int incY);

template <>
bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const double* alpha, const double* A, const int lda,
          const double* X, const int incX, const double* beta, double* Y, const int incY);

template <>
bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const Complex64* alpha, const Complex64* A, const int lda,
          const Complex64* X, const int incX, const Complex64* beta, Complex64* Y, const int incY);

template <>
bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const Complex128* alpha, const Complex128* A, const int lda,
          const Complex128* X, const int incX, const Complex128* beta, Complex128* Y, const int incY);

template <typename DType, typename IType>
void numbmm(const unsigned int n, const unsigned int m, const IType* ia, const IType* ja, const DType* a, const bool diaga,
            const IType* ib, const IType* jb, const DType* b, const bool diagb, IType* ic, IType* jc, DType* c, const bool diagc);

template <typename IType>
void symbmm(const unsigned int n, const unsigned int m, const IType* ia, const IType* ja, const bool diaga,
            const IType* ib, const IType* jb, const bool diagb, IType* ic, const bool diagc);

template <typename DType, typename IType>
void smmp_sort_columns(const unsigned int n, const IType* ia, IType* ja, DType* a);

#endif // MATH_H
