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

/*
 * Project Includes
 */

#include "types.h"

/*
 * Macros
 */

/*
 * Types
 */

/*
 * Data
 */

/*
 * Functions
 */
template <typename LT, typename T>
bool gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const T* alpha, const T* A, const int lda, const T* B, const int ldb, const T* beta, T* C, const int ldc);

template <typename LT, typename T>
bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const T* alpha, const T* A, const int lda,
          const T* X, const int incX, const T* beta, T* Y, const int incY);

template <typename DType, typename IType>
void numbmm(const unsigned int n, const unsigned int m, const IType* ia, const IType* ja, const DType* a, const bool diaga,
            const IType* ib, const IType* jb, const DType* b, const bool diagb, IType* ic, IType* jc, DType* c, const bool diagc);

template <typename IType>
void symbmm(const unsigned int n, const unsigned int m, const IType* ia, const IType* ja, const bool diaga,
            const IType* ib, const IType* jb, const bool diagb, IType* ic, const bool diagc);

template <typename DType, typename IType>
void smmp_sort_columns(const unsigned int n, const IType* ia, IType* ja, DType* a);

#endif // MATH_H
