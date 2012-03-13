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
// == cblas.c
//
// Functions in this file allow us to call CBLAS functions using
// arrays of function pointers, by ensuring that each has the same
// signature.

#ifndef CBLAS_C
# define CBLAS_C

#include <cblas.h>
#include "nmatrix.h"

//extern const enum CBLAS_ORDER;
//extern const enum CBLAS_TRANSPOSE;


inline void cblas_r32gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  if (Order == CblasColMajor) r32gemm(TransA, TransB, p.M, p.N, p.K, p.alpha.r[0], p.A, p.lda, p.B, p.ldb, p.beta.r[0], p.C, p.ldc);
  else                        r32gemm(TransB, TransA, p.N, p.M, p.K, p.alpha.r[0], p.B, p.ldb, p.A, p.lda, p.beta.r[0], p.C, p.ldc);
}

inline void cblas_r64gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  if (Order == CblasColMajor) r64gemm(TransA, TransB, p.M, p.N, p.K, p.alpha.ra[0], p.A, p.lda, p.B, p.ldb, p.beta.ra[0], p.C, p.ldc);
  else                        r64gemm(TransB, TransA, p.N, p.M, p.K, p.alpha.ra[0], p.B, p.ldb, p.A, p.lda, p.beta.ra[0], p.C, p.ldc);
}

inline void cblas_r128gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  if (Order == CblasColMajor) r128gemm(TransA, TransB, p.M, p.N, p.K, p.alpha.rat, p.A, p.lda, p.B, p.ldb, p.beta.rat, p.C, p.ldc);
  else                        r128gemm(TransB, TransA, p.N, p.M, p.K, p.alpha.rat, p.B, p.ldb, p.A, p.lda, p.beta.rat, p.C, p.ldc);
}

inline void cblas_bgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  bgemv(TransA, p.M, p.N, p.alpha.b[0], p.A, p.lda, p.B, p.ldb, p.beta.b[0], p.C, p.ldc);
}

inline void cblas_bgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  if (Order == CblasColMajor) bgemm(TransA, TransB, p.M, p.N, p.K, p.alpha.b[0], p.A, p.lda, p.B, p.ldb, p.beta.b[0], p.C, p.ldc);
  else                        bgemm(TransB, TransA, p.N, p.M, p.K, p.alpha.b[0], p.B, p.ldb, p.A, p.lda, p.beta.b[0], p.C, p.ldc);
}

inline void cblas_i8gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  i8gemv(TransA, p.M, p.N, p.alpha.i[0], p.A, p.lda, p.B, p.ldb, p.beta.i[0], p.C, p.ldc);
}

inline void cblas_i8gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  if (Order == CblasColMajor) i8gemm(TransA, TransB, p.M, p.N, p.K, p.alpha.i[0], p.A, p.lda, p.B, p.ldb, p.beta.i[0], p.C, p.ldc);
  else                        i8gemm(TransB, TransA, p.N, p.M, p.K, p.alpha.i[0], p.B, p.ldb, p.A, p.lda, p.beta.i[0], p.C, p.ldc);
}

inline void cblas_i16gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  i16gemv(TransA, p.M, p.N, p.alpha.i[0], p.A, p.lda, p.B, p.ldb, p.beta.i[0], p.C, p.ldc);
}

inline void cblas_i16gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  if (Order == CblasColMajor) i16gemm(TransA, TransB, p.M, p.N, p.K, p.alpha.i[0], p.A, p.lda, p.B, p.ldb, p.beta.i[0], p.C, p.ldc);
  else                        i16gemm(TransB, TransA, p.N, p.M, p.K, p.alpha.i[0], p.B, p.ldb, p.A, p.lda, p.beta.i[0], p.C, p.ldc);
}

inline void cblas_i32gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  i32gemv(TransA, p.M, p.N, p.alpha.i[0], p.A, p.lda, p.B, p.ldb, p.beta.i[0], p.C, p.ldc);
}

inline void cblas_i32gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  if (Order == CblasColMajor) i32gemm(TransA, TransB, p.M, p.N, p.K, p.alpha.i[0], p.A, p.lda, p.B, p.ldb, p.beta.i[0], p.C, p.ldc);
  else                        i32gemm(TransB, TransA, p.N, p.M, p.K, p.alpha.i[0], p.B, p.ldb, p.A, p.lda, p.beta.i[0], p.C, p.ldc);
}

inline void cblas_i64gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  i64gemv(TransA, p.M, p.N, p.alpha.i[0], p.A, p.lda, p.B, p.ldb, p.beta.i[0], p.C, p.ldc);
}

inline void cblas_i64gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  if (Order == CblasColMajor) i64gemm(TransA, TransB, p.M, p.N, p.K, p.alpha.i[0], p.A, p.lda, p.B, p.ldb, p.beta.i[0], p.C, p.ldc);
  else                        i64gemm(TransB, TransA, p.N, p.M, p.K, p.alpha.i[0], p.B, p.ldb, p.A, p.lda, p.beta.i[0], p.C, p.ldc);
}

inline void cblas_vgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  if (Order == CblasColMajor) vgemm(TransA, TransB, p.M, p.N, p.K, p.alpha.v[0], p.A, p.lda, p.B, p.ldb, p.beta.v[0], p.C, p.ldc);
  else                        vgemm(TransB, TransA, p.N, p.M, p.K, p.alpha.v[0], p.B, p.ldb, p.A, p.lda, p.beta.v[0], p.C, p.ldc);
}

inline void cblas_sgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  cblas_sgemv(Order, TransA, p.M, p.N, p.alpha.d[0], p.A, p.lda, p.B, p.ldb, p.beta.d[0], p.C, p.ldc);
}

inline void cblas_sgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  cblas_sgemm(Order, TransA, TransB, p.M, p.N, p.K, p.alpha.d[0], p.A, p.lda, p.B, p.ldb, p.beta.d[0], p.C, p.ldc);
}

inline void cblas_dgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  cblas_dgemv(Order, TransA, p.M, p.N, p.alpha.d[0], p.A, p.lda, p.B, p.ldb, p.beta.d[0], p.C, p.ldc);
}

inline void cblas_dgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  cblas_dgemm(Order, TransA, TransB, p.M, p.N, p.K, p.alpha.d[0], p.A, p.lda, p.B, p.ldb, p.beta.d[0], p.C, p.ldc);
}

inline void cblas_cgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  cblas_cgemv(Order, TransA, p.M, p.N, &(p.alpha.c), p.A, p.lda, p.B, p.ldb, &(p.beta.c), p.C, p.ldc);
}

inline void cblas_cgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  cblas_cgemm(Order, TransA, TransB, p.M, p.N, p.K, &(p.alpha.c), p.A, p.lda, p.B, p.ldb, &(p.beta.c), p.C, p.ldc);
}

inline void cblas_zgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  cblas_zgemv(Order, TransA, p.M, p.N, &(p.alpha.z), p.A, p.lda, p.B, p.ldb, &(p.beta.z), p.C, p.ldc);
}

inline void cblas_zgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  cblas_zgemm(Order, TransA, TransB, p.M, p.N, p.K, &(p.alpha.z), p.A, p.lda, p.B, p.ldb, &(p.beta.z), p.C, p.ldc);
}


#endif