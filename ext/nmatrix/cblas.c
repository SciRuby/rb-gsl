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

//extern enum CBLAS_ORDER;
//extern enum CBLAS_TRANSPOSE;

inline DENSE_PARAM init_cblas_params_for_nm_multiply_matrix(int8_t dtype) {
  DENSE_PARAM p;

  // memset(&(p.beta), 0, sizeof(nm_size128_t)); // set beta to 0

  switch(dtype) {
  case NM_FLOAT32:
  case NM_FLOAT64:
    p.alpha.d[0] = 1.0;
    p.beta.d[0] = 0.0;
    break;

  case NM_COMPLEX64:
    p.alpha.c[0].r = 1.0;
    p.alpha.c[0].i = 0.0;
    p.beta.c[0].r = 0.0;
    p.beta.c[0].i = 0.0;
    break;

  case NM_COMPLEX128:
    p.alpha.z.r = 1.0;
    p.alpha.z.i = 0.0;
    p.beta.z.r = 0.0;
    p.beta.z.i = 0.0;
    break;
  }

  return p;
}


inline void cblas_sgemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  cblas_sgemv(Order, TransA, p.M, p.N, p.alpha.d[0], p.A, p.lda, p.B, p.ldb, p.beta.d[0], p.C, p.ldc);
}

inline void cblas_sgemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  cblas_sgemm(Order, TransA, TransB, p.M, p.N, p.K, p.alpha.d[0], p.A, p.lda, p.B, p.ldb, p.beta.d[0], p.C, p.ldc);
}

inline void cblas_dgemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  cblas_dgemv(Order, TransA, p.M, p.N, p.alpha.d[0], p.A, p.lda, p.B, p.ldb, p.beta.d[0], p.C, p.ldc);
}

inline void cblas_dgemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  cblas_dgemm(Order, TransA, TransB, p.M, p.N, p.K, p.alpha.d[0], p.A, p.lda, p.B, p.ldb, p.beta.d[0], p.C, p.ldc);
}

inline void cblas_cgemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  cblas_cgemv(Order, TransA, p.M, p.N, &(p.alpha.c), p.A, p.lda, p.B, p.ldb, &(p.beta.c), p.C, p.ldc);
}

inline void cblas_cgemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  cblas_cgemm(Order, TransA, TransB, p.M, p.N, p.K, &(p.alpha.c), p.A, p.lda, p.B, p.ldb, &(p.beta.c), p.C, p.ldc);
}

inline void cblas_zgemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p) {
  cblas_zgemv(Order, TransA, p.M, p.N, &(p.alpha.z), p.A, p.lda, p.B, p.ldb, &(p.beta.z), p.C, p.ldc);
}

inline void cblas_zgemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p) {
  cblas_zgemm(Order, TransA, TransB, p.M, p.N, p.K, &(p.alpha.z), p.A, p.lda, p.B, p.ldb, &(p.beta.z), p.C, p.ldc);
}


#endif