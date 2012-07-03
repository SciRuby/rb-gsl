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
// == cblas.h
//
// Singly-linked list implementation

#ifndef NMATRIX_CBLAS_H
#define NMATRIX_CBLAS_H

/*
 * Standard Includes
 */
 
#include <cblas.h>

/*
 * Project Includes
 */

#include "nmatrix.h"

/*
 * Macros
 */

// BLAS functions
#define SMMP_MAX_THREE(a,b,c) ((a)>(b) ? ( (a)>(c) ? (a) : (c) ) : ( (b)>(c) ? (b) : (c) ))
#define SMMP_MIN(a,b) ((a)>(b) ? (b) : (a))
#define SMMP_MAX(a,b) ((a)>(b) ? (a) : (b))

/*
 * Types
 */

// For calling cblas_gemm functions.
typedef struct {
  int M, N, K, lda, ldb, ldc;
  void *A, *B, *C;
  nm_size128_t alpha, beta;
} DENSE_PARAM;


typedef struct {
  void *ia, *ja, *a;
  bool diag;
} YALE_PARAM;
																																	// dense transpose
typedef void     (*nm_gemv_t[NUM_DTYPES])();																																						// general matrix/vector multiply
typedef void     (*nm_smmp_t[NUM_DTYPES][NM_INDEX_TYPES])();																														// sparse (yale) multiply
typedef void     (*nm_smmp_transpose_t[NUM_DTYPES][NM_INDEX_TYPES])(y_size_t, y_size_t, YALE_PARAM, YALE_PARAM, bool);	// sparse (yale) transpose

/*
 * Data
 */

extern int (*Gemm[15])(const enum CBLAS_TRANSPOSE, const enum CBLAS_TRANSPOSE, const int, const int, const int, const void *, const void *, const int, const void *, const int, const void *, void *, const int);
extern int (*Gemv[15])(const enum CBLAS_TRANSPOSE, const int, const int, const void *, const void *, const int, const void *, const int, const void *, void *, const int);
extern void (*Symbmm[7])(const unsigned int, const unsigned int, const void *, const void *, const bool, const void *, const void *, const bool, void *, const bool);
extern void (*Numbmm[15][7])(const unsigned int, const unsigned int, const void *, const void *, const void *, const bool, const void *, const void *, const void *, const bool, void *, void *, void *, const bool);

/*
 * Functions
 */

DENSE_PARAM init_cblas_params_for_nm_multiply_matrix(dtype_t dtype);
void cblas_bgemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_bgemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_i8gemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_i8gemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_i16gemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_i16gemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_i32gemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_i32gemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_i64gemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_i64gemv_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_sgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_sgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_dgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_dgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_cgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_cgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_zgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_zgemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_r32gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_r32gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_r64gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_r64gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_r128gemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);
void cblas_r128gemv_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, DENSE_PARAM p);
void cblas_vgemm_(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, DENSE_PARAM p);

#endif
