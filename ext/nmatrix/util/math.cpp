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
// == math.cpp
//
// Ruby-exposed BLAS functions.

/*
 * Project Includes
 */

#include "math.h"

#include "nmatrix.hpp"
#include "ruby_constants.h"

/*
 * Forward Declarations
 */

static VALUE nm_cblas_gemm(VALUE self, VALUE trans_a, VALUE trans_b, VALUE m, VALUE n, VALUE k, VALUE vAlpha,
                           VALUE a, VALUE lda, VALUE b, VALUE ldb, VALUE vBeta, VALUE c, VALUE ldc);

static VALUE nm_cblas_gemv(VALUE self, VALUE trans_a, VALUE m, VALUE n, VALUE vAlpha, VALUE a, VALUE lda,
                           VALUE x, VALUE incx, VALUE vBeta, VALUE y, VALUE incy);


////////////////////
// Math Functions //
////////////////////

/*
 * Calculate the determinant for a dense matrix (A [elements]) of size 2 or 3. Return the result.
 */
template <typename DType>
void det_exact_template(const int M, const void* A_elements, const int lda, void* result_arg) {
  DType* result  = reinterpret_cast<DType*>(result_arg);
  const DType* A = reinterpret_cast<const DType*>(A_elements);

  typename LongDType<DType>::type x, y;

  if (M == 2) {
    *result = A[0] * A[lda+1] - A[1] * A[lda];

  } else if (M == 3) {
    x = A[lda+1] * A[2*lda+2] - A[lda+2] * A[2*lda+1]; // ei - fh
    y = A[lda] * A[2*lda+2] -   A[lda+2] * A[2*lda];   // fg - di
    x = A[0]*x - A[1]*y ; // a*(ei-fh) - b*(fg-di)

    y = A[lda] * A[2*lda+1] - A[lda+1] * A[2*lda];    // dh - eg
    *result = A[2]*y + x; // c*(dh-eg) + _
  } else if (M < 2) {
    rb_raise(rb_eArgError, "can only calculate exact determinant of a square matrix of size 2 or larger");
  } else {
    rb_raise(rb_eNotImpError, "exact determinant calculation needed for matrices larger than 3x3");
  }
}


void det_exact(const int M, const void* elements, const int lda, dtype_t dtype, void* result) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, det_exact_template, void, const int M, const void* A_elements, const int lda, void* result_arg);

  ttable[dtype](M, elements, lda, result);
}


/*
 * Transpose an array of elements that represent a row-major dense matrix. Does not allocate anything, only does an memcpy.
 */
void transpose_generic(const size_t M, const size_t N, const void* A, const int lda, void* B, const int ldb, size_t element_size) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {

      memcpy(reinterpret_cast<char*>(B) + (i*ldb+j)*element_size,
             reinterpret_cast<const char*>(A) + (j*lda+i)*element_size,
             element_size);

    }
  }
}


///////////////////
// Ruby Bindings //
///////////////////

void Init_blas() {
  cBLAS = rb_define_module_under(cNMatrix, "BLAS");

	rb_define_singleton_method(cBLAS, "cblas_gemm", (METHOD)nm_cblas_gemm, 13);
	rb_define_singleton_method(cBLAS, "cblas_gemv", (METHOD)nm_cblas_gemv, 11);
}


/* Interprets cblas argument which could be any of false/:no_transpose, :transpose, or :complex_conjugate,
 * into an enum recognized by cblas.
 *
 * Called by nm_cblas_gemm -- basically inline.
 *
 */
static enum CBLAS_TRANSPOSE gemm_op_sym(VALUE op) {
  if (op == Qfalse || rb_to_id(op) == rbsym_no_transpose) return CblasNoTrans;
  else if (rb_to_id(op) == rbsym_transpose) return CblasTrans;
  else if (rb_to_id(op) == rbsym_complex_conjugate) return CblasConjTrans;
  else rb_raise(rb_eArgError, "Expected false, :transpose, or :complex_conjugate");
  return CblasNoTrans;
}


/* Call any of the cblas_xgemm functions as directly as possible.
 *
 * The cblas_xgemm functions (dgemm, sgemm, cgemm, and zgemm) define the following operation:
 *
 *    C = alpha*op(A)*op(B) + beta*C
 *
 * where op(X) is one of <tt>op(X) = X</tt>, <tt>op(X) = X**T</tt>, or the complex conjugate of X.
 *
 * Note that this will only work for dense matrices that are of types :float32, :float64, :complex64, and :complex128.
 * Other types are not implemented in BLAS, and while they exist in NMatrix, this method is intended only to
 * expose the ultra-optimized ATLAS versions.
 *
 * == Arguments
 * See: http://www.netlib.org/blas/dgemm.f
 *
 * You probably don't want to call this function. Instead, why don't you try cblas_gemm, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
template <typename DType>
inline static bool nm_cblas_gemm_template(dtype_t dtype,
                                   const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_TRANSPOSE trans_b,
                                   int m, int n, int k,
                                   void* alpha,
                                   void* a, int lda,
                                   void* b, int ldb,
                                   void* beta,
                                   void* c, int ldc)
{
  return  gemm<DType>(trans_a, trans_b, n, m, k, reinterpret_cast<DType*>(alpha),
                      reinterpret_cast<DType*>(b), ldb,
                      reinterpret_cast<DType*>(a), lda, reinterpret_cast<DType*>(beta),
                      reinterpret_cast<DType*>(c), ldc);
}


static VALUE nm_cblas_gemm(VALUE self,
                           VALUE trans_a, VALUE trans_b,
                           VALUE m, VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm_cblas_gemm_template, bool, dtype_t dtype, const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_TRANSPOSE trans_b, int m, int n, int k, void* alpha, void* a, int lda, void* b, int ldb, void* beta, void* c, int ldc);

  dtype_t dtype = NM_DTYPE(a);

  void *pAlpha = ALLOCA_N(char, DTYPE_SIZES[dtype]),
       *pBeta  = ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, pAlpha);
  rubyval_to_cval(beta, dtype, pBeta);

  return ttable[dtype](dtype, gemm_op_sym(trans_a), gemm_op_sym(trans_b), FIX2INT(m), FIX2INT(n), FIX2INT(k), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb), pBeta, NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc)) ? Qtrue : Qfalse;
}


/* Call any of the cblas_xgemv functions as directly as possible.
 *
 * The cblas_xgemv functions (dgemv, sgemv, cgemv, and zgemv) define the following operation:
 *
 *    y = alpha*op(A)*x + beta*y
 *
 * where op(A) is one of <tt>op(A) = A</tt>, <tt>op(A) = A**T</tt>, or the complex conjugate of A.
 *
 * Note that this will only work for dense matrices that are of types :float32, :float64, :complex64, and :complex128.
 * Other types are not implemented in BLAS, and while they exist in NMatrix, this method is intended only to
 * expose the ultra-optimized ATLAS versions.
 *
 * == Arguments
 * See: http://www.netlib.org/blas/dgemm.f
 *
 * You probably don't want to call this function. Instead, why don't you try cblas_gemv, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
template <typename DType>
inline static bool nm_cblas_gemv_template(dtype_t dtype,
                                    const enum CBLAS_TRANSPOSE trans_a,
                                    int m, int n,
                                    void* alpha,
                                    void* a, int lda,
                                    void* x, int incx,
                                    void* beta,
                                    void* y, int incy)
{
  return gemv<DType>(trans_a,
                     m, n, reinterpret_cast<DType*>(alpha),
                     reinterpret_cast<DType*>(a), lda,
                     reinterpret_cast<DType*>(x), incx, reinterpret_cast<DType*>(beta),
                     reinterpret_cast<DType*>(y), incy);
}

static VALUE nm_cblas_gemv(VALUE self,
                           VALUE trans_a,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE x, VALUE incx,
                           VALUE beta,
                           VALUE y, VALUE incy)
{
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm_cblas_gemv_template, bool, dtype_t dtype, const enum CBLAS_TRANSPOSE trans_a, int m, int n, void* alpha, void* a, int lda, void* x, int incx, void* beta, void* y, int incy);

  dtype_t dtype = NM_DTYPE(a);

  void *pAlpha = ALLOCA_N(char, DTYPE_SIZES[dtype]),
       *pBeta  = ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, pAlpha);
  rubyval_to_cval(beta, dtype, pBeta);

  return ttable[dtype](dtype, gemm_op_sym(trans_a), FIX2INT(m), FIX2INT(n), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx), pBeta, NM_STORAGE_DENSE(y)->elements, FIX2INT(incy)) ? Qtrue : Qfalse;
}



