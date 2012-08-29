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
//
// === Procedure for adding LAPACK or CBLAS functions to math.cpp/math.h:
//
// This procedure is written as if for a fictional function with double
// version dbacon, which we'll say is from LAPACK.
//
// 1. Write a default templated version which probably returns a boolean.
//    Call it bacon, and put it in math.h.
//
//    Order will always be row-major, so we don't need to pass that.
//    CBLAS_TRANSPOSE-type arguments, however, should be passed.
//
//    Otherwise, arguments should look like those in cblas.h or clapack.h:
//
//    template <typename DType>
//    bool bacon(const CBLAS_TRANSPOSE trans, const int M, const int N, DType* A, ...) {
//      rb_raise(rb_eNotImpError, "only implemented for ATLAS types (float32, float64, complex64, complex128)");
//    }
//
// 2. In math.cpp, add a templated inline static version of the function which takes
//    only void* pointers and uses reinterpret_cast to convert them to the
//    proper dtype.
//
//    This function may also need to switch m and n if these arguments are given.
//
//    For an example, see cblas_gemm. This function should do nothing other than cast
//    appropriately. If clapack_dbacon, clapack_sbacon, clapack_cbacon, and clapack_zbacon
//    all take void* only, and no other pointers that vary between functions, you can skip
//    this particular step -- as we can call them directly using a custom function pointer
//    array (same function signature!).
//
//    This version of the function will be the one exposed through NMatrix::LAPACK. We
//    want it to be as close to the actual LAPACK version of the function as possible,
//    and with as few checks as possible.
//
//    You will probably need a forward declaration in the extern "C" block.
//
//    Note: In that case, the function you wrote in Step 1 should also take exactly the
//    same arguments as clapack_xbacon. Otherwise Bad Things will happen.
//
// 3. In math.cpp, add inline specialized versions of bacon for the different ATLAS types.
//
//    You could do this with a macro, if the arguments are all similar (see #define LAPACK_GETRF).
//    Or you may prefer to do it by hand:
//
//    template <>
//    inline bool bacon(const CBLAS_TRANSPOSE trans, const int M, const int N, float* A, ...) {
//      clapack_sbacon(trans, M, N, A, ...);
//      return true;
//    }
//
//    Make sure these functions are in the namespace nm::math.
//
//    Note that you should do everything in your power here to parse any return values
//    clapack_sbacon may give you. We're not trying very hard in this example, but you might
//    look at getrf to see how it might be done.
//
// 4. Expose the function in nm_math_init_blas(), in math.cpp:
//
//    rb_define_singleton_method(cNMatrix_LAPACK, "clapack_bacon", (METHOD)nm_lapack_bacon, 5);
//
//    Here, we're telling Ruby that nm_lapack_bacon takes five arguments as a Ruby function.
//
// 5. In blas.rb, write a bacon function which accesses clapack_bacon, but does all the
//    sanity checks we left out in step 2.
//
// 6. Write tests for NMatrix::LAPACK::getrf, confirming that it works for the ATLAS dtypes.
//
// 7. After you get it working properly with ATLAS, download dbacon.f from NETLIB, and use
//    f2c to convert it to C. Clean it up so it's readable. Remove the extra indices -- f2c
//    inserts a lot of unnecessary stuff.
//
//    Copy and paste the output into the default templated function you wrote in Step 1.
//    Fix it so it works as a template instead of just for doubles.
//
// 8. Write tests to confirm that it works for integers, rationals, and Ruby objects.
//
// 9. See about adding a Ruby-like interface, such as matrix_matrix_multiply for cblas_gemm,
//    or matrix_vector_multiply for cblas_gemv. This step is not mandatory.
//
// 10. Pull request!

/*
 * Project Includes
 */

#include "math.h"
#include "lapack.h"

#include "nmatrix.h"
#include "ruby_constants.h"

/*
 * Forward Declarations
 */

extern "C" {

  static VALUE nm_cblas_gemm(VALUE self, VALUE trans_a, VALUE trans_b, VALUE m, VALUE n, VALUE k, VALUE vAlpha,
                             VALUE a, VALUE lda, VALUE b, VALUE ldb, VALUE vBeta, VALUE c, VALUE ldc);

  static VALUE nm_cblas_gemv(VALUE self, VALUE trans_a, VALUE m, VALUE n, VALUE vAlpha, VALUE a, VALUE lda,
                             VALUE x, VALUE incx, VALUE vBeta, VALUE y, VALUE incy);

  static VALUE nm_cblas_trsm(VALUE self, VALUE side, VALUE uplo, VALUE trans_a, VALUE diag, VALUE m, VALUE n,
                             VALUE vAlpha, VALUE a, VALUE lda, VALUE b, VALUE ldb);

  static VALUE nm_clapack_getrf(VALUE self, VALUE m, VALUE n, VALUE a, VALUE lda);

  static VALUE nm_clapack_scal(VALUE self, VALUE n, VALUE scale, VALUE vector, VALUE incx);

} // end of extern "C" block

////////////////////
// Math Functions //
////////////////////

namespace nm { namespace math {

/*
 * Calculate the determinant for a dense matrix (A [elements]) of size 2 or 3. Return the result.
 */
template <typename DType>
void det_exact(const int M, const void* A_elements, const int lda, void* result_arg) {
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




/*
 * Function signature conversion for calling CBLAS' gemm functions as directly as possible.
 *
 * For documentation: http://www.netlib.org/blas/dgemm.f
 */
template <typename DType>
inline static bool cblas_gemm(const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_TRANSPOSE trans_b,
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


/*
 * Function signature conversion for calling CBLAS's gemv functions as directly as possible.
 *
 * For documentation: http://www.netlib.org/lapack/double/dgetrf.f
 */
template <typename DType>
inline static bool cblas_gemv(const enum CBLAS_TRANSPOSE trans_a,
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


/*
 * Function signature conversion for calling CBLAS' trsm functions as directly as possible.
 *
 * For documentation: http://www.netlib.org/blas/dtrsm.f
 */
template <typename DType>
inline static void cblas_trsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                               const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_DIAG diag,
                               const int m, const int n, const void* alpha, const void* a,
                               const int lda, void* b, const int ldb)
{
  trsm<DType>(order, side, uplo, trans_a, diag, m, n, *reinterpret_cast<const DType*>(alpha),
              reinterpret_cast<const DType*>(a), lda, reinterpret_cast<DType*>(b), ldb);
}


}} // end of namespace nm::math


extern "C" {

///////////////////
// Ruby Bindings //
///////////////////

void nm_math_init_blas() {
	cNMatrix_LAPACK = rb_define_module_under(cNMatrix, "LAPACK");

  rb_define_singleton_method(cNMatrix_LAPACK, "clapack_getrf", (METHOD)nm_clapack_getrf, 4);
  rb_define_singleton_method(cNMatrix_LAPACK, "clapack_scal", (METHOD)nm_clapack_scal, 4);

  cNMatrix_BLAS = rb_define_module_under(cNMatrix, "BLAS");

	rb_define_singleton_method(cNMatrix_BLAS, "cblas_gemm", (METHOD)nm_cblas_gemm, 13);
	rb_define_singleton_method(cNMatrix_BLAS, "cblas_gemv", (METHOD)nm_cblas_gemv, 11);
	rb_define_singleton_method(cNMatrix_BLAS, "cblas_trsm", (METHOD)nm_cblas_trsm, 11);
}


/* Interprets cblas argument which could be any of false/:no_transpose, :transpose, or :complex_conjugate,
 * into an enum recognized by cblas.
 *
 * Called by nm_cblas_gemm -- basically inline.
 *
 */
static inline enum CBLAS_TRANSPOSE blas_transpose_sym(VALUE op) {
  if (op == Qfalse || rb_to_id(op) == nm_rb_no_transpose) return CblasNoTrans;
  else if (rb_to_id(op) == nm_rb_transpose) return CblasTrans;
  else if (rb_to_id(op) == nm_rb_complex_conjugate) return CblasConjTrans;
  else rb_raise(rb_eArgError, "Expected false, :transpose, or :complex_conjugate");
  return CblasNoTrans;
}

/*
 * Interprets cblas argument which could be :left or :right
 *
 * Called by nm_cblas_trsm -- basically inline
 */
static inline enum CBLAS_SIDE blas_side_sym(VALUE op) {
  ID op_id = rb_to_id(op);
  if (op_id == nm_rb_left)  return CblasLeft;
  if (op_id == nm_rb_right) return CblasRight;
  rb_raise(rb_eArgError, "Expected :left or :right for side argument");
  return CblasLeft;
}

/*
 * Interprets cblas argument which could be :upper or :lower
 *
 * Called by nm_cblas_trsm -- basically inline
 */
static inline enum CBLAS_UPLO blas_uplo_sym(VALUE op) {
  ID op_id = rb_to_id(op);
  if (op_id == nm_rb_upper) return CblasUpper;
  if (op_id == nm_rb_lower) return CblasLower;
  rb_raise(rb_eArgError, "Expected :upper or :lower for uplo argument");
  return CblasUpper;
}


/*
 * Interprets cblas argument which could be :unit (true) or :nonunit (false or anything other than true/:unit)
 *
 * Called by nm_cblas_trsm -- basically inline
 */
static inline enum CBLAS_DIAG blas_diag_sym(VALUE op) {
  if (rb_to_id(op) == nm_rb_unit || op == Qtrue) return CblasUnit;
  return CblasNonUnit;
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
static VALUE nm_cblas_gemm(VALUE self,
                           VALUE trans_a, VALUE trans_b,
                           VALUE m, VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::cblas_gemm, bool, const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_TRANSPOSE trans_b, int m, int n, int k, void* alpha, void* a, int lda, void* b, int ldb, void* beta, void* c, int ldc);

  dtype_t dtype = NM_DTYPE(a);

  void *pAlpha = ALLOCA_N(char, DTYPE_SIZES[dtype]),
       *pBeta  = ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, pAlpha);
  rubyval_to_cval(beta, dtype, pBeta);

  return ttable[dtype](blas_transpose_sym(trans_a), blas_transpose_sym(trans_b), FIX2INT(m), FIX2INT(n), FIX2INT(k), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb), pBeta, NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc)) ? Qtrue : Qfalse;
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
static VALUE nm_cblas_gemv(VALUE self,
                           VALUE trans_a,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE x, VALUE incx,
                           VALUE beta,
                           VALUE y, VALUE incy)
{
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::cblas_gemv, bool, const enum CBLAS_TRANSPOSE trans_a, int m, int n, void* alpha, void* a, int lda, void* x, int incx, void* beta, void* y, int incy);

  dtype_t dtype = NM_DTYPE(a);

  void *pAlpha = ALLOCA_N(char, DTYPE_SIZES[dtype]),
       *pBeta  = ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, pAlpha);
  rubyval_to_cval(beta, dtype, pBeta);

  return ttable[dtype](blas_transpose_sym(trans_a), FIX2INT(m), FIX2INT(n), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx), pBeta, NM_STORAGE_DENSE(y)->elements, FIX2INT(incy)) ? Qtrue : Qfalse;
}


static VALUE nm_cblas_trsm(VALUE self,
                           VALUE side, VALUE uplo,
                           VALUE trans_a, VALUE diag,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb)
{
  static void (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const enum CBLAS_SIDE, const enum CBLAS_UPLO,
                                        const enum CBLAS_TRANSPOSE, const enum CBLAS_DIAG,
                                        const int, const int, const void* alpha, const void* a,
                                        const int lda, void* b, const int ldb) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::cblas_trsm<float>,
      nm::math::cblas_trsm<double>,
      cblas_ctrsm, cblas_ztrsm, // call directly, same function signature!
      nm::math::cblas_trsm<nm::Rational32>,
      nm::math::cblas_trsm<nm::Rational64>,
      nm::math::cblas_trsm<nm::Rational128>,
      nm::math::cblas_trsm<nm::RubyObject>
  };

  dtype_t dtype = NM_DTYPE(a);

  void *pAlpha = ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, pAlpha);

  ttable[dtype](CblasRowMajor, blas_side_sym(side), blas_uplo_sym(uplo), blas_transpose_sym(trans_a), blas_diag_sym(diag), FIX2INT(m), FIX2INT(n), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));

  return Qtrue;
}


/*
 * Based on LAPACK's dscal function, but for any dtype.
 *
 * In-place modification; returns the modified vector as well.
 */
static VALUE nm_clapack_scal(VALUE self, VALUE n, VALUE scale, VALUE vector, VALUE incx) {
  dtype_t dtype = NM_DTYPE(vector);

  void* da      = ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(scale, dtype, da);

  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::clapack_scal, void, const int n, const void* da, void* dx, const int incx);

  ttable[dtype](FIX2INT(n), da, NM_STORAGE_DENSE(vector)->elements, FIX2INT(incx));

  return vector;
}


/* Call any of the clpack_xgetrf functions as directly as possible.
 *
 * The clapack_getrf functions (dgetrf, sgetrf, cgetrf, and zgetrf) compute an LU factorization of a general M-by-N
 * matrix A using partial pivoting with row interchanges.
 *
 * The factorization has the form:
 *    A = P * L * U
 * where P is a permutation matrix, L is lower triangular with unit diagonal elements (lower trapezoidal if m > n),
 * and U is upper triangular (upper trapezoidal if m < n).
 *
 * This is the right-looking level 3 BLAS version of the algorithm.
 *
 * == Arguments
 * See: http://www.netlib.org/lapack/double/dgetrf.f
 * (You don't need argument 5; this is the value returned by this function.)
 *
 * You probably don't want to call this function. Instead, why don't you try clapack_getrf, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 *
 * Returns an array giving the pivot indices (normally these are argument #5).
 */
static VALUE nm_clapack_getrf(VALUE self, VALUE m, VALUE n, VALUE a, VALUE lda) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::clapack_getrf, bool, const int m, const int n, void* a, const int lda, int* ipiv);

  int M = FIX2INT(m),
      N = FIX2INT(n);

  // Allocate the pivot index array, which is of size MIN(M, N).
  size_t ipiv_size = std::min(M,N);
  int* ipiv = ALLOCA_N(int, ipiv_size);

  // Call either our version of getrf or the LAPACK version.
  ttable[NM_DTYPE(a)](FIX2INT(m), FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), ipiv);

  // Result will be stored in a. We return ipiv as an array.
  VALUE ipiv_array = rb_ary_new2(ipiv_size);
  for (size_t i = 0; i < ipiv_size; ++i) {
    rb_ary_store(ipiv_array, i, INT2FIX(ipiv[i]));
  }

  return ipiv_array;
}


/*
 * C accessor for calculating an exact determinant.
 */
void nm_math_det_exact(const int M, const void* elements, const int lda, dtype_t dtype, void* result) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::det_exact, void, const int M, const void* A_elements, const int lda, void* result_arg);

  ttable[dtype](M, elements, lda, result);
}


/*
 * Transpose an array of elements that represent a row-major dense matrix. Does not allocate anything, only does an memcpy.
 */
void nm_math_transpose_generic(const size_t M, const size_t N, const void* A, const int lda, void* B, const int ldb, size_t element_size) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {

      memcpy(reinterpret_cast<char*>(B) + (i*ldb+j)*element_size,
             reinterpret_cast<const char*>(A) + (j*lda+i)*element_size,
             element_size);

    }
  }
}


} // end of extern "C" block
