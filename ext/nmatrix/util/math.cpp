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
#ifdef HAVE_CLAPACK_H
  #include <clapack.h>
#endif

  static VALUE nm_cblas_rot(VALUE self, VALUE n, VALUE x, VALUE incx, VALUE y, VALUE incy, VALUE c, VALUE s);
  static VALUE nm_cblas_rotg(VALUE self, VALUE ab);

  static VALUE nm_cblas_gemm(VALUE self, VALUE order, VALUE trans_a, VALUE trans_b, VALUE m, VALUE n, VALUE k, VALUE vAlpha,
                             VALUE a, VALUE lda, VALUE b, VALUE ldb, VALUE vBeta, VALUE c, VALUE ldc);
  static VALUE nm_cblas_gemv(VALUE self, VALUE trans_a, VALUE m, VALUE n, VALUE vAlpha, VALUE a, VALUE lda,
                             VALUE x, VALUE incx, VALUE vBeta, VALUE y, VALUE incy);
  static VALUE nm_cblas_trsm(VALUE self, VALUE order, VALUE side, VALUE uplo, VALUE trans_a, VALUE diag, VALUE m, VALUE n,
                             VALUE vAlpha, VALUE a, VALUE lda, VALUE b, VALUE ldb);
  static VALUE nm_cblas_trmm(VALUE self, VALUE order, VALUE side, VALUE uplo, VALUE trans_a, VALUE diag, VALUE m, VALUE n,
                             VALUE alpha, VALUE a, VALUE lda, VALUE b, VALUE ldb);
  static VALUE nm_cblas_herk(VALUE self, VALUE order, VALUE uplo, VALUE trans, VALUE n, VALUE k, VALUE alpha, VALUE a,
                             VALUE lda, VALUE beta, VALUE c, VALUE ldc);
  static VALUE nm_cblas_syrk(VALUE self, VALUE order, VALUE uplo, VALUE trans, VALUE n, VALUE k, VALUE alpha, VALUE a,
                             VALUE lda, VALUE beta, VALUE c, VALUE ldc);

  static VALUE nm_clapack_getrf(VALUE self, VALUE order, VALUE m, VALUE n, VALUE a, VALUE lda);
  static VALUE nm_clapack_potrf(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda);
  static VALUE nm_clapack_getrs(VALUE self, VALUE order, VALUE trans, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE ipiv, VALUE b, VALUE ldb);
  static VALUE nm_clapack_potrs(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE b, VALUE ldb);
  static VALUE nm_clapack_getri(VALUE self, VALUE order, VALUE n, VALUE a, VALUE lda, VALUE ipiv);
  static VALUE nm_clapack_potri(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda);
  static VALUE nm_clapack_laswp(VALUE self, VALUE n, VALUE a, VALUE lda, VALUE k1, VALUE k2, VALUE ipiv, VALUE incx);
  static VALUE nm_clapack_scal(VALUE self, VALUE n, VALUE scale, VALUE vector, VALUE incx);
  static VALUE nm_clapack_lauum(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda);

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
inline static void cblas_gemm(const enum CBLAS_ORDER order,
                              const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_TRANSPOSE trans_b,
                              int m, int n, int k,
                              void* alpha,
                              void* a, int lda,
                              void* b, int ldb,
                              void* beta,
                              void* c, int ldc)
{
  gemm<DType>(order, trans_a, trans_b, m, n, k, reinterpret_cast<DType*>(alpha),
              reinterpret_cast<DType*>(a), lda,
              reinterpret_cast<DType*>(b), ldb, reinterpret_cast<DType*>(beta),
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


/*
 * Function signature conversion for calling CBLAS' trmm functions as directly as possible.
 *
 * For documentation: http://www.netlib.org/blas/dtrmm.f
 */
template <typename DType>
inline static void cblas_trmm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side, const enum CBLAS_UPLO uplo,
                              const enum CBLAS_TRANSPOSE ta, const enum CBLAS_DIAG diag, const int m, const int n, const void* alpha,
                              const void* A, const int lda, void* B, const int ldb)
{
  trmm<DType>(order, side, uplo, ta, diag, m, n, reinterpret_cast<const DType*>(alpha),
              reinterpret_cast<const DType*>(A), lda, reinterpret_cast<DType*>(B), ldb);
}


/*
 * Function signature conversion for calling CBLAS' syrk functions as directly as possible.
 *
 * For documentation: http://www.netlib.org/blas/dsyrk.f
 */
template <typename DType>
inline static void cblas_syrk(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const enum CBLAS_TRANSPOSE trans,
                              const int n, const int k, const void* alpha,
                              const void* A, const int lda, const void* beta, void* C, const int ldc)
{
  syrk<DType>(order, uplo, trans, n, k, reinterpret_cast<const DType*>(alpha),
              reinterpret_cast<const DType*>(A), lda, reinterpret_cast<const DType*>(beta), reinterpret_cast<DType*>(C), ldc);
}




}} // end of namespace nm::math


extern "C" {

///////////////////
// Ruby Bindings //
///////////////////

void nm_math_init_blas() {
	cNMatrix_LAPACK = rb_define_module_under(cNMatrix, "LAPACK");

  rb_define_singleton_method(cNMatrix_LAPACK, "clapack_getrf", (METHOD)nm_clapack_getrf, 5);
  rb_define_singleton_method(cNMatrix_LAPACK, "clapack_potrf", (METHOD)nm_clapack_potrf, 5);
  rb_define_singleton_method(cNMatrix_LAPACK, "clapack_getrs", (METHOD)nm_clapack_getrs, 9);
  rb_define_singleton_method(cNMatrix_LAPACK, "clapack_potrs", (METHOD)nm_clapack_potrs, 8);
  rb_define_singleton_method(cNMatrix_LAPACK, "clapack_getri", (METHOD)nm_clapack_getri, 5);
  rb_define_singleton_method(cNMatrix_LAPACK, "clapack_potri", (METHOD)nm_clapack_potri, 5);
  rb_define_singleton_method(cNMatrix_LAPACK, "clapack_laswp", (METHOD)nm_clapack_laswp, 7);
  rb_define_singleton_method(cNMatrix_LAPACK, "clapack_scal",  (METHOD)nm_clapack_scal,  4);
  rb_define_singleton_method(cNMatrix_LAPACK, "clapack_lauum", (METHOD)nm_clapack_lauum, 5);

  cNMatrix_BLAS = rb_define_module_under(cNMatrix, "BLAS");

  rb_define_singleton_method(cNMatrix_BLAS, "cblas_rot",  (METHOD)nm_cblas_rot,  7);
  rb_define_singleton_method(cNMatrix_BLAS, "cblas_rotg", (METHOD)nm_cblas_rotg, 1);

	rb_define_singleton_method(cNMatrix_BLAS, "cblas_gemm", (METHOD)nm_cblas_gemm, 14);
	rb_define_singleton_method(cNMatrix_BLAS, "cblas_gemv", (METHOD)nm_cblas_gemv, 11);
	rb_define_singleton_method(cNMatrix_BLAS, "cblas_trsm", (METHOD)nm_cblas_trsm, 12);
	rb_define_singleton_method(cNMatrix_BLAS, "cblas_trmm", (METHOD)nm_cblas_trmm, 12);
	rb_define_singleton_method(cNMatrix_BLAS, "cblas_syrk", (METHOD)nm_cblas_syrk, 11);
	rb_define_singleton_method(cNMatrix_BLAS, "cblas_herk", (METHOD)nm_cblas_herk, 11);
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

/*
 * Interprets cblas argument which could be :row or :col
 */
static inline enum CBLAS_ORDER blas_order_sym(VALUE op) {
  if (rb_to_id(op) == rb_intern("row") || rb_to_id(op) == rb_intern("row_major")) return CblasRowMajor;
  else if (rb_to_id(op) == rb_intern("col") || rb_to_id(op) == rb_intern("col_major") ||
           rb_to_id(op) == rb_intern("column") || rb_to_id(op) == rb_intern("column_major")) return CblasColMajor;
  rb_raise(rb_eArgError, "Expected :row or :col for order argument");
  return CblasRowMajor;
}


/*
 * Call any of the cblas_xrotg functions as directly as possible.
 *
 * xROTG computes the elements of a Givens plane rotation matrix such that:
 *
 *  |  c s |   | a |   | r |
 *  | -s c | * | b | = | 0 |
 *
 * where r = +- sqrt( a**2 + b**2 ) and c**2 + s**2 = 1.
 *
 * The Givens plane rotation can be used to introduce zero elements into a matrix selectively.
 *
 * This function differs from most of the other raw BLAS accessors. Instead of providing a, b, c, s as arguments, you
 * should only provide a and b (the inputs), and you should provide them as a single NVector (or the first two elements
 * of any dense NMatrix or NVector type, specifically).
 *
 * The outputs [c,s] will be returned in a Ruby Array at the end; the input NVector will also be modified in-place.
 *
 * If you provide rationals, be aware that there's a high probability of an error, since rotg includes a square root --
 * and most rationals' square roots are irrational. You're better off converting to Float first.
 *
 * This function, like the other cblas_ functions, does minimal type-checking.
 */
static VALUE nm_cblas_rotg(VALUE self, VALUE ab) {
  static void (*ttable[nm::NUM_DTYPES])(void* a, void* b, void* c, void* s) = {
      NULL, NULL, NULL, NULL, NULL, // can't represent c and s as integers, so no point in having integer operations.
      nm::math::cblas_rotg<float>,
      nm::math::cblas_rotg<double>,
      nm::math::cblas_rotg<nm::Complex64>,
      nm::math::cblas_rotg<nm::Complex128>,
      nm::math::cblas_rotg<nm::Rational32>,
      nm::math::cblas_rotg<nm::Rational64>,
      nm::math::cblas_rotg<nm::Rational128>,
      nm::math::cblas_rotg<nm::RubyObject>
  };

  dtype_t dtype = NM_DTYPE(ab);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
    return Qnil;

  } else {
    void *pC = ALLOCA_N(char, DTYPE_SIZES[dtype]),
         *pS = ALLOCA_N(char, DTYPE_SIZES[dtype]);

    // extract A and B from the NVector (first two elements)
    void* pA = NM_STORAGE_DENSE(ab)->elements;
    void* pB = (char*)(NM_STORAGE_DENSE(ab)->elements) + DTYPE_SIZES[dtype];
    // c and s are output

    ttable[dtype](pA, pB, pC, pS);

    VALUE result = rb_ary_new2(2);
    rb_ary_store(result, 0, rubyobj_from_cval(pC, dtype).rval);
    rb_ary_store(result, 1, rubyobj_from_cval(pS, dtype).rval);

    return result;
  }
}


/*
 * Call any of the cblas_xrot functions as directly as possible.
 *
 * xROT is a BLAS level 1 routine (taking two vectors) which applies a plane rotation.
 *
 * It's tough to find documentation on xROT. Here are what we think the arguments are for:
 *  * n     :: number of elements to consider in x and y
 *  * x     :: a vector (expects an NVector)
 *  * incx  :: stride of x
 *  * y     :: a vector (expects an NVector)
 *  * incy  :: stride of y
 *  * c     :: cosine of the angle of rotation
 *  * s     :: sine of the angle of rotation
 *
 * Note that c and s will be the same dtype as x and y, except when x and y are complex. If x and y are complex, c and s
 * will be float for Complex64 or double for Complex128.
 *
 * You probably don't want to call this function. Instead, why don't you try rot, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_cblas_rot(VALUE self, VALUE n, VALUE x, VALUE incx, VALUE y, VALUE incy, VALUE c, VALUE s) {
  static void (*ttable[nm::NUM_DTYPES])(const int N, void*, const int, void*, const int, const void*, const void*) = {
      NULL, NULL, NULL, NULL, NULL, // can't represent c and s as integers, so no point in having integer operations.
      nm::math::cblas_rot<float,float>,
      nm::math::cblas_rot<double,double>,
      nm::math::cblas_rot<nm::Complex64,float>,
      nm::math::cblas_rot<nm::Complex128,double>,
      nm::math::cblas_rot<nm::Rational32,nm::Rational32>,
      nm::math::cblas_rot<nm::Rational64,nm::Rational64>,
      nm::math::cblas_rot<nm::Rational128,nm::Rational128>,
      nm::math::cblas_rot<nm::RubyObject,nm::RubyObject>
  };

  dtype_t dtype = NM_DTYPE(x);


  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
    return Qfalse;
  } else {
    void *pC, *pS;

    // We need to ensure the cosine and sine arguments are the correct dtype -- which may differ from the actual dtype.
    if (dtype == COMPLEX64) {
      pC = ALLOCA_N(float,1);
      pS = ALLOCA_N(float,1);
      rubyval_to_cval(c, FLOAT32, pC);
      rubyval_to_cval(s, FLOAT32, pS);
    } else if (dtype == COMPLEX128) {
      pC = ALLOCA_N(double,1);
      pS = ALLOCA_N(double,1);
      rubyval_to_cval(c, FLOAT64, pC);
      rubyval_to_cval(s, FLOAT64, pS);
    } else {
      pC = ALLOCA_N(char, DTYPE_SIZES[dtype]);
      pS = ALLOCA_N(char, DTYPE_SIZES[dtype]);
      rubyval_to_cval(c, dtype, pC);
      rubyval_to_cval(s, dtype, pS);
    }


    ttable[dtype](FIX2INT(n), NM_STORAGE_DENSE(x)->elements, FIX2INT(incx), NM_STORAGE_DENSE(y)->elements, FIX2INT(incy), pC, pS);

    return Qtrue;
  }
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
 * You probably don't want to call this function. Instead, why don't you try gemm, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 */
static VALUE nm_cblas_gemm(VALUE self,
                           VALUE order,
                           VALUE trans_a, VALUE trans_b,
                           VALUE m, VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, nm::math::cblas_gemm, void, const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE trans_a, const enum CBLAS_TRANSPOSE trans_b, int m, int n, int k, void* alpha, void* a, int lda, void* b, int ldb, void* beta, void* c, int ldc);

  dtype_t dtype = NM_DTYPE(a);

  void *pAlpha = ALLOCA_N(char, DTYPE_SIZES[dtype]),
       *pBeta  = ALLOCA_N(char, DTYPE_SIZES[dtype]);
  rubyval_to_cval(alpha, dtype, pAlpha);
  rubyval_to_cval(beta, dtype, pBeta);

  ttable[dtype](blas_order_sym(order), blas_transpose_sym(trans_a), blas_transpose_sym(trans_b), FIX2INT(m), FIX2INT(n), FIX2INT(k), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb), pBeta, NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc));

  return c;
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
                           VALUE order,
                           VALUE side, VALUE uplo,
                           VALUE trans_a, VALUE diag,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb)
{
  static void (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const enum CBLAS_SIDE, const enum CBLAS_UPLO,
                                        const enum CBLAS_TRANSPOSE, const enum CBLAS_DIAG,
                                        const int m, const int n, const void* alpha, const void* a,
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

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    void *pAlpha = ALLOCA_N(char, DTYPE_SIZES[dtype]);
    rubyval_to_cval(alpha, dtype, pAlpha);

    ttable[dtype](blas_order_sym(order), blas_side_sym(side), blas_uplo_sym(uplo), blas_transpose_sym(trans_a), blas_diag_sym(diag), FIX2INT(m), FIX2INT(n), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));
  }

  return Qtrue;
}


static VALUE nm_cblas_trmm(VALUE self,
                           VALUE order,
                           VALUE side, VALUE uplo,
                           VALUE trans_a, VALUE diag,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb)
{
  static void (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER,
                                        const enum CBLAS_SIDE, const enum CBLAS_UPLO,
                                        const enum CBLAS_TRANSPOSE, const enum CBLAS_DIAG,
                                        const int m, const int n, const void* alpha, const void* a,
                                        const int lda, void* b, const int ldb) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::cblas_trmm<float>,
      nm::math::cblas_trmm<double>,
      cblas_ctrmm, cblas_ztrmm // call directly, same function signature!
      /*
      nm::math::cblas_trmm<nm::Rational32>,
      nm::math::cblas_trmm<nm::Rational64>,
      nm::math::cblas_trmm<nm::Rational128>,
      nm::math::cblas_trmm<nm::RubyObject>*/
  };

  dtype_t dtype = NM_DTYPE(a);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this matrix operation not yet defined for non-BLAS dtypes");
  } else {
    void *pAlpha = ALLOCA_N(char, DTYPE_SIZES[dtype]);
    rubyval_to_cval(alpha, dtype, pAlpha);

    ttable[dtype](blas_order_sym(order), blas_side_sym(side), blas_uplo_sym(uplo), blas_transpose_sym(trans_a), blas_diag_sym(diag), FIX2INT(m), FIX2INT(n), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));
  }

  return b;
}


static VALUE nm_cblas_syrk(VALUE self,
                           VALUE order,
                           VALUE uplo,
                           VALUE trans,
                           VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{
  static void (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const enum CBLAS_UPLO, const enum CBLAS_TRANSPOSE,
                                        const int n, const int k, const void* alpha, const void* a,
                                        const int lda, const void* beta, void* c, const int ldc) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::cblas_syrk<float>,
      nm::math::cblas_syrk<double>,
      cblas_csyrk, cblas_zsyrk// call directly, same function signature!
      /*nm::math::cblas_trsm<nm::Rational32>,
      nm::math::cblas_trsm<nm::Rational64>,
      nm::math::cblas_trsm<nm::Rational128>,
      nm::math::cblas_trsm<nm::RubyObject>*/
  };

  dtype_t dtype = NM_DTYPE(a);

  if (!ttable[dtype]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    void *pAlpha = ALLOCA_N(char, DTYPE_SIZES[dtype]),
         *pBeta = ALLOCA_N(char, DTYPE_SIZES[dtype]);
    rubyval_to_cval(alpha, dtype, pAlpha);
    rubyval_to_cval(beta, dtype, pBeta);

    ttable[dtype](blas_order_sym(order), blas_uplo_sym(uplo), blas_transpose_sym(trans), FIX2INT(n), FIX2INT(k), pAlpha, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), pBeta, NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc));
  }

  return Qtrue;
}


static VALUE nm_cblas_herk(VALUE self,
                           VALUE order,
                           VALUE uplo,
                           VALUE trans,
                           VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{

  dtype_t dtype = NM_DTYPE(a);

  if (dtype == COMPLEX64) {
    cblas_cherk(blas_order_sym(order), blas_uplo_sym(uplo), blas_transpose_sym(trans), FIX2INT(n), FIX2INT(k), NUM2DBL(alpha), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NUM2DBL(beta), NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc));
  } else if (dtype == COMPLEX128) {
    cblas_zherk(blas_order_sym(order), blas_uplo_sym(uplo), blas_transpose_sym(trans), FIX2INT(n), FIX2INT(k), NUM2DBL(alpha), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), NUM2DBL(beta), NM_STORAGE_DENSE(c)->elements, FIX2INT(ldc));
  } else
    rb_raise(rb_eNotImpError, "this matrix operation undefined for non-complex dtypes");


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


static VALUE nm_clapack_lauum(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const enum CBLAS_UPLO, const int n, void* a, const int lda) = {
      /*nm::math::clapack_lauum<uint8_t, false>,
      nm::math::clapack_lauum<int8_t, false>,
      nm::math::clapack_lauum<int16_t, false>,
      nm::math::clapack_lauum<uint32_t, false>,
      nm::math::clapack_lauum<uint64_t, false>,*/
      NULL, NULL, NULL, NULL, NULL,
      nm::math::clapack_lauum<false, float>,
      nm::math::clapack_lauum<false, double>,
#ifdef HAVE_CLAPACK_H
      clapack_clauum, clapack_zlauum, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::clapack_lauum<true, nm::Complex64>,
      nm::math::clapack_lauum<true, nm::Complex128>,
#endif
/*
      nm::math::clapack_lauum<nm::Rational32, false>,
      nm::math::clapack_lauum<nm::Rational64, false>,
      nm::math::clapack_lauum<nm::Rational128, false>,
      nm::math::clapack_lauum<nm::RubyObject, false>

*/
  };

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(rb_eNotImpError, "does not yet work for non-BLAS dtypes (needs herk, syrk, trmm)");
  } else {
    // Call either our version of lauum or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), blas_uplo_sym(uplo), FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda));
  }

  return a;
}


/* Call any of the clapack_xgetrf functions as directly as possible.
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
static VALUE nm_clapack_getrf(VALUE self, VALUE order, VALUE m, VALUE n, VALUE a, VALUE lda) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const int m, const int n, void* a, const int lda, int* ipiv) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::clapack_getrf<float>,
      nm::math::clapack_getrf<double>,
#ifdef HAVE_CLAPACK_H
      clapack_cgetrf, clapack_zgetrf, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::clapack_getrf<nm::Complex64>,
      nm::math::clapack_getrf<nm::Complex128>,
#endif
      nm::math::clapack_getrf<nm::Rational32>,
      nm::math::clapack_getrf<nm::Rational64>,
      nm::math::clapack_getrf<nm::Rational128>,
      nm::math::clapack_getrf<nm::RubyObject>
  };

  int M = FIX2INT(m),
      N = FIX2INT(n);

  // Allocate the pivot index array, which is of size MIN(M, N).
  size_t ipiv_size = std::min(M,N);
  int* ipiv = ALLOCA_N(int, ipiv_size);

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    // Call either our version of getrf or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), M, N, NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), ipiv);
  }

  // Result will be stored in a. We return ipiv as an array.
  VALUE ipiv_array = rb_ary_new2(ipiv_size);
  for (size_t i = 0; i < ipiv_size; ++i) {
    rb_ary_store(ipiv_array, i, INT2FIX(ipiv[i]));
  }

  return ipiv_array;
}


/* Call any of the clapack_xpotrf functions as directly as possible.
 *
 * You probably don't want to call this function. Instead, why don't you try clapack_potrf, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 *
 * Returns an array giving the pivot indices (normally these are argument #5).
 */
static VALUE nm_clapack_potrf(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda) {
#ifndef HAVE_CLAPACK_H
  rb_raise(rb_eNotImpError, "potrf currently requires LAPACK");
#endif

  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const enum CBLAS_UPLO, const int n, void* a, const int lda) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::clapack_potrf<float>,
      nm::math::clapack_potrf<double>,
#ifdef HAVE_CLAPACK_H
      clapack_cpotrf, clapack_zpotrf, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::clapack_potrf<nm::Complex64>,
      nm::math::clapack_potrf<nm::Complex128>,
#endif
      NULL, NULL, NULL, NULL /*
      nm::math::clapack_potrf<nm::Rational32>,
      nm::math::clapack_potrf<nm::Rational64>,
      nm::math::clapack_potrf<nm::Rational128>,
      nm::math::clapack_potrf<nm::RubyObject> */
  };

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
    // FIXME: Once BLAS dtypes are implemented, replace error above with the error below.
    //rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    // Call either our version of potrf or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), blas_uplo_sym(uplo), FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda));
  }

  return a;
}


/*
 * Call any of the clapack_xgetrs functions as directly as possible.
 */
static VALUE nm_clapack_getrs(VALUE self, VALUE order, VALUE trans, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE ipiv, VALUE b, VALUE ldb) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans, const int N,
                                       const int NRHS, const void* A, const int lda, const int* ipiv, void* B,
                                       const int ldb) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::clapack_getrs<float>,
      nm::math::clapack_getrs<double>,
#ifdef HAVE_CLAPACK_H
      clapack_cgetrs, clapack_zgetrs, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::clapack_getrs<nm::Complex64>,
      nm::math::clapack_getrs<nm::Complex128>,
#endif
      nm::math::clapack_getrs<nm::Rational32>,
      nm::math::clapack_getrs<nm::Rational64>,
      nm::math::clapack_getrs<nm::Rational128>,
      nm::math::clapack_getrs<nm::RubyObject>
  };

  // Allocate the C version of the pivot index array
  // TODO: Allow for an NVector here also, maybe?
  int* ipiv_;
  if (TYPE(ipiv) != T_ARRAY) {
    rb_raise(rb_eArgError, "ipiv must be of type Array");
  } else {
    ipiv_ = ALLOCA_N(int, RARRAY_LEN(ipiv));
    for (int index = 0; index < RARRAY_LEN(ipiv); ++index) {
      ipiv_[index] = FIX2INT( RARRAY_PTR(ipiv)[index] );
    }
  }

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {

    // Call either our version of getrs or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), blas_transpose_sym(trans), FIX2INT(n), FIX2INT(nrhs), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda),
                        ipiv_, NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));
  }

  // b is both returned and modified directly in the argument list.
  return b;
}


/*
 * Call any of the clapack_xpotrs functions as directly as possible.
 */
static VALUE nm_clapack_potrs(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE nrhs, VALUE a, VALUE lda, VALUE b, VALUE ldb) {
  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const int N,
                                       const int NRHS, const void* A, const int lda, void* B, const int ldb) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::clapack_potrs<float,false>,
      nm::math::clapack_potrs<double,false>,
#ifdef HAVE_CLAPACK_H
      clapack_cpotrs, clapack_zpotrs, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::clapack_potrs<nm::Complex64,true>,
      nm::math::clapack_potrs<nm::Complex128,true>,
#endif
      nm::math::clapack_potrs<nm::Rational32,false>,
      nm::math::clapack_potrs<nm::Rational64,false>,
      nm::math::clapack_potrs<nm::Rational128,false>,
      nm::math::clapack_potrs<nm::RubyObject,false>
  };


  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {

    // Call either our version of potrs or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), blas_uplo_sym(uplo), FIX2INT(n), FIX2INT(nrhs), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda),
                        NM_STORAGE_DENSE(b)->elements, FIX2INT(ldb));
  }

  // b is both returned and modified directly in the argument list.
  return b;
}


/* Call any of the clapack_xgetri functions as directly as possible.
 *
 * You probably don't want to call this function. Instead, why don't you try clapack_getri, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 *
 * Returns an array giving the pivot indices (normally these are argument #5).
 */
static VALUE nm_clapack_getri(VALUE self, VALUE order, VALUE n, VALUE a, VALUE lda, VALUE ipiv) {
#ifndef HAVE_CLAPACK_H
  rb_raise(rb_eNotImpError, "getri currently requires LAPACK");
#endif

  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const int n, void* a, const int lda, const int* ipiv) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::clapack_getri<float>,
      nm::math::clapack_getri<double>,
#ifdef HAVE_CLAPACK_H
      clapack_cgetri, clapack_zgetri, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::clapack_getri<nm::Complex64>,
      nm::math::clapack_getri<nm::Complex128>,
#endif
      NULL, NULL, NULL, NULL /*
      nm::math::clapack_getri<nm::Rational32>,
      nm::math::clapack_getri<nm::Rational64>,
      nm::math::clapack_getri<nm::Rational128>,
      nm::math::clapack_getri<nm::RubyObject> */
  };

  // Allocate the C version of the pivot index array
  // TODO: Allow for an NVector here also, maybe?
  int* ipiv_;
  if (TYPE(ipiv) != T_ARRAY) {
    rb_raise(rb_eArgError, "ipiv must be of type Array");
  } else {
    ipiv_ = ALLOCA_N(int, RARRAY_LEN(ipiv));
    for (int index = 0; index < RARRAY_LEN(ipiv); ++index) {
      ipiv_[index] = FIX2INT( RARRAY_PTR(ipiv)[index] );
    }
  }

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
    // FIXME: Once BLAS dtypes are implemented, replace error above with the error below.
    //rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    // Call either our version of getri or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), ipiv_);
  }

  return a;
}


/* Call any of the clapack_xpotri functions as directly as possible.
 *
 * You probably don't want to call this function. Instead, why don't you try clapack_potri, which is more flexible
 * with its arguments?
 *
 * This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
 * handling, so you can easily crash Ruby!
 *
 * Returns an array giving the pivot indices (normally these are argument #5).
 */
static VALUE nm_clapack_potri(VALUE self, VALUE order, VALUE uplo, VALUE n, VALUE a, VALUE lda) {
#ifndef HAVE_CLAPACK_H
  rb_raise(rb_eNotImpError, "getri currently requires LAPACK");
#endif

  static int (*ttable[nm::NUM_DTYPES])(const enum CBLAS_ORDER, const enum CBLAS_UPLO, const int n, void* a, const int lda) = {
      NULL, NULL, NULL, NULL, NULL, // integers not allowed due to division
      nm::math::clapack_potri<float>,
      nm::math::clapack_potri<double>,
#ifdef HAVE_CLAPACK_H
      clapack_cpotri, clapack_zpotri, // call directly, same function signature!
#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::clapack_potri<nm::Complex64>,
      nm::math::clapack_potri<nm::Complex128>,
#endif
      NULL, NULL, NULL, NULL /*
      nm::math::clapack_getri<nm::Rational32>,
      nm::math::clapack_getri<nm::Rational64>,
      nm::math::clapack_getri<nm::Rational128>,
      nm::math::clapack_getri<nm::RubyObject> */
  };

  if (!ttable[NM_DTYPE(a)]) {
    rb_raise(rb_eNotImpError, "this operation not yet implemented for non-BLAS dtypes");
    // FIXME: Once BLAS dtypes are implemented, replace error above with the error below.
    //rb_raise(nm_eDataTypeError, "this matrix operation undefined for integer matrices");
  } else {
    // Call either our version of getri or the LAPACK version.
    ttable[NM_DTYPE(a)](blas_order_sym(order), blas_uplo_sym(uplo), FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda));
  }

  return a;
}


/*
 * Call any of the clapack_xlaswp functions as directly as possible.
 *
 * Note that LAPACK's xlaswp functions accept a column-order matrix, but NMatrix uses row-order. Thus, n should be the
 * number of rows and lda should be the number of columns, no matter what it says in the documentation for dlaswp.f.
 */
static VALUE nm_clapack_laswp(VALUE self, VALUE n, VALUE a, VALUE lda, VALUE k1, VALUE k2, VALUE ipiv, VALUE incx) {
  static void (*ttable[nm::NUM_DTYPES])(const int n, void* a, const int lda, const int k1, const int k2, const int* ipiv, const int incx) = {
      nm::math::clapack_laswp<uint8_t>,
      nm::math::clapack_laswp<int8_t>,
      nm::math::clapack_laswp<int16_t>,
      nm::math::clapack_laswp<int32_t>,
      nm::math::clapack_laswp<int64_t>,
      nm::math::clapack_laswp<float>,
      nm::math::clapack_laswp<double>,
//#ifdef HAVE_CLAPACK_H // laswp doesn't actually exist in clapack.h!
//      clapack_claswp, clapack_zlaswp, // call directly, same function signature!
//#else // Especially important for Mac OS, which doesn't seem to include the ATLAS clapack interface.
      nm::math::clapack_laswp<nm::Complex64>,
      nm::math::clapack_laswp<nm::Complex128>,
//#endif
      nm::math::clapack_laswp<nm::Rational32>,
      nm::math::clapack_laswp<nm::Rational64>,
      nm::math::clapack_laswp<nm::Rational128>,
      nm::math::clapack_laswp<nm::RubyObject>
  };

  // Allocate the C version of the pivot index array
  // TODO: Allow for an NVector here also, maybe?
  int* ipiv_;
  if (TYPE(ipiv) != T_ARRAY) {
    rb_raise(rb_eArgError, "ipiv must be of type Array");
  } else {
    ipiv_ = ALLOCA_N(int, RARRAY_LEN(ipiv));
    for (int index = 0; index < RARRAY_LEN(ipiv); ++index) {
      ipiv_[index] = FIX2INT( RARRAY_PTR(ipiv)[index] );
    }
  }

  // Call either our version of laswp or the LAPACK version.
  ttable[NM_DTYPE(a)](FIX2INT(n), NM_STORAGE_DENSE(a)->elements, FIX2INT(lda), FIX2INT(k1), FIX2INT(k2), ipiv_, FIX2INT(incx));

  // a is both returned and modified directly in the argument list.
  return a;
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
