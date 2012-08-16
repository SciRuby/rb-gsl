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
#include <limits> // std::numeric_limits

/*
 * Project Includes
 */
#include "data/data.h"

/*
 * Macros
 */
#ifndef NM_MAX
# define NM_MAX(a,b) (((a)>(b))?(a):(b))
# define NM_MIN(a,b) (((a)>(b))?(b):(a))
#endif

/*
 * Types
 */

// These allow an increase in precision for intermediate values of gemm and gemv.
// See also: http://stackoverflow.com/questions/11873694/how-does-one-increase-precision-in-c-templates-in-a-template-typename-dependen
template <typename DType> struct LongDType;
template <> struct LongDType<uint8_t> { typedef int16_t type; };
template <> struct LongDType<int8_t> { typedef int16_t type; };
template <> struct LongDType<int16_t> { typedef int32_t type; };
template <> struct LongDType<int32_t> { typedef int64_t type; };
template <> struct LongDType<int64_t> { typedef int64_t type; };
template <> struct LongDType<float> { typedef double type; };
template <> struct LongDType<double> { typedef double type; };
template <> struct LongDType<Complex64> { typedef Complex128 type; };
template <> struct LongDType<Complex128> { typedef Complex128 type; };
template <> struct LongDType<Rational32> { typedef Rational128 type; };
template <> struct LongDType<Rational64> { typedef Rational128 type; };
template <> struct LongDType<Rational128> { typedef Rational128 type; };
template <> struct LongDType<RubyObject> { typedef RubyObject type; };


/*
 * Data
 */

/*
 * Functions
 */
void Init_blas(void);

/*
 * GEneral Matrix Multiplication: based on dgemm.f from Netlib.
 *
 * This is an extremely inefficient algorithm. Recommend using ATLAS' version instead.
 *
 * Template parameters: LT -- long version of type T. Type T is the matrix dtype.
 */
template <typename DType>
inline bool gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const DType* alpha, const DType* A, const int lda, const DType* B, const int ldb, const DType* beta, DType* C, const int ldc) {
  int num_rows_a, /*num_cols_a,*/ num_rows_b; // nrowa, ncola, nrowb

  typename LongDType<DType>::type temp;

  // %%= if [:rational,:complex,:value].include?(dtype.type); "#{dtype.long_dtype.sizeof} temp1, temp2;"; end%%
  int i, j, l;

  if (TransA == CblasNoTrans) num_rows_a = M;
  else                        num_rows_a = K;

  if (TransB == CblasNoTrans) num_rows_b = K;
  else                        num_rows_b = N;

  // Test the input parameters
  if (TransA < 111 || TransA > 113) {
    rb_raise(rb_eArgError, "GEMM: TransA must be CblasNoTrans, CblasTrans, or CblasConjTrans");
    return false;
  } else if (TransB < 111 || TransB > 113) {
    rb_raise(rb_eArgError, "GEMM: TransB must be CblasNoTrans, CblasTrans, or CblasConjTrans");
    return false;
  } else if (M < 0) {
    rb_raise(rb_eArgError, "GEMM: Expected M >= 0");
    return false;
  } else if (N < 0) {
    rb_raise(rb_eArgError, "GEMM: Expected N >= 0");
    return false;
  } else if (K < 0) {
    rb_raise(rb_eArgError, "GEMM: Expected K >= 0");
    return false;
  } else if (lda < NM_MAX(1, num_rows_a)) {
    fprintf(stderr, "GEMM: num_rows_a = %d; got lda=%d\n", num_rows_a, lda);
    rb_raise(rb_eArgError, "GEMM: Expected lda >= max(1, num_rows_a)");
    return false;
  } else if (ldb < NM_MAX(1, num_rows_b)) {
    fprintf(stderr, "GEMM: num_rows_b = %d; got ldb=%d\n", num_rows_b, ldb);
    rb_raise(rb_eArgError, "GEMM: Expected ldb >= max(1, num_rows_b)");
    return false;
  } else if (ldc < NM_MAX(1,M)) {
    fprintf(stderr, "GEMM: M=%d; got ldc=%d\n", M, ldc);
    rb_raise(rb_eArgError, "GEMM: Expected ldc >= max(1,M)");
    return false;
  }

  // Quick return if possible
  if (!M or !N or ((*alpha == 0 or !K) and *beta == 1)) return true;

  // For alpha = 0
  if (*alpha == 0) {
    if (*beta == 0) {
      for (j = 0; j < N; ++j)
        for (i = 0; i < M; ++i) {
          C[i+j*ldc] = 0;
        }
    } else {
      for (j = 0; j < N; ++j)
        for (i = 0; i < M; ++i) {
          C[i+j*ldc] *= *beta;
        }
    }
    return false;
  }

  // Start the operations
  if (TransB == CblasNoTrans) {
    if (TransA == CblasNoTrans) {
      // C = alpha*A*B+beta*C
      for (j = 0; j < N; ++j) {
        if (*beta == 0) {
          for (i = 0; i < M; ++i) {
            C[i+j*ldc] = 0;
          }
        } else if (*beta != 1) {
          for (i = 0; i < M; ++i) {
            C[i+j*ldc] *= *beta;
          }
        }

        for (l = 0; l < K; ++l) {
          if (B[l+j*ldb] != 0) {
            temp = *alpha * B[l+j*ldb];
            for (i = 0; i < M; ++i) {
              C[i+j*ldc] += A[i+l*lda] * temp;
            }
          }
        }
      }

    } else {

      // C = alpha*A**DType*B + beta*C
      for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {
          temp = 0;
          for (l = 0; l < K; ++l) {
            temp += A[l+i*lda] * B[l+j*ldb];
          }

          if (*beta == 0) {
            C[i+j*ldc] = *alpha*temp;
          } else {
            C[i+j*ldc] = *alpha*temp + *beta*C[i+j*ldc];
          }
        }
      }

    }

  } else if (TransA == CblasNoTrans) {

    // C = alpha*A*B**T + beta*C
    for (j = 0; j < N; ++j) {
      if (*beta == 0) {
        for (i = 0; i < M; ++i) {
          C[i+j*ldc] = 0;
        }
      } else if (*beta != 1) {
        for (i = 0; i < M; ++i) {
          C[i+j*ldc] *= *beta;
        }
      }

      for (l = 0; l < K; ++l) {
        if (B[j+l*ldb] != 0) {
          temp = *alpha * B[j+l*ldb];
          for (i = 0; i < M; ++i) {
            C[i+j*ldc] += A[i+l*lda] * temp;
          }
        }
      }

    }

  } else {

    // C = alpha*A**DType*B**T + beta*C
    for (j = 0; j < N; ++j) {
      for (i = 0; i < M; ++i) {
        temp = 0;
        for (l = 0; l < K; ++l) {
          temp += A[l+i*lda] * B[j+l*ldb];
        }

        if (*beta == 0) {
          C[i+j*ldc] = *alpha*temp;
        } else {
          C[i+j*ldc] = *alpha*temp + *beta*C[i+j*ldc];
        }
      }
    }

  }

  return true;
}


template <>
inline bool gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const float* alpha, const float* A, const int lda, const float* B, const int ldb, const float* beta, float* C, const int ldc) {
  cblas_sgemm(CblasRowMajor, TransA, TransB, N, M, K, *alpha, B, ldb, A, lda, *beta, C, ldc);
  return true;
}

template <>
inline bool gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const double* alpha, const double* A, const int lda, const double* B, const int ldb, const double* beta, double* C, const int ldc) {
  cblas_dgemm(CblasRowMajor, TransA, TransB, N, M, K, *alpha, B, ldb, A, lda, *beta, C, ldc);
  return true;
}

template <>
inline bool gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const Complex64* alpha, const Complex64* A, const int lda, const Complex64* B, const int ldb, const Complex64* beta, Complex64* C, const int ldc) {
  cblas_cgemm(CblasRowMajor, TransA, TransB, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
  return true;
}

template <>
inline bool gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
          const Complex128* alpha, const Complex128* A, const int lda, const Complex128* B, const int ldb, const Complex128* beta, Complex128* C, const int ldc) {
  cblas_zgemm(CblasRowMajor, TransA, TransB, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
  return true;
}


/*
 * GEneral Matrix-Vector multiplication: based on dgemv.f from Netlib.
 *
 * This is an extremely inefficient algorithm. Recommend using ATLAS' version instead.
 *
 * Template parameters: LT -- long version of type T. Type T is the matrix dtype.
 */
template <typename DType>
inline bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const DType* alpha, const DType* A, const int lda,
          const DType* X, const int incX, const DType* beta, DType* Y, const int incY) {
  int lenX, lenY, i, j;
  int kx, ky, iy, jx, jy, ix;

  typename LongDType<DType>::type temp;

  // Test the input parameters
  if (Trans < 111 || Trans > 113) {
    rb_raise(rb_eArgError, "GEMV: TransA must be CblasNoTrans, CblasTrans, or CblasConjTrans");
    return false;
  } else if (lda < NM_MAX(1, N)) {
    fprintf(stderr, "GEMV: N = %d; got lda=%d", N, lda);
    rb_raise(rb_eArgError, "GEMV: Expected lda >= max(1, N)");
    return false;
  } else if (incX == 0) {
    rb_raise(rb_eArgError, "GEMV: Expected incX != 0\n");
    return false;
  } else if (incY == 0) {
    rb_raise(rb_eArgError, "GEMV: Expected incY != 0\n");
    return false;
  }

  // Quick return if possible
  if (!M or !N or (*alpha == 0 and *beta == 1)) return true;

  if (Trans == CblasNoTrans) {
    lenX = N;
    lenY = M;
  } else {
    lenX = M;
    lenY = N;
  }

  if (incX > 0) kx = 0;
  else          kx = (lenX - 1) * -incX;

  if (incY > 0) ky = 0;
  else          ky =  (lenY - 1) * -incY;

  // Start the operations. In this version, the elements of A are accessed sequentially with one pass through A.
  if (*beta != 1) {
    if (incY == 1) {
      if (*beta == 0) {
        for (i = 0; i < lenY; ++i) {
          Y[i] = 0;
        }
      } else {
        for (i = 0; i < lenY; ++i) {
          Y[i] *= *beta;
        }
      }
    } else {
      iy = ky;
      if (*beta == 0) {
        for (i = 0; i < lenY; ++i) {
          Y[iy] = 0;
          iy += incY;
        }
      } else {
        for (i = 0; i < lenY; ++i) {
          Y[iy] *= *beta;
          iy += incY;
        }
      }
    }
  }

  if (*alpha == 0) return false;

  if (Trans == CblasNoTrans) {

    // Form  y := alpha*A*x + y.
    jx = kx;
    if (incY == 1) {
      for (j = 0; j < N; ++j) {
        if (X[jx] != 0) {
          temp = *alpha * X[jx];
          for (i = 0; i < M; ++i) {
            Y[i] += A[j+i*lda] * temp;
          }
        }
        jx += incX;
      }
    } else {
      for (j = 0; j < N; ++j) {
        if (X[jx] != 0) {
          temp = *alpha * X[jx];
          iy = ky;
          for (i = 0; i < M; ++i) {
            Y[iy] += A[j+i*lda] * temp;
            iy += incY;
          }
        }
        jx += incX;
      }
    }

  } else { // TODO: Check that indices are correct! They're switched for C.

    // Form  y := alpha*A**DType*x + y.
    jy = ky;

    if (incX == 1) {
      for (j = 0; j < N; ++j) {
        temp = 0;
        for (i = 0; i < M; ++i) {
          temp += A[j+i*lda]*X[j];
        }
        Y[jy] += *alpha * temp;
        jy += incY;
      }
    } else {
      for (j = 0; j < N; ++j) {
        temp = 0;
        ix = kx;
        for (i = 0; i < M; ++i) {
          temp += A[j+i*lda] * X[ix];
          ix += incX;
        }

        Y[jy] += *alpha * temp;
        jy += incY;
      }
    }
  }

  return true;
}  // end of GEMV

template <>
inline bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const float* alpha, const float* A, const int lda,
          const float* X, const int incX, const float* beta, float* Y, const int incY) {
  cblas_sgemv(CblasRowMajor, Trans, M, N, *alpha, A, lda, X, incX, *beta, Y, incY);
  return true;
}

template <>
inline bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const double* alpha, const double* A, const int lda,
          const double* X, const int incX, const double* beta, double* Y, const int incY) {
  cblas_dgemv(CblasRowMajor, Trans, M, N, *alpha, A, lda, X, incX, *beta, Y, incY);
  return true;
}

template <>
inline bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const Complex64* alpha, const Complex64* A, const int lda,
          const Complex64* X, const int incX, const Complex64* beta, Complex64* Y, const int incY) {
  cblas_cgemv(CblasRowMajor, Trans, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  return true;
}

template <>
inline bool gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const Complex128* alpha, const Complex128* A, const int lda,
          const Complex128* X, const int incX, const Complex128* beta, Complex128* Y, const int incY) {
  cblas_zgemv(CblasRowMajor, Trans, M, N, alpha, A, lda, X, incX, beta, Y, incY);
  return true;
}


// Yale: numeric matrix multiply c=a*b
template <typename DType, typename IType>
inline void numbmm(const unsigned int n, const unsigned int m, const IType* ia, const IType* ja, const DType* a, const bool diaga,
            const IType* ib, const IType* jb, const DType* b, const bool diagb, IType* ic, IType* jc, DType* c, const bool diagc) {
  IType next[m];
  DType sums[m];

  DType v;

  IType head, length, temp, ndnz = 0;
  IType jj_start, jj_end, kk_start, kk_end;
  IType i, j, k, kk, jj;
  IType minmn = NM_MIN(m,n);

  for (i = 0; i < m; ++i) { // initialize scratch arrays
    next[i] = std::numeric_limits<IType>::max();
    sums[i] = 0;
  }

  for (i = 0; i < n; ++i) { // walk down the rows
    head = std::numeric_limits<IType>::max()-1; // head gets assigned as whichever column of B's row j we last visited
    length = 0;

    jj_start = ia[i];
    jj_end   = ia[i+1];

    for (jj = jj_start; jj <= jj_end; ++jj) { // walk through entries in each row

      if (jj == jj_end) { // if we're in the last entry for this row:
        if (!diaga || i >= minmn) continue;
        j   = i;      // if it's a new Yale matrix, and last entry, get the diagonal position (j) and entry (ajj)
        v   = a[i];
      } else {
        j   = ja[jj]; // if it's not the last entry for this row, get the column (j) and entry (ajj)
        v   = a[jj];
      }

      kk_start = ib[j];   // Find the first entry of row j of matrix B
      kk_end   = ib[j+1];
      for (kk = kk_start; kk <= kk_end; ++kk) {

        if (kk == kk_end) { // Get the column id for that entry
          if (!diagb || j >= minmn) continue;
          k  = j;
          sums[k] += v*b[k];
        } else {
          k  = jb[kk];
          sums[k] += v*b[kk];
        }

        if (next[k] == std::numeric_limits<IType>::max()) {
          next[k] = head;
          head    = k;
          ++length;
        }
      }
    }

    for (jj = 0; jj < length; ++jj) {
      if (sums[head] != 0) {
        if (diagc && head == i) {
          c[head] = sums[head];
        } else {
          jc[n+1+ndnz] = head;
          c[n+1+ndnz]  = sums[head];
          ++ndnz;
        }
      }

      temp = head;
      head = next[head];

      next[temp] = std::numeric_limits<IType>::max();
      sums[temp] = 0;
    }

    ic[i+1] = n+1+ndnz;
  }
} /* numbmm_ */



// Yale: Symbolic matrix multiply c=a*b
template <typename IType>
inline void symbmm(const unsigned int n, const unsigned int m, const IType* ia, const IType* ja, const bool diaga,
            const IType* ib, const IType* jb, const bool diagb, IType* ic, const bool diagc) {
  IType mask[m];
  IType j, k, ndnz = n; /* Local variables */


  for (j = 0; j < m; ++j)
    mask[j] = std::numeric_limits<IType>::max();

  if (diagc)  ic[0] = n+1;
  else        ic[0] = 0;

  IType minmn = NM_MIN(m,n);

  for (IType i = 0; i < n; ++i) { // MAIN LOOP: through rows

    for (IType jj = ia[i]; jj <= ia[i+1]; ++jj) { // merge row lists, walking through columns in each row

      // j <- column index given by JA[jj], or handle diagonal.
      if (jj == ia[i+1]) { // Don't really do it the last time -- just handle diagonals in a new yale matrix.
        if (!diaga || i >= minmn) continue;
        j = i;
      } else j = ja[jj];

      for (IType kk = ib[j]; kk <= ib[j+1]; ++kk) { // Now walk through columns of row J in matrix B.
        if (kk == ib[j+1]) {
          if (!diagb || j >= minmn) continue;
          k = j;
        } else k = jb[kk];

        if (mask[k] != i) {
          mask[k] = i;
          ++ndnz;
        }
      }
    }

    if (diagc && !mask[i]) --ndnz;

    ic[i+1] = ndnz;
  }
} /* symbmm_ */


//TODO: More efficient sorting algorithm than selection sort would be nice, probably.
// Remember, we're dealing with unique keys, which simplifies things.
// Doesn't have to be in-place, since we probably just multiplied and that wasn't in-place.
template <typename DType, typename IType>
inline void smmp_sort_columns(const unsigned int n, const IType* ia, IType* ja, DType* a) {
  IType i, jj, jj_start, jj_end, min, min_jj;
  DType temp_val;

  for (i = 0; i < n; ++i) {
    // No need to sort if there are 0 or 1 entries in the row
    if (ia[i+1] - ia[i] < 2) continue;

    jj_end = ia[i+1];
    for (jj_start = ia[i]; jj_start < jj_end; ++jj_start) {

      // If the previous min is just current-1, this key/value pair is already in sorted order.
      // This follows from the unique condition on our column keys.
      if (jj_start > ia[i] && min+1 == ja[jj_start]) {
        min    = ja[jj_start];
        continue;
      }

      // find the minimum key (column index) between jj_start and jj_end
      min    = ja[jj_start];
      min_jj = jj_start;
      for (jj = jj_start+1; jj < jj_end; ++jj) {
        if (ja[jj] < min) {
          min_jj = jj;
          min    = ja[jj];
        }
      }

      // if min is already first, skip this iteration
      if (min_jj == jj_start) continue;

      for (jj = jj_start; jj < jj_end; ++jj) {
        // swap minimum key/value pair with key/value pair in the first position.
        if (min_jj != jj) {
          // min already = ja[min_jj], so use this as temp_key
          temp_val = a[min_jj];

          ja[min_jj] = ja[jj];
          a[min_jj] = a[jj];

          ja[jj] = min;
          a[jj] = temp_val;
        }
      }
    }
  }
}


/*
 * Transposes a generic Yale matrix (old or new). Specify new by setting diaga = true.
 *
 * Based on transp from SMMP (same as symbmm and numbmm).
 *
 * This is not named in the same way as most yale_storage functions because it does not act on a YALE_STORAGE
 * object.
 */
template <typename DType, typename IType>
void transpose_yale_template(const size_t n, const size_t m, const void* ia_, const void* ja_, const void* a_,
                             const bool diaga, void* ib_, void* jb_, void* b_, const bool move)
{
  const IType *ia = reinterpret_cast<const IType*>(ia_),
              *ja = reinterpret_cast<const IType*>(ja_);
  const DType *a  = reinterpret_cast<const DType*>(a_);

  IType *ib = reinterpret_cast<IType*>(ib_),
        *jb = reinterpret_cast<IType*>(jb_);
  DType *b  = reinterpret_cast<DType*>(b_);



  size_t index;

  // Clear B
  for (size_t i = 0; i < m+1; ++i) ib[i] = 0;

  if (move)
    for (size_t i = 0; i < m+1; ++i) b[i] = 0;

  if (diaga) ib[0] = m + 1;
  else       ib[0] = 0;

  /* count indices for each column */

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = ia[i]; j < ia[i+1]; ++j) {
      ++(ib[ja[j]+1]);
    }
  }

  for (size_t i = 0; i < m; ++i) {
    ib[i+1] = ib[i] + ib[i+1];
  }

  /* now make jb */

  for (size_t i = 0; i < n; ++i) {

    for (size_t j = ia[i]; j < ia[i+1]; ++j) {
      index = ja[j];
      jb[ib[index]] = i;

      if (move)
        b[ib[index]] = a[j];

      ++(ib[index]);
    }
  }

  /* now fixup ib */

  for (size_t i = m; i >= 1; --i) {
    ib[i] = ib[i-1];
  }


  if (diaga) {
    if (move) {
      size_t j = NM_MIN(n,m);

      for (size_t i = 0; i < j; ++i) {
        b[i] = a[i];
      }
    }
    ib[0] = m + 1;

  } else {
    ib[0] = 0;
  }
}



void det_exact(const int M, const void* elements, const int lda, dtype_t dtype, void* result);
void transpose_generic(const size_t M, const size_t N, const void* A, const int lda, void* B, const int ldb, size_t element_size);

#endif // MATH_H
