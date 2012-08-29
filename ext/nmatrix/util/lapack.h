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
// == lapack.h
//
// Templated versions of LAPACK functions, in C++.

#ifndef LAPACK_H
  #define LAPACK_H

#include <cmath> // std::round

#include "math.h"

namespace nm { namespace math {



template <typename DType>
extern bool gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K, const DType* alpha,
                 const DType* A, const int lda, const DType* B, const int ldb,
                 const DType* beta, DType* C, const int ldc);

namespace lapack {



/*  -- LAPACK auxiliary routine (version 3.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2006 */

/*     .. Scalar Arguments .. */

/*  Purpose */
/*  ======= */

/*       This program sets problem and machine dependent parameters */
/*       useful for xHSEQR and its subroutines. It is called whenever */
/*       ILAENV is called with 12 <= ISPEC <= 16 */

/*  Arguments */
/*  ========= */

/*       ISPEC  (input) int scalar */
/*              ISPEC specifies which tunable parameter IPARMQ should */
/*              return. */

/*              ISPEC=12: (INMIN)  Matrices of order nmin or less */
/*                        are sent directly to xLAHQR, the implicit */
/*                        double shift QR algorithm.  NMIN must be */
/*                        at least 11. */

/*              ISPEC=13: (INWIN)  Size of the deflation window. */
/*                        This is best set greater than or equal to */
/*                        the number of simultaneous shifts NS. */
/*                        Larger matrices benefit from larger deflation */
/*                        windows. */

/*              ISPEC=14: (INIBL) Determines when to stop nibbling and */
/*                        invest in an (expensive) multi-shift QR sweep. */
/*                        If the aggressive early deflation subroutine */
/*                        finds LD converged eigenvalues from an order */
/*                        NW deflation window and LD.GT.(NW*NIBBLE)/100, */
/*                        then the next QR sweep is skipped and early */
/*                        deflation is applied immediately to the */
/*                        remaining active diagonal block.  Setting */
/*                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a */
/*                        multi-shift QR sweep whenever early deflation */
/*                        finds a converged eigenvalue.  Setting */
/*                        IPARMQ(ISPEC=14) greater than or equal to 100 */
/*                        prevents TTQRE from skipping a multi-shift */
/*                        QR sweep. */

/*              ISPEC=15: (NSHFTS) The number of simultaneous shifts in */
/*                        a multi-shift QR iteration. */

/*              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the */
/*                        following meanings. */
/*                        0:  During the multi-shift QR sweep, */
/*                            xLAQR5 does not accumulate reflections and */
/*                            does not use matrix-matrix multiply to */
/*                            update the far-from-diagonal matrix */
/*                            entries. */
/*                        1:  During the multi-shift QR sweep, */
/*                            xLAQR5 and/or xLAQRaccumulates reflections and uses */
/*                            matrix-matrix multiply to update the */
/*                            far-from-diagonal matrix entries. */
/*                        2:  During the multi-shift QR sweep. */
/*                            xLAQR5 accumulates reflections and takes */
/*                            advantage of 2-by-2 block structure during */
/*                            matrix-matrix multiplies. */
/*                        (If xTRMM is slower than xGEMM, then */
/*                        IPARMQ(ISPEC=16)=1 may be more efficient than */
/*                        IPARMQ(ISPEC=16)=2 despite the greater level of */
/*                        arithmetic work implied by the latter choice.) */

/*       NAME    (input) character string */
/*               Name of the calling subroutine */

/*       OPTS    (input) character string */
/*               This is a concatenation of the string arguments to */
/*               TTQRE. */

/*       N       (input) int scalar */
/*               N is the order of the Hessenberg matrix H. */

/*       ILO     (input) INTEGER */
/*       IHI     (input) INTEGER */
/*               It is assumed that H is already upper triangular */
/*               in rows and columns 1:ILO-1 and IHI+1:N. */

/*       LWORK   (input) int scalar */
/*               The amount of workspace available. */

/*  Further Details */
/*  =============== */

/*       Little is known about how best to choose these parameters. */
/*       It is possible to use different values of the parameters */
/*       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR. */

/*       It is probably best to choose different parameters for */
/*       different matrices and different parameters at different */
/*       times during the iteration, but this has not been */
/*       implemented --- yet. */


/*       The best choices of most of the parameters depend */
/*       in an ill-understood way on the relative execution */
/*       rate of xLAQR3 and xLAQR5 and on the nature of each */
/*       particular eigenvalue problem.  Experiment may be the */
/*       only practical way to determine which choices are most */
/*       effective. */

/*       Following is a list of default values supplied by IPARMQ. */
/*       These defaults may be adjusted in order to attain better */
/*       performance in any particular computational environment. */

/*       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point. */
/*                        Default: 75. (Must be at least 11.) */

/*       IPARMQ(ISPEC=13) Recommended deflation window size. */
/*                        This depends on ILO, IHI and NS, the */
/*                        number of simultaneous shifts returned */
/*                        by IPARMQ(ISPEC=15).  The default for */
/*                        (IHI-ILO+1).LE.500 is NS.  The default */
/*                        for (IHI-ILO+1).GT.500 is 3*NS/2. */

/*       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14. */

/*       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS. */
/*                        a multi-shift QR iteration. */

/*                        If IHI-ILO+1 is ... */

/*                        greater than      ...but less    ... the */
/*                        or equal to ...      than        default is */

/*                                0               30       NS =   2+ */
/*                               30               60       NS =   4+ */
/*                               60              150       NS =  10 */
/*                              150              590       NS =  ** */
/*                              590             3000       NS =  64 */
/*                             3000             6000       NS = 128 */
/*                             6000             infinity   NS = 256 */

/*                    (+)  By default matrices of this order are */
/*                         passed to the implicit double shift routine */
/*                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These */
/*                         values of NS are used only in case of a rare */
/*                         xLAHQR failure. */

/*                    (**) The asterisks (**) indicate an ad-hoc */
/*                         function increasing from 10 to 64. */

/*       IPARMQ(ISPEC=16) Select structured matrix multiply. */
/*                        (See ISPEC=16 above for details.) */
/*                        Default: 3. */

/*     ================================================================ */
inline int iparmq(int ispec, int ilo, int ihi) {

  const int INMIN = 12, INWIN = 13, INIBL = 14, ISHFTS = 15, IACC22 = 16;
  const int NMIN = 75, K22MIN = 14, KACMIN = 14, NIBBLE = 14, KNWSWP = 500;

  int ns = 2, nh = ihi - ilo + 1;

  if (ispec == ISHFTS || ispec == INWIN|| ispec == IACC22) {

    /*        ==== Set the number of simultaneous shifts ==== */
    if (nh >= 30)	  ns = 4;
    if (nh >= 60)   ns = 10;
    if (nh >= 150)  ns = std::max(10, (int)(nh / std::round(std::log((float) (ihi - ilo + 1)) / log(2.f)))); /* Computing MAX */
    if (nh >= 590)  ns = 64;
    if (nh >= 3000) ns = 128;
    if (nh >= 6000) ns = 256;
    ns = std::max(2,ns - ns % 2); /* Computing MAX */
  }

  if (ispec == INMIN) {
    /*        ===== Matrices of order smaller than NMIN get sent */
    /*        .     to xLAHQR, the classic double shift algorithm. */
    /*        .     This must be at least 11. ==== */
    return NMIN;

  } else if (ispec == INIBL) {

    /*        ==== INIBL: skip a multi-shift qr iteration and */
    /*        .    whenever aggressive early deflation finds */
    /*        .    at least (NIBBLE*(window size)/100) deflations. ==== */

    return NIBBLE;

  } else if (ispec == ISHFTS) {

    /*        ==== NSHFTS: The number of simultaneous shifts ===== */
    return ns;

  } else if (ispec == INWIN) {

    /*        ==== NW: deflation window size.  ==== */

    if (nh <= KNWSWP)  return ns;
    else               return ns * 3 / 2;

  } else if (ispec == 16) {

    /*        ==== IACC22: Whether to accumulate reflections */
    /*        .     before updating the far-from-diagonal elements */
    /*        .     and whether to use 2-by-2 block structure while */
    /*        .     doing it.  A small amount of work could be saved */
    /*        .     by making this choice dependent also upon the */
    /*        .     NH=IHI-ILO+1. */

    if (ns >= KACMIN) return 1;
    if (ns >= K22MIN) return 2;

  }

  return -1;
} /* iparmq_ */




/*  Purpose */
/*  ======= */

/*  DGER   performs the rank 1 operation */

/*     A := alpha*x*y**T + A, */

/*  where alpha is a scalar, x is an m element vector, y is an n element */
/*  vector and A is an m by n matrix. */

/*  Arguments */
/*  ========== */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - DOUBLE PRECISION. */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  X      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( m - 1 )*abs( INCX ) ). */
/*           Before entry, the incremented array X must contain the m */
/*           element vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  Y      - DOUBLE PRECISION array of dimension at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ). */
/*           Before entry, the incremented array Y must contain the n */
/*           element vector y. */
/*           Unchanged on exit. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */

/*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. On exit, A is */
/*           overwritten by the updated matrix. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, m ). */
/*           Unchanged on exit. */

/*  Further Details */
/*  =============== */

/*  Level 2 Blas routine. */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */

/*  ===================================================================== */

template <typename DType>
inline int ger(int m, int n, DType alpha, DType* x, int incx, DType* y, int incy, DType* a, int lda) {

  // FIXME: Call BLAS ger if available

  if (m < 0) {
    return 1;
  } else if (n < 0) {
    return 2;
  } else if (incx == 0) {
    return 5;
  } else if (incy == 0) {
    return 7;
  } else if (lda < std::max(1,m)) {
    return 9;
  }

  if (m == 0 || n == 0 || alpha == 0) return 0; /* Quick return if possible. */

  /*     Start the operations. In this version the elements of A are */
  /*     accessed sequentially with one pass through A. */

  // FIXME: These have been unrolled in a way that the compiler can handle. Collapse into a single case, or optimize
  // FIXME: in a more modern way.

  int jy = incy > 0 ? 0 :  -(n-1) * incy;

  if (incx == 1) {

	  for (size_t j = 0; j < n; ++j, jy += incy) {
	    if (y[jy] != 0) {
		    DType temp = alpha * y[jy];
		    for (size_t i = 0; i < m; ++i) {
		      a[i + j * lda] += x[i] * temp;
		    }
	    }
	  }

  } else {

    int kx = incx > 0 ? 0 : -(m-1) * incx;

	  for (size_t j = 0; j < n; ++j, jy += incy) {
	    if (y[jy] != 0) {
    		DType temp = alpha * y[jy];

    		for (size_t i = 0, ix = kx; i < m; ++i, ix += incx) {
          a[i + j * lda] += x[ix] * temp;
    		}
	    }
	  }

  }

  return 0;

/*     End of DGER  . */

} /* dger_ */


/*  Purpose */
/*  ======= */

/*     DSCAL scales a vector by a constant. */
/*     uses unrolled loops for increment equal to one. */

/*  Further Details */
/*  =============== */

/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */

/*  ===================================================================== */

template <typename DType>
inline void scal(const int n, const DType da, DType* dx,	const int incx) {

  // This used to have unrolled loops, like dswap. They were in the way.

  if (n <= 0 || incx <= 0) return;

  for (int i = 0; incx < 0 ? i > n*incx : i < n*incx; i += incx) {
    dx[i] = da * dx[i];
  }
} /* scal */


/*  Purpose */
/*  ======= */

/*     interchanges two vectors. */
/*     uses unrolled loops for increments equal one. */

/*  Further Details */
/*  =============== */

/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */

/*  ===================================================================== */
// Formerly dswap
template <typename DType>
inline void swap(int n, DType *dx, int incx, DType *dy, int incy) {

    /* Function Body */
    if (n <= 0) return;

    /*
     * The NETLIB version of dswap has loops manually unrolled, per commented code below.
     * This doesn't make sense with modern compilers, which know much more about arch-
     * itectures than we do. Combine that with our use of templates, and it's much more
     * efficient to let the compiler do the unrolling in most cases.
     */
/*
    if (incx == 1 && incy == 1) { // if both increments are 1


      m = n % 3;
      if (m) {
        for (size_t i = 0; i < m; ++i) { // If number is not divisible by three, swap just one or two singly.
          dtemp = dx[i];
          dx[i] = dy[i];
          dy[i] = dtemp;
        }
        if (n < 3) return;
      }

      for (i = m; i < n; i += 3) { // Why does it swap three at a time? -- John 8/27/12
        DType dtemp = dx[i];
        dx[i]       = dy[i];
        dy[i]       = dtemp;

        dtemp = dx[i + 1];
        dx[i + 1] = dy[i + 1];
        dy[i + 1] = dtemp;

        dtemp = dx[i + 2];
        dx[i + 2] = dy[i + 2];
        dy[i + 2] = dtemp;
      }

    } else { // when any increment is not 1
*/

  // For negative increments, start at the end of the array.
  int ix = incx < 0 ? (-n+1)*incx : 0,
      iy = incy < 0 ? (-n+1)*incy : 0;

  if (incx < 0) ix = (-n + 1) * incx;
  if (incy < 0) iy = (-n + 1) * incy;

  for (size_t i = 0; i < n; ++i, ix += incx, iy += incy) {
    DType dtemp = dx[ix];
    dx[ix]      = dy[iy];
    dy[iy]      = dtemp;
  }
  /*} */
  return;
} /* dswap */




/*  Purpose */
/*  ======= */

/*     IDAMAX finds the index of element having max. absolute value. */

/*  Further Details */
/*  =============== */

/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */

/*  ===================================================================== */

template <typename DType>
inline int idamax(int n, DType *dx, int incx) {

  /* Function Body */
  if (n < 1 || incx <= 0) return -1;
  if (n == 1)             return 0;

  DType dmax;
  size_t imax = 0;

  if (incx == 1) { // if incrementing by 1

    dmax = abs(dx[0]);

    for (size_t i = 1; i < n; ++i) {
      if (std::abs(dx[i]) > dmax) {
        imax = i;
        dmax = std::abs(dx[i]);
      }
    }

  } else { // if incrementing by more than 1

    dmax = std::abs(dx[0]);

    for (size_t i = 1, ix = incx; i < n; ++i, ix += incx) {
      if (std::abs(dx[ix]) > dmax) {
        imax = i;
        dmax = std::abs(dx[ix]);
      }
    }
  }
  return imax;
} /* idamax_ */




/* > \brief \b DLASWP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASWP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaswp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaswp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaswp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, K1, K2, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASWP performs a series of row interchanges on the matrix A. */
/* > One row interchange is initiated for each of rows K1 through K2 of A. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the matrix of column dimension N to which the row */
/* >          interchanges will be applied. */
/* >          On exit, the permuted matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* > \endverbatim */
/* > */
/* > \param[in] K1 */
/* > \verbatim */
/* >          K1 is INTEGER */
/* >          The first element of IPIV for which a row interchange will */
/* >          be done. */
/* > \endverbatim */
/* > */
/* > \param[in] K2 */
/* > \verbatim */
/* >          K2 is INTEGER */
/* >          The last element of IPIV for which a row interchange will */
/* >          be done. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (K2*abs(INCX)) */
/* >          The vector of pivot indices.  Only the elements in positions */
/* >          K1 through K2 of IPIV are accessed. */
/* >          IPIV(K) = L implies rows K and L are to be interchanged. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between successive values of IPIV.  If IPIV */
/* >          is negative, the pivots are applied in reverse order. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup DTypeOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Modified by */
/* >   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */
template <typename DType>
inline void laswp(int n, DType *a, int lda, int k1, int k2, int *ipiv, int incx) {
    /* System generated locals */
    int a_dim1, a_offset;


/*  -- LAPACK auxiliary routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Interchange row I with row IPIV(I) for each of rows K1 through K2. */

  /* Parameter adjustments */
  // a_offset = 1 + lda;
  // a -= a_offset;
  // --ipiv;

  int i1, i2, inc, ix0;

  /* Function Body */
  if (incx > 0) {
    ix0 = k1;
    i1 = k1;
    i2 = k2;
    inc = 1;
  } else if (incx < 0) {
    ix0 = (1 - k2) * incx + 1;
    i1 = k2;
    i2 = k1;
    inc = -1;
  } else {
    return;
  }

  int n32 = (n / 32) * 32;
  if (n32 != 0) {
    for (size_t j = 0; j < n32; j += 32) { // did change indices
      int ix = ix0;
      for (int i = i1; inc < 0 ? i >= i2 : i <= i2; i += inc) { // did not change indices
        int ip = ipiv[ix];
        if (ip != i) {
          for (size_t k = j; k < j+32; ++k) {
            DType temp = a[i + k * lda];
            a[i + k * lda] = a[ip + k * lda];
            a[ip + k * lda] = temp;
          }
        }
        ix += incx;
      }
    }
  }

  if (n32 != n) {
    ++n32;
    int ix = ix0;
    for (int i = i1; inc < 0 ? i >= i2 : i <= i2; i += inc) { // did not change indices
      int ip = ipiv[ix];
      if (ip != i) {

        for (size_t k = n32; k < n; ++k) { // did change indices
          DType temp = a[i + k * lda];
          a[i + k * lda] = a[ip + k * lda];
          a[ip + k * lda] = temp;
        }
      }
      ix += incx;
    }
  }

  return;

/*     End of LASWP */

} /* laswp */


/*  -- LAPACK auxiliary routine (version 3.3.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*  -- April 2011                                                      -- */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  IEEECK is called from the ILAENV to verify that Infinity and */
/*  possibly NaN arithmetic is safe (i.e. will not trap). */
// FIXME: Can we use std::numeric_limits::traps for this?

/*  Arguments */
/*  ========= */

/*  ISPEC   (input) INTEGER */
/*          Specifies whether to test just for inifinity arithmetic */
/*          or whether to test for infinity and NaN arithmetic. */
/*          = 0: Verify infinity arithmetic only. */
/*          = 1: Verify infinity and NaN arithmetic. */

/*  ZERO    (input) REAL */
/*          Must contain the value 0.0 */
/*          This is passed to prevent the compiler from optimizing */
/*          away this code. */

/*  ONE     (input) REAL */
/*          Must contain the value 1.0 */
/*          This is passed to prevent the compiler from optimizing */
/*          away this code. */

/*  RETURN VALUE:  INTEGER */
/*          = 0:  Arithmetic failed to produce the correct answers */
/*          = 1:  Arithmetic produced the correct answers */

/*  ===================================================================== */

/*
 * Note from John: This seems totally unnecessary in modern C++.
 * FIXME: Remove this after testing that on modern systems this always returns 1.
 */

inline int ieeeck(bool ispec) {

    float posinf = 1.0 / 0.0;
    if (posinf <= 1.0) return 0;

    float neginf = -1.0 / 0.0;
    if (neginf >= 0.0) return 0;

    float negzro = 1.0 / (neginf + 1.0);
    if (negzro != 0.0) return 0;

    neginf = 1.0 / negzro;
    if (neginf >= 0.0) return 0;

    float newzro = negzro + 0.0;
    if (newzro != 0.0) return 0;

    posinf = 1.0 / newzro;
    if (posinf <= 1.0) return 0;

    neginf *= posinf;
    if (neginf >= 0.0) return 0;

    posinf *= posinf;
    if (posinf <= 1.0) return 0;


/*     Return if we were only asked to check infinity arithmetic */

    if (!ispec) return 1;

    float nan1 = posinf + neginf;
    float nan2 = posinf / neginf;
    float nan3 = posinf / posinf;
    float nan4 = posinf * 0.0;
    float nan5 = neginf * negzro;
    float nan6 = nan5 * 0.0;

    if (nan1 == nan1) return 0;
    if (nan2 == nan2) return 0;
    if (nan3 == nan3) return 0;
    if (nan4 == nan4) return 0;
    if (nan5 == nan5) return 0;
    if (nan6 == nan6) return 0;

    return 1;
} /* ieeeck_ */




inline int ilaenv_block_size(int n2, int n4, const std::string& c2, const std::string& c3, const std::string& c4, bool sname, bool cname) {
    if (c2 == "GE") { //(s_cmp(c2, "GE", (size_t)2, (size_t)2) == 0) {
	    if (c3 == "TRF") { //if (s_cmp(c3, "TRF", (size_t)3, (size_t)3) == 0) {
	      if (sname)  return 64;
	      else        return 64;
	    } else if (c3 == "QRF" || c3 == "RQF" || c3 == "LQF" || c3 == "QLF") { //(s_cmp(c3, "QRF", (size_t)3, (size_t)3) == 0 || s_cmp(c3, "RQF", (size_t)3, (size_t)3) == 0 || s_cmp(c3, "LQF", (size_t) 3, (size_t)3) == 0 || s_cmp(c3, "QLF", (size_t)3, (size_t)3) == 0) {
	      if (sname)	return 32;
	      else    		return 32;
	    } else if (c3 == "HRD") {
	      if (sname)	return 32;
	      else    		return 32;
	    } else if (c3 == "BRD") {
	      if (sname)  return 32;
	      else    		return 32;
    	} else if (c3 == "TRI") {
	      if (sname) 	return 64;
        else    		return 64;
	    }
    } else if (c2 == "PO") {
	    if (c3 == "TRF") {
	      if (sname) 	return 64;
	      else    		return 64;
	    }
    } else if (c2 == "SY") {
	    if (c3 == "TRF") {
	      if (sname) 	return 64;
	      else    		return 64;
	    } else if (sname && c3 == "TRD") {
	      return 32;
	    } else if (sname && c3 == "GST") {
	      return 64;
	    }
    } else if (cname && c2 == "HE") {
	    if (c3 == "TRF")        return 64;
	    else if (c3 == "TRD")   return 32;
	    else if (c3 == "GST")   return 64;
    } else if (sname && c2 == "OR") {
	    if (c3.at(0) == 'G') {
        if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") return 32;
	    } else if (c3.at(0) == 'M') {
	      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") return 32;
	    }
    } else if (cname && c2 == "UN") {
	    if (c3.at(0) == 'G') {
        if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") return 32;
	    } else if (c3.at(0) == 'M') {
	      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") return 32;
	    }
    } else if (c2 == "GB") {
      if (c3 == "TRF") {
        if (sname) {
          if (n4 <= 64)     return 1;
          else              return 32;
	      } else {
          if (n4 <= 64)     return 1;
          else              return 32;
        }
      }
	  } else if (c2 == "PB") {
	    if (c3 == "TRF") {
	      if (sname) {
		      if (n2 <= 64)     return 1;
		      else      		    return 32;
	      } else {
		      if (n2 <= 64)     return 1;
		      else      		    return 32;
	      }
	    }
    } else if (c2 == "TR") {
	    if (c3 == "TRI") {
	      if (sname)  return 64;
	      else	      return 64;
	    }
    } else if (c2 == "LA") {
	    if (c3 == "UUM") {
	      if (sname)  return 64;
	      else        return 64;
	    }
    } else if (sname && c2 == "ST") {
	    if (c3 == "EBZ") return 1;
    }
    return 1;
}


inline int ilaenv_min_block_size(const std::string& c2, const std::string& c3, const std::string& c4, bool sname, bool cname) {

  if (c2 == "GE") {
    if (c3 == "QRF" || c3 == "RQF" || c3 == "LQF" || c3 == "QLF") {
      if (sname) {
        return 2;
      } else {
        return 2;
      }
    } else if (c3 == "HRD") {
      if (sname) {
        return 2;
      } else {
        return 2;
      }
    } else if (c3 == "BRD") {
      if (sname) {
        return 2;
      } else {
        return 2;
      }
    } else if (c3 == "TRI") {
      if (sname) {
        return 2;
      } else {
        return 2;
      }
    }
  } else if (c2 == "SY") {
    if (c3 == "TRF") {
      if (sname) {
        return 8;
      } else {
        return 8;
      }
    } else if (sname && c3 == "TRD") {
      return 2;
    }
  } else if (cname && c2 == "HE") {
    if (c3 == "TRD") {
      return 2;
    }
  } else if (sname && c2 == "OR") {
    if (c3.at(0) == 'G') {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
        return 2;
      }
    } else if (c3.at(0) == 'M') {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
        return 2;
      }
    }
  } else if (cname && c2 == "UN") {
    if (c3.at(0) == 'G') {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
        return 2;
      }
    } else if (c3.at(0) == 'M') {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
        return 2;
      }
    }
  }
  return 2;
}


inline int ilaenv_crossover_point(const std::string& c2, const std::string& c3, const std::string& c4, bool sname, bool cname) {
  if (c2 == "GE") {
    if (c3 == "QRF" || c3 == "RQF" || c3 == "LQF" || c3 == "QLF") {
      if (sname) {
        return 128;
      } else {
        return 128;
      }
    } else if (c3 == "HRD") {
      if (sname) {
        return 128;
      } else {
        return 128;
      }
    } else if (c3 == "BRD") {
      if (sname) {
        return 128;
      } else {
        return 128;
      }
    }
  } else if (c2 == "SY") {
    if (sname && c3 == "TRD") {
      return 32;
    }
  } else if (cname && c2 == "HE") {
    if (c3 == "TRD") {
      return 32;
    }
  } else if (sname && c2 == "OR") {
    if (c3.at(0) == 'G') {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
        return 128;
      }
    }
  } else if (cname && c2 == "UN") {
    if (c3.at(0) == 'G') {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR")  {
        return 128;
      }
    }
  }
  return 0;
}


/*  -- LAPACK auxiliary routine (version 3.2.1)                        -- */

/*  -- April 2009                                                      -- */

/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ILAENV is called from the LAPACK routines to choose problem-dependent */
/*  parameters for the local environment.  See ISPEC for a description of */
/*  the parameters. */

/*  ILAENV returns an INTEGER */
/*  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC */
/*  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value. */

/*  This version provides a set of parameters which should give good, */
/*  but not optimal, performance on many of the currently available */
/*  computers.  Users are encouraged to modify this subroutine to set */
/*  the tuning parameters for their particular machine using the option */
/*  and problem size information in the arguments. */

/*  This routine will not function correctly if it is converted to all */
/*  lower case.  Converting it to all upper case is allowed. */

/*  Arguments */
/*  ========= */

/*  ISPEC   (input) INTEGER */
/*          Specifies the parameter to be returned as the value of */
/*          ILAENV. */
/*          = 1: the optimal blocksize; if this value is 1, an unblocked */
/*               algorithm will give the best performance. */
/*          = 2: the minimum block size for which the block routine */
/*               should be used; if the usable block size is less than */
/*               this value, an unblocked routine should be used. */
/*          = 3: the crossover point (in a block routine, for N less */
/*               than this value, an unblocked routine should be used) */
/*          = 4: the number of shifts, used in the nonsymmetric */
/*               eigenvalue routines (DEPRECATED) */
/*          = 5: the minimum column dimension for blocking to be used; */
/*               rectangular blocks must have dimension at least k by m, */
/*               where k is given by ILAENV(2,...) and m by ILAENV(5,...) */
/*          = 6: the crossover point for the SVD (when reducing an m by n */
/*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds */
/*               this value, a QR factorization is used first to reduce */
/*               the matrix to a triangular form.) */
/*          = 7: the number of processors */
/*          = 8: the crossover point for the multishift QR method */
/*               for nonsymmetric eigenvalue problems (DEPRECATED) */
/*          = 9: maximum size of the subproblems at the bottom of the */
/*               computation tree in the divide-and-conquer algorithm */
/*               (used by xGELSD and xGESDD) */
/*          =10: ieee NaN arithmetic can be trusted not to trap */
/*          =11: infinity arithmetic can be trusted not to trap */
/*          12 <= ISPEC <= 16: */
/*               xHSEQR or one of its subroutines, */
/*               see IPARMQ for detailed explanation */

/*  NAME    (input) CHARACTER*(*) */
/*          The name of the calling subroutine, in either upper case or */
/*          lower case. */

/*  OPTS    (input) CHARACTER*(*) */
/*          The character options to the subroutine NAME, concatenated */
/*          into a single character string.  For example, UPLO = 'U', */
/*          TRANS = 'T', and DIAG = 'N' for a triangular routine would */
/*          be specified as OPTS = 'UTN'. */

/*  N1      (input) INTEGER */
/*  N2      (input) INTEGER */
/*  N3      (input) INTEGER */
/*  N4      (input) INTEGER */
/*          Problem dimensions for the subroutine NAME; these may not all */
/*          be required. */

/*  Further Details */
/*  =============== */

/*  The following conventions have been used when calling ILAENV from the */
/*  LAPACK routines: */
/*  1)  OPTS is a concatenation of all of the character options to */
/*      subroutine NAME, in the same order that they appear in the */
/*      argument list for NAME, even if they are not used in determining */
/*      the value of the parameter specified by ISPEC. */
/*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order */
/*      that they appear in the argument list for NAME.  N1 is used */
/*      first, N2 second, and so on, and unused problem dimensions are */
/*      passed a value of -1. */
/*  3)  The parameter value returned by ILAENV is checked for validity in */
/*      the calling subroutine.  For example, ILAENV is used to retrieve */
/*      the optimal blocksize for STRTRI as follows: */

/*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 ) */
/*      IF( NB.LE.1 ) NB = MAX( 1, N ) */

/*  ===================================================================== */
inline int ilaenv(int ispec, const std::string& name, int n1, int n2, int n3, int n4) {

  if (ispec < 1 || ispec > 3) {
    switch (ispec) {
      case 4:  return 6;  /* ISPEC = 4:  number of shifts (used by xHSEQR) */
      case 5:  return 2;  /* ISPEC = 5:  minimum column dimension (not used) */
      case 6:             /* ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */
        return (int) ((float)std::min(n1, n2) * 1.6f);
      case 7:  return 1;  /* ISPEC = 7:  number of processors (not used) */
      case 8:  return 50; /* ISPEC = 8:  crossover point for multishift (used by xHSEQR) */
      case 9:  return 25; /* ISPEC = 9:  maximum size of the subproblems at the bottom of the */
                          /*             computation tree in the divide-and-conquer algorithm */
                          /*             (used by xGELSD and xGESDD) */
      case 10:            /* ISPEC = 10: ieee NaN arithmetic can be trusted not to trap */
        return ieeeck(1);

      case 11:            /* ISPEC = 11: infinity arithmetic can be trusted not to trap */
         return ieeeck(0);

      default:
        if (ispec >= 12 && ispec <= 16) { /* 12 <= ISPEC <= 16: xHSEQR or one of its subroutines. */
          return iparmq(ispec, n2, n3);
        } else {
          return -1; /* Invalid value for ISPEC */
        }
      }
    }


/*     Convert NAME to upper case if the first character is lower case. */

    std::string subnam(name);
    std::transform(subnam.begin(), subnam.end(), subnam.begin(), ::toupper);
    std::string c1(subnam);

    bool sname = c1.at(0) == 'S' || c1.at(0) == 'D',
         cname = c1.at(0) == 'C' || c1.at(0) == 'Z';

    if (! (cname || sname)) return 1;

    std::string c2(subnam.substr(1, 2)),
                c3(subnam.substr(3, 3)),
                c4(c3.substr(1, 2));

    if (ispec == 2) return ilaenv_min_block_size(c2, c3, c4, sname, cname);
    if (ispec == 3) return ilaenv_crossover_point(c2, c3, c4, sname, cname);
    return ilaenv_block_size(n2, n4, c2, c3, c4, sname, cname);

} /* ilaenv_ */



/* > \brief \b DGETF2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGETF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetf2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetf2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetf2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGETF2 computes an LU factorization of a general m-by-n matrix A */
/* > using partial pivoting with row interchanges. */
/* > */
/* > The factorization has the form */
/* >    A = P * L * U */
/* > where P is a permutation matrix, L is lower triangular with unit */
/* > diagonal elements (lower trapezoidal if m > n), and U is upper */
/* > triangular (upper trapezoidal if m < n). */
/* > */
/* > This is the right-looking Level 2 BLAS version of the algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the m by n matrix to be factored. */
/* >          On exit, the factors L and U from the factorization */
/* >          A = P*L*U; the unit diagonal elements of L are not stored. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (min(M,N)) */
/* >          The pivot indices; for 1 <= i <= min(M,N), row i of the */
/* >          matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -k, the k-th argument had an illegal value */
/* >          > 0: if INFO = k, U(k,k) is exactly zero. The factorization */
/* >               has been completed, but the factor U is exactly */
/* >               singular, and division by zero will occur if it is used */
/* >               to solve a system of equations. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleGEcomputational */

/*  ===================================================================== */

template <typename DType>
inline int getf2(const int m, const int n, DType* a, const int lda, int *ipiv) {

  /* Function Body */
  if (m < 0)                      return -1; // error
  else if (n < 0)                 return -2; // error
  else if (lda < std::max(1,m))   return -4; // error


  if (m == 0 || n == 0)     return 0;   /* Quick return if possible */

  for (size_t j = 0; j < std::min(m,n); ++j) { // changed

    /* Find pivot and test for singularity. */

    int jp = j - 1 + idamax<DType>(m-j+1, &a[j + j * lda], 1);

    ipiv[j] = jp;


    if (a[jp + j*lda] != 0) {

      /* Apply the interchange to columns 1:N. */
      // (Don't swap two columns that are the same.)
      if (jp != j) swap<DType>(n, &a[j], lda, &a[jp], lda);

      /* Compute elements J+1:M of J-th column. */

	    if (j < m-1) {
        if (std::abs(a[j+j*lda]) >= std::numeric_limits<DType>::min()) {
          scal<DType>(m-j, 1.0 / a[j+j*lda], &a[j+1+j*lda], 1);
		    } else {
		      for (size_t i = 0; i < m-j; ++i) { // changed
			      a[j+i+j*lda] /= a[j+j*lda];
		      }
		    }
	    }

    } else { // singular matrix
      return j; // U(j,j) is exactly zero, div by zero if answer is used to solve a system of equations.
    }

    if (j < std::min(m,n)-1) /*           Update trailing submatrix. */
      ger<DType>(m-j, n-j, -1.0, &a[j+1+j*lda], 1, &a[j+(j+1)*lda], lda, &a[j+1+(j+1)*lda], lda);

  }
  return 0;
} /* dgetf2_ */



/*
 * Templated version of getrf, derived from dgetrf.c (from LAPACK 3.2) using f2c.
 */
template <typename DType>
inline int getrf(int m, int n, DType *a, int lda, int *ipiv) {

  int info = 0;
  if (m < 0)                        return -1;
  else if (n < 0)                   return -2;
  else if (lda < std::max(1,m))     return -4;

  /*     Quick return if possible */

  if (m == 0 || n == 0) return 0;

  /*     Determine the block size for this environment. */

  int nb = ilaenv(1, "DGETRF", m, n, -1, -1);

  if (nb <= 1 || nb >= std::min(m,n)) { /*        Use unblocked code. */

    info = getf2<DType>(m, n, a, lda, ipiv);

  } else { /*        Use blocked code. */

    for (int j = 0; nb < 0 ? j > std::min(m,n) : j < std::min(m,n); j += nb) {

      int jb = std::min(std::min(m,n) - j + 1, nb);

      /*           Factor diagonal and subdiagonal blocks and test for exact */
      /*           singularity. */

      int iinfo = getf2<DType>(m-j+1, jb, &a[j + j * lda], lda, &ipiv[j]);

      /*           Adjust INFO and the pivot indices. */

      if (!info && iinfo > 0) {
        info = iinfo + j - 1;
      }

      /* Computing MIN */
      for (size_t i = j; i < std::min(m,j+jb-1); ++i) {
        ipiv[i] = j - 1 + ipiv[i];
      }

      /*           Apply interchanges to columns 1:J-1. */

      laswp<DType>(j-1, a, lda, j, j+jb-1, ipiv, 1);

      if (j + jb < n) {

        /* Apply interchanges to columns J+JB:N. */

        laswp<DType>(n-j-jb+1, &a[(j + jb) * lda], lda, j, j+jb-1, ipiv, 1);

        /* void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                            const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                            const enum CBLAS_DIAG Diag, const int M, const int N,
                            const DType alpha, const DType *A, const int lda,
                            DType *B, const int ldb); */
        nm::math::trsm<DType>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, jb, n-j-jb+1, 1.0, // M, N, alpha
          &a[j+j*lda], lda, &a[j+(j+jb)*lda], lda); // *A, lda, *B, ldb

        if (j + jb < m) {

          /*                 Update trailing submatrix. */
          gemm<DType>(CblasRowMajor, CblasNoTrans, CblasNoTrans, m-j-jb+1, n-j-jb+1, jb,
                      -1.0, &a[j+jb+j*lda], lda, &a[j+(j+jb)*lda], lda, 1.0, &a[j+jb+(j+jb)*lda], lda);
        }
      }
    }
  }
  return 0;
} /* dgetrf_ */




} // end namespace lapack

/*
 * Function signature conversion for LAPACK's scal function.
 */
template <typename DType>
inline void clapack_scal(const int n, const void* da, void* dx, const int incx) {
  // FIXME: See if we can call the clapack version instead of our C++ version.
  nm::math::lapack::scal<DType>(n, *reinterpret_cast<const DType*>(da), reinterpret_cast<DType*>(dx), incx);
}


}}

#endif