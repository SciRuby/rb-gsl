
int %%INT_ABBREV%%gemv(enum CBLAS_TRANSPOSE Trans, const size_t M, const size_t N, const %%INT%% alpha,
  const %%INT%%* A, const size_t lda, const %%INT%%* X, const int incX, const %%INT%% beta, %%INT%%* Y, const int incY)
{
  size_t lenX, lenY, i, j;
  int kx, ky, iy, jx, jy, ix;
  rational128 temp, temp2, temp3;

  // Test the input parameters
  if (Trans < 111 || Trans > 113) {
    fprintf(stderr, "IGEMV: TransA must be CblasNoTrans, CblasTrans, or CblasConjTrans\n");
    return 0;
  } else if (lda < NM_MAX(1, N)) {
    fprintf(stderr, "IGEMV: Expected lda >= max(1, N), with N = %d; got lda=%d\n", N, lda);
    return 0;
  } else if (incX == 0) {
    fprintf(stderr, "IGEMV: Expected incX != 0\n");
    return 0;
  } else if (incY == 0) {
    fprintf(stderr, "IGEMV: Expected incY != 0\n");
    return 0;
  }

  // Quick return if possible
  if (!M || !N || !alpha.n && beta.n == beta.d) return 0;

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
  if (beta.n != beta.d) {
    if (incY == 1) {
      if (!beta.n) {
        for (i = 0; i < lenY; ++i) {
          Y[i].n = 0;
          Y[i].d = 1;
        }
      } else {
        for (i = 0; i < lenY; ++i)
          Y[i] = %%INT_ABBREV%%_muldiv(beta.n, beta.d, Y[i].n, Y[i].d, '*');
      }
    } else {
      iy = ky;
      if (!beta.n) {
        for (i = 0; i < lenY; ++i) {
          Y[iy].n = 0;
          Y[iy].d = 1;
          iy += incY;
        }
      } else {
        for (i = 0; i < lenY; ++i) {
          Y[iy] = %%INT_ABBREV%%_muldiv(beta.n, beta.d, Y[iy].n, Y[iy].d, '*');
          iy += incY;
        }
      }
    }
  }

  if (!alpha.n) return 0;

  if (Trans == CblasNoTrans) {

    // Form  y := alpha*A*x + y.
    jx = kx;
    if (incY == 1) {
      for (j = 0; j < N; ++j) {
        if (X[jx].n != 0) {
          temp = r128_muldiv(alpha.n, alpha.d, X[jx].n, X[jx].d, '*');
          for (i = 0; i < M; ++i) {
            temp2 = r128_muldiv(A[j+i*lda].n, A[j+i*lda].d, temp.n, temp.d, '*');
            Y[i]  = %%INT_ABBREV%%_addsub(Y[i].n, Y[i].d, temp2.n, temp2.d, '+');
          }
        }
        jx += incX;
      }
    } else {
      for (j = 0; j < N; ++j) {
        if (X[jx].n != 0) {
          temp = r128_muldiv(alpha.n, alpha.d, X[jx].n, X[jx].d, '*');
          iy = ky;
          for (i = 0; i < M; ++i) {
            temp2 = r128_muldiv(A[j+i*lda].n, A[j+i*lda].d, temp.n, temp.d, '*');
            Y[iy] = %%INT_ABBREV%%_addsub(Y[iy].n, Y[iy].d, temp2.n, temp2.d, '*');
            iy += incY;
          }
        }
        jx += incX;
      }
    }

  } else { // TODO: Check that indices are correct! They're switched for C.

    // Form  y := alpha*A**T*x + y.

    jy = ky;

    if (incX == 1) {
      for (j = 0; j < N; ++j) {
        temp.n = 0;
        temp.d = 1;
        for (i = 0; i < M; ++i) {
          temp2 = r128_muldiv(A[j+i*lda].n, A[j+i*lda].d, X[j].n, X[j].d, '*');
          temp  = r128_addsub(temp.n, temp.d, temp2.n, temp2.d, '+');
        }
        temp3 = r128_muldiv(alpha.n, alpha.d, temp.n, temp.d, '*');
        Y[jy] = %%INT_ABBREV%%_addsub(temp3.n, temp3.d, Y[jy].n, Y[jy].d, '+');
        jy += incY;
      }
    } else {
      for (j = 0; j < N; ++j) {
        temp.n = 0;
        temp.d = 1;
        ix = kx;
        for (i = 0; i < M; ++i) {
          temp2 = r128_muldiv(A[j+i*lda].n, A[j+i*lda].d, X[ix].n, X[ix].d, '*');
          temp  = r128_addsub(temp2.n, temp2.d, temp.n, temp.d, '+');
          ix += incX;
        }
        temp3 = r128_muldiv(alpha.n, alpha.d, temp.n, temp.d, '*');
        Y[jy] = %%INT_ABBREV%%_addsub(temp3.n, temp3.d, Y[jy].n, Y[jy].d, '+');
        jy += incY;
      }
    }
  }

  return 0;
  // end of GEMV
}
