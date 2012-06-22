
int gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const TYPE* alpha,
  const TYPE* A, const int lda, const TYPE* X, const int incX, const TYPE* beta, TYPE* Y, const int incY)
{
  int lenX, lenY, i, j;
  int kx, ky, iy, jx, jy, ix;
  LONG_TYPE temp;

  // Test the input parameters
  if (Trans < 111 || Trans > 113) {
    rb_raise(rb_eArgError, "GEMV: TransA must be CblasNoTrans, CblasTrans, or CblasConjTrans");
    return 0;
  } else if (lda < NM_MAX(1, N)) {
    fprintf(stderr, "GEMV: N = %d; got lda=%d", N, lda);
    rb_raise(rb_eArgError, "GEMV: Expected lda >= max(1, N)");
    return 0;
  } else if (incX == 0) {
    rb_raise(rb_eArgError, "GEMV: Expected incX != 0\n");
    return 0;
  } else if (incY == 0) {
    rb_raise(rb_eArgError, "GEMV: Expected incY != 0\n");
    return 0;
  }

  // Quick return if possible
  if (!M || !N || *alpha == 0 && *beta == 1) return 0;

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

  if (*alpha == 0) return 0;

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

    // Form  y := alpha*A**T*x + y.
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

  return 0;
}  // end of GEMV
