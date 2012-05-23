
int gemm(int TransA, int TransB, const int M, const int N, const int K, const TYPE alpha,
  const TYPE* A, const int lda, const TYPE* B, const int ldb,
  const TYPE beta, TYPE* C, const int ldc)
{
  int num_rows_a, /*num_cols_a,*/ num_rows_b; // nrowa, ncola, nrowb

  // use longest possible type for intermediate value storage:
  LONG_TYPE temp;
  // %%= if [:rational,:complex,:value].include?(dtype.type); "#{dtype.long_dtype.sizeof} temp1, temp2;"; end%%
  int i, j, l;

  if (TransA == CblasNoTrans) num_rows_a = M;
  else                        num_rows_a = K;

  if (TransB == CblasNoTrans) num_rows_b = K;
  else                        num_rows_b = N;

  // Test the input parameters
  if (TransA < 111 || TransA > 113) {
    fprintf(stderr, "GEMM: TransA must be CblasNoTrans, CblasTrans, or CblasConjTrans\n");
    return 0;
  } else if (TransB < 111 || TransB > 113) {
    fprintf(stderr, "GEMM: TransB must be CblasNoTrans, CblasTrans, or CblasConjTrans\n");
    return 0;
  } else if (M < 0) {
    fprintf(stderr, "GEMM: Expected M >= 0\n");
    return 0;
  } else if (N < 0) {
    fprintf(stderr, "GEMM: Expected N >= 0\n");
    return 0;
  } else if (K < 0) {
    fprintf(stderr, "GEMM: Expected K >= 0\n");
    return 0;
  } else if (lda < NM_MAX(1, num_rows_a)) {
    fprintf(stderr, "GEMM: Expected lda >= max(1, num_rows_a), with num_rows_a = %d; got lda=%d\n", num_rows_a, lda);
    return 0;
  } else if (ldb < NM_MAX(1, num_rows_b)) {
    fprintf(stderr, "GEMM: Expected ldb >= max(1, num_rows_b), with num_rows_b = %d; got ldb=%d\n", num_rows_b, ldb);
    return 0;
  } else if (ldc < NM_MAX(1,M)) {
    fprintf(stderr, "GEMM: Expected ldc >= max(1,M) with M=%d; got ldc=%d\n", M, ldc);
    return 0;
  }

  // Quick return if possible
  if (!M || !N || (alpha == 0 || !K) && beta == 1) return 0;

  // For alpha = 0
  if (alpha == 0) {
    if (beta == 0) {
      for (j = 0; j < N; ++j)
        for (i = 0; i < M; ++i) {
          C[i+j*ldc] = 0;
        }
    } else {
      for (j = 0; j < N; ++j)
        for (i = 0; i < M; ++i) {
          C[i+j*ldc] *= beta;
        }
    }
    return 0;
  }

  // Start the operations
  if (TransB == CblasNoTrans) {
    if (TransA == CblasNoTrans) {
      // C = alpha*A*B+beta*C
      for (j = 0; j < N; ++j) {
        if (beta == 0) {
          for (i = 0; i < M; ++i) {
            C[i+j*ldc] = 0;
          }
        } else if (beta != 1) {
          for (i = 0; i < M; ++i) {
            C[i+j*ldc] *= beta;
          }
        }

        for (l = 0; l < K; ++l) {
          if (B[l+j*ldb] != 0) {
            temp = alpha * B[l+j*ldb];
            for (i = 0; i < M; ++i) {
              C[i+j*ldc] += A[i+l*lda] * temp;
            }
          }
        }
      }

    } else {

      // C = alpha*A**T*B + beta*C
      for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {
          temp = 0;
          for (l = 0; l < K; ++l) {
            temp += A[l+i*lda] * B[l+j*ldb];
          }

          if (beta == 0) {
            C[i+j*ldc] = alpha*temp;
          } else {
            C[i+j*ldc] = alpha*temp + beta*C[i+j*ldc];
          }
        }
      }

    }

  } else if (TransA == CblasNoTrans) {

    // C = alpha*A*B**T + beta*C
    for (j = 0; j < N; ++j) {
      if (beta == 0) {
        for (i = 0; i < M; ++i) {
          C[i+j*ldc] = 0;
        }
      } else if (beta != 1) {
        for (i = 0; i < M; ++i) {
          C[i+j*ldc] *= beta;
        }
      }

      for (l = 0; l < K; ++l) {
        if (B[j+l*ldb] != 0) {
          temp = alpha * B[j+l*ldb];
          for (i = 0; i < M; ++i) {
            C[i+j*ldc] += A[i+l*lda] * temp;
          }
        }
      }

    }

  } else {

    // C = alpha*A**T*B**T + beta*C
    for (j = 0; j < N; ++j) {
      for (i = 0; i < M; ++i) {
        temp = 0;
        for (l = 0; l < K; ++l) {
          temp += A[l+i*lda] * B[j+l*ldb];
        }

        if (beta == 0) {
          C[i+j*ldc] = alpha*temp;
        } else {
          C[i+j*ldc] = alpha*temp + beta*C[i+j*ldc];
        }
      }
    }

  }

  return 0;
}
