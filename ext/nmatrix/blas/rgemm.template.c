
int %%INT_ABBREV%%gemm(enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB,
  const int M, const int N, const int K, const %%INT%% alpha,
  const %%INT%%* A, const int lda,
  const %%INT%%* B, const int ldb, const %%INT%% beta,
  %%INT%%* C, const int ldc)
{
  int num_rows_a, /*num_cols_a,*/ num_rows_b; // nrowa, ncola, nrowb
  rational128 temp, temp2, temp3; // use int64 since we're accumulating here, and it may well exceed the max int size.
  int i, j, l;

  if (TransA == CblasNoTrans) {
    num_rows_a = M;
    //num_cols_a = K;
  } else {
    num_rows_a = K;
    //num_cols_a = M;
  }

  if (TransB == CblasNoTrans) num_rows_b = K;
  else                        num_rows_b = N;

  // Test the input parameters
  if (TransA < 111 || TransA > 113) {
    fprintf(stderr, "IGEMM: TransA must be CblasNoTrans, CblasTrans, or CblasConjTrans\n");
    return 0;
  } else if (TransB < 111 || TransB > 113) {
    fprintf(stderr, "IGEMM: TransB must be CblasNoTrans, CblasTrans, or CblasConjTrans\n");
    return 0;
  } else if (M < 0) {
    fprintf(stderr, "IGEMM: Expected M >= 0\n");
    return 0;
  } else if (N < 0) {
    fprintf(stderr, "IGEMM: Expected N >= 0\n");
    return 0;
  } else if (K < 0) {
    fprintf(stderr, "IGEMM: Expected K >= 0\n");
    return 0;
  } else if (lda < NM_MAX(1, num_rows_a)) {
    fprintf(stderr, "IGEMM: Expected lda >= max(1, num_rows_a), with num_rows_a = %d; got lda=%d\n", num_rows_a, lda);
    return 0;
  } else if (ldb < NM_MAX(1, num_rows_b)) {
    fprintf(stderr, "IGEMM: Expected ldb >= max(1, num_rows_b), with num_rows_b = %d; got ldb=%d\n", num_rows_b, ldb);
    return 0;
  } else if (ldc < NM_MAX(1,M)) {
    fprintf(stderr, "IGEMM: Expected ldc >= max(1,M) with M=%d; got ldc=%d\n", M, ldc);
    return 0;
  }

  // Quick return if possible
  if (!M || !N || (!alpha.n || !K) && (beta.n == 1 && beta.d == 1)) return 0;

  // For alpha = 0
  if (alpha.n == 0) {
    if (beta.n == 0) {
      for (j = 0; j < N; ++j)
        for (i = 0; i < M; ++i) {
          C[i+j*ldc].n = 0;
          C[i+j*ldc].d = 1;
        }
    } else {
      for (j = 0; j < N; ++j)
        for (i = 0; i < M; ++i)
          C[i+j*ldc] = %%INT_ABBREV%%_muldiv(beta.n, beta.d, C[i+j*ldc].n, C[i+j*ldc].d, '*');
    }
    return 0;
  }

  // Start the operations
  if (TransB == CblasNoTrans) {
    if (TransA == CblasNoTrans) {
      // C = alpha*A*B+beta*C
      for (j = 0; j < N; ++j) {
        if (beta.n == 0) {
          for (i = 0; i < M; ++i) {
            C[i+j*ldc].n = 0;
            C[i+j*ldc].d = 1;
          }
        } else if (!(beta.n == 1 && beta.d == 1)) {
          for (i = 0; i < M; ++i)
            C[i+j*ldc] = %%INT_ABBREV%%_muldiv(beta.n, beta.d, C[i+j*ldc].n, C[i+j*ldc].d, '*');
        }

        for (l = 0; l < K; ++l) {
          if (B[l+j*ldb].n != 0) {
            temp = r128_muldiv(alpha.n, alpha.d, B[l+j*ldb].n, B[l+j*ldb].d, '*');
            for (i = 0; i < M; ++i)
              C[i+j*ldc] = %%INT_ABBREV%%_muldiv(temp.n, temp.d, A[i+l*lda].n, A[i+l*lda].d, '*');
          }
        }
      }

    } else {

      // C = alpha*A**T*B + beta*C
      for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {

          temp.n = 0;
          temp.d = 1;

          for (l = 0; l < K; ++l) { //temp += A[l+i*lda] * B[l+j*ldb];
            temp2 = r128_muldiv(A[l+i*lda].n, A[l+i*lda].d, B[l+j*ldb].n, B[l+j*ldb].d, '*');
            temp  = r128_addsub(temp.n, temp.d, temp2.n, temp2.d, '+');
          }

          if (beta.n == 0)  C[i+j*ldc] = %%INT_ABBREV%%_muldiv(alpha.n, alpha.d, temp.n, temp.d, '*');
          else { // C[i+j*ldc] = alpha*temp + beta*C[i+j*ldc];
            temp2 = r128_muldiv(alpha.n, alpha.d, temp.n, temp.d, '*');
            temp3 = r128_muldiv(beta.n, beta.d, C[i+j*ldc].n, C[i+j*ldc].d, '*');
            C[i+j*ldc] = %%INT_ABBREV%%_addsub(temp2.n, temp2.d, temp3.n, temp3.d, '+');
          }
        }
      }

    }

  } else if (TransA == CblasNoTrans) {

    // C = alpha*A*B**T + beta*C
    for (j = 0; j < N; ++j) {
      if (beta.n == 0) {
        for (i = 0; i < M; ++i) {
          C[i+j*ldc].n = 0;
          C[i+j*ldc].d = 1;
        }
      } else if (!(beta.n == 1 && beta.d == 1)) {
        for (i = 0; i < M; ++i)
          C[i+j*ldc] = %%INT_ABBREV%%_muldiv(beta.n, beta.d, C[i+j*ldc].n, C[i+j*ldc].d, '*');
      }

      for (l = 0; l < K; ++l) {
        if (B[j+l*ldb].n != 0) {
          //temp = alpha * B[j+l*ldb];
          temp = r128_muldiv(alpha.n, alpha.d, B[j+l*ldb].n, B[j+l*ldb].d, '*');
          for (i = 0; i < M; ++i) {
            temp2 = r128_muldiv(A[i+l*lda].n, A[i+l*lda].d, temp.n, temp.d, '*');
            C[i+j*ldc] = %%INT_ABBREV%%_addsub(C[i+j*ldc].n, C[i+j*ldc].d, temp2.n, temp2.d, '+');
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
          temp2 = r128_muldiv(A[l+i*lda].n, A[l+i*lda].d, B[j+l*ldb].n, B[j+l*ldb].d, '*');
          temp  = r128_addsub(temp.n, temp.d, temp2.n, temp2.d, '+');
        }

        if (beta.n == 0) C[i+j*ldc] = %%INT_ABBREV%%_muldiv(alpha.n, alpha.d, temp.n, temp.d, '*'); //alpha*temp;
        else {
          temp2 = r128_muldiv(alpha.n, alpha.d, temp.n, temp.d, '*');
          temp3 = r128_muldiv(beta.n, beta.d, C[i+j*ldc].n, C[i+j*ldc].d, '*');
          C[i+j*ldc] = %%INT_ABBREV%%_addsub(temp2.n, temp2.d, temp3.n, temp3.d, '+');
        }
      }
    }
  }

  return 0;
}

