
int gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const TYPE* alpha,
  const TYPE* A, const int lda, const TYPE* B, const int ldb,
  const TYPE* beta, TYPE* C, const int ldc) {

  if (sizeof(TYPE) == sizeof(float))
    cblas_sgemm(CblasRowMajor, TransB, TransA, N, M, K, *alpha, B, ldb, A, lda, *beta, C, ldc);
  else if (sizeof(TYPE) == sizeof(double))
    cblas_dgemm(CblasRowMajor, TransB, TransA, N, M, K, *alpha, B, ldb, A, lda, *beta, C, ldc);
  else
    rb_raise(rb_eStandardError, "unrecognized type");

  return 0;
  // CblasRowMajor, TransA, TransB, M, N, K, *alpha, A, lda, B, ldb, *beta, C, ldc
}
