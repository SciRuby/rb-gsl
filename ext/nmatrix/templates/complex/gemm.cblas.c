
int gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const TYPE* alpha,
  const TYPE* A, const int lda, const TYPE* B, const int ldb,
  const TYPE* beta, TYPE* C, const int ldc) {

  if (sizeof(TYPE) == sizeof(float)*2)
    cblas_cgemm(CblasRowMajor, TransB, TransA, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
  else if (sizeof(TYPE) == sizeof(double)*2)
    cblas_zgemm(CblasRowMajor, TransB, TransA, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc);
  else
    rb_raise(rb_eStandardError, "unrecognized type");
  // CblasRowMajor, TransB, TransA, N, M, K, alpha, B, ldb, A, lda, beta, C, ldc
  return 0;
}
