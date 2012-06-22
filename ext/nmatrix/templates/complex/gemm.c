
int gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const TYPE* alpha,
  const TYPE* A, const int lda, const TYPE* B, const int ldb,
  const TYPE* beta, TYPE* C, const int ldc) {

  if (sizeof(TYPE) == sizeof(float)*2)
    return cblas_cgemm(CblasRowMajor, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  else if (sizeof(TYPE) == sizeof(double)*2)
    return cblas_zgemm(CblasRowMajor, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
