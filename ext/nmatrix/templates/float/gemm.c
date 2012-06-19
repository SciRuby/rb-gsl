
int gemm(const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const TYPE* alpha,
  const TYPE* A, const int lda, const TYPE* B, const int ldb,
  const TYPE* beta, TYPE* C, const int ldc) {

  if (sizeof(TYPE) == sizeof(float))
    return cblas_sgemm(CblasRowMajor, TransA, TransB, M, N, K, *alpha, A, lda, B, ldb, *beta, C, ldc);
  else if (sizeof(TYPE) == sizeof(double))
    return cblas_dgemm(CblasRowMajor, TransA, TransB, M, N, K, *alpha, A, lda, B, ldb, *beta, C, ldc);
}
