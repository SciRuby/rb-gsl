
int gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const TYPE* alpha,
  const TYPE* A, const int lda, const TYPE* X, const int incX, const TYPE* beta, TYPE* Y, const int incY)
{
  if (sizeof(TYPE) == sizeof(float))
    return cblas_sgemv(CblasRowMajor, TransA, TransB, M, N, *alpha, A, lda, B, ldb, *beta, C, ldc);
  else if (sizeof(TYPE) == sizeof(double))
    return cblas_dgemv(CblasRowMajor, TransA, TransB, M, N, *alpha, A, lda, B, ldb, *beta, C, ldc);
}