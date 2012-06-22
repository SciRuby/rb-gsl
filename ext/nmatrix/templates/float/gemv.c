
int gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const TYPE* alpha,
  const TYPE* A, const int lda, const TYPE* X, const int incX, const TYPE* beta, TYPE* Y, const int incY)
{
  if (sizeof(TYPE) == sizeof(float)*2)
    return cblas_cgemv(CblasRowMajor, TransA, TransB, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
  else if (sizeof(TYPE) == sizeof(double)*2)
    return cblas_zgemv(CblasRowMajor, TransA, TransB, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}
