
int gemv(const enum CBLAS_TRANSPOSE Trans, const int M, const int N, const TYPE* alpha,
  const TYPE* A, const int lda, const TYPE* X, const int incX, const TYPE* beta, TYPE* Y, const int incY)
{
  if (sizeof(TYPE) == sizeof(float))
    cblas_sgemv(CblasRowMajor, Trans, M, N, *alpha, A, lda, X, incX, *beta, Y, incY);
  else if (sizeof(TYPE) == sizeof(double))
    cblas_dgemv(CblasRowMajor, Trans, M, N, *alpha, A, lda, X, incX, *beta, Y, incY);
  else
    rb_raise(rb_eStandardError, "unrecognized type");

  return 0;
}
