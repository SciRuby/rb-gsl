#ifndef CBLAS_C
# define CBLAS_C

#include <cblas.h>
#include "types.h"

//extern enum CBLAS_ORDER;
//extern enum CBLAS_TRANSPOSE;

inline void cblas_sgemm_( const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                         const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                         const int K, const double alpha, const void* A,
                         const int lda, const void* B, const int ldb,
                         const double beta, void* C, const int ldc)
{
  cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

/*inline void cblas_dgemm_( const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                   const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const double alpha, const void* A,
                   const int lda, const void* B, const int ldb,
                   const double beta, void* C, const int ldc)
{
  cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}*/


inline void cblas_cgemm_( const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                   const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const complex128 alpha, const void* A,
                   const int lda, const void* B, const int ldb,
                   const complex128 beta, void* C, const int ldc)
{
  complex64 alpha_, beta_;

  alpha_.r = (float)(alpha.r);
  alpha_.i = (float)(alpha.i);

  beta_.r  = (float)(beta.r);
  beta_.i  = (float)(beta.i);

  cblas_cgemm(Order, TransA, TransB, M, N, K, &alpha_, A, lda, B, ldb, &beta_, C, ldc);
}


inline void cblas_zgemm_( const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                   const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                   const int K, const complex128 alpha, const void* A,
                   const int lda, const void* B, const int ldb,
                   const complex128 beta, void* C, const int ldc)
{
  cblas_zgemm(Order, TransA, TransB, M, N, K, &alpha, A, lda, B, ldb, &beta, C, ldc);
}


#endif