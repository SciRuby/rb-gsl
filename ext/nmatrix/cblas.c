#ifndef CBLAS_C
# define CBLAS_C

#include <cblas.h>
#include "nmatrix.h"

//extern enum CBLAS_ORDER;
//extern enum CBLAS_TRANSPOSE;

inline CBLAS_PARAM init_cblas_params_for_nm_multiply_matrix(int8_t dtype) {
  CBLAS_PARAM p;

  // memset(&(p.beta), 0, sizeof(nm_size128_t)); // set beta to 0

  switch(dtype) {
  case NM_FLOAT32:
  case NM_FLOAT64:
    p.alpha.d[0] = 1.0;
    p.beta.d[0] = 0.0;
    break;

  case NM_COMPLEX64:
    p.alpha.c[0].r = 1.0;
    p.alpha.c[0].i = 0.0;
    p.beta.c[0].r = 0.0;
    p.beta.c[0].i = 0.0;
    break;

  case NM_COMPLEX128:
    p.alpha.z.r = 1.0;
    p.alpha.z.i = 0.0;
    p.beta.z.r = 0.0;
    p.beta.z.i = 0.0;
    break;
  }

  return p;
}


inline void cblas_sgemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, CBLAS_PARAM p) {
  cblas_sgemm(Order, TransA, TransB, p.M, p.N, p.K, p.alpha.d[0], p.A, p.lda, p.B, p.ldb, p.beta.d[0], p.C, p.ldc);
}

inline void cblas_dgemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, CBLAS_PARAM p) {
  cblas_dgemm(Order, TransA, TransB, p.M, p.N, p.K, p.alpha.d[0], p.A, p.lda, p.B, p.ldb, p.beta.d[0], p.C, p.ldc);
}

inline void cblas_cgemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, CBLAS_PARAM p) {
  cblas_cgemm(Order, TransA, TransB, p.M, p.N, p.K, &(p.alpha.c), p.A, p.lda, p.B, p.ldb, &(p.beta.c), p.C, p.ldc);
}

inline void cblas_zgemm_(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, CBLAS_PARAM p) {
  cblas_zgemm(Order, TransA, TransB, p.M, p.N, p.K, &(p.alpha.z), p.A, p.lda, p.B, p.ldb, &(p.beta.z), p.C, p.ldc);
}


#endif