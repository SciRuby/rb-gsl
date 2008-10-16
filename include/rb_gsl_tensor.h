#ifndef ___RB_GSL_TENSOR_H___
#define ___RB_GSL_TENSOR_H___

#ifdef HAVE_GSL_TENSOR_GSL_TENSOR_H
#include "rb_gsl.h"
#include "rb_gsl_common.h"
#include <gsl_tensor/gsl_tensor.h>

EXTERN VALUE cgsl_tensor, cgsl_tensor_int;

enum {
  GSL_TENSOR_ADD,
  GSL_TENSOR_SUB,
  GSL_TENSOR_MUL_ELEMENTS,
  GSL_TENSOR_DIV_ELEMENTS,
  GSL_TENSOR_SCALE,
  GSL_TENSOR_ADD_CONSTANT,
  GSL_TENSOR_ADD_DIAGONAL,
  GSL_TENSOR_PRODUCT,
  GSL_TENSOR_CONTRACT,
};

typedef gsl_permutation gsl_tensor_indices;

typedef struct __rbgsl_tensor {
  gsl_tensor *tensor;
  gsl_tensor_indices *indices;
} rbgsl_tensor;

typedef struct __rbgsl_tensor_int {
  gsl_tensor_int *tensor;
  gsl_tensor_indices *indices;
} rbgsl_tensor_int;

rbgsl_tensor* rbgsl_tensor_alloc(const unsigned int rank, const size_t dimension);
rbgsl_tensor_int* rbgsl_tensor_int_alloc(const unsigned int rank, const size_t dimension);
void rbgsl_tensor_free(rbgsl_tensor*);
void rbgsl_tensor_int_free(rbgsl_tensor_int*);
#endif
#endif
