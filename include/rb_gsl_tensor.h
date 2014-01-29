#ifndef ___RB_TENSOR_H___
#define ___RB_TENSOR_H___

#ifdef HAVE_TENSOR_TENSOR_H
#include "rb_gsl.h"
#include <tensor/tensor.h>

EXTERN VALUE cgsl_tensor, cgsl_tensor_int;

enum {
  TENSOR_ADD,
  TENSOR_SUB,
  TENSOR_MUL_ELEMENTS,
  TENSOR_DIV_ELEMENTS,
  TENSOR_SCALE,
  TENSOR_ADD_CONSTANT,
  TENSOR_ADD_DIAGONAL,
  TENSOR_PRODUCT,
  TENSOR_CONTRACT,
};

typedef gsl_permutation tensor_indices;

typedef struct __rbgsl_tensor {
  /*  tensor *tensor;
      tensor_indices *indices;*/
  tensor *tensor;
  tensor_indices *indices;
} rbgsl_tensor;

typedef struct __rbgsl_tensor_int {
  /*  tensor_int *tensor;
      tensor_indices *indices;*/
  tensor_int *tensor;
  tensor_indices *indices;
} rbgsl_tensor_int;

rbgsl_tensor* rbgsl_tensor_alloc(const unsigned int rank, const size_t dimension);
rbgsl_tensor_int* rbgsl_tensor_int_alloc(const unsigned int rank, const size_t dimension);
void rbgsl_tensor_free(rbgsl_tensor*);
void rbgsl_tensor_int_free(rbgsl_tensor_int*);
#endif
#endif
