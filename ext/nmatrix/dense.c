#ifndef DENSE_C
#define DENSE_C

#include <ruby.h>

#include "nmatrix.h"


/* Calculate the number of elements in the dense storage structure, based on shape and rank */
size_t count_dense_storage_elements(DENSE_STORAGE* s) {
  size_t i;
  size_t count = 1;
  for (i = 0; i < s->rank; ++i) count *= s->shape[i];
  return count;
}


size_t dense_storage_pos(DENSE_STORAGE* s, size_t* coords) {
  size_t k, l;
  size_t inner, outer = 0;
  for (k = 0; k < s->rank; ++k) {
    inner = coords[k];
    for (l = k+1; l < s->rank; ++l) {
      inner *= s->shape[l];
    }
    outer += inner;
  }
  return outer;
}


void* dense_storage_get(DENSE_STORAGE* s, size_t* coords, size_t elem_size) {
  return (char*)(s->elements) + dense_storage_pos(s, coords) * elem_size;
}


/* Does not free passed-in value! Different from list_storage_insert. */
void dense_storage_set(DENSE_STORAGE* s, size_t* coords, void* val, size_t elem_size) {
  memcpy((char*)(s->elements) + dense_storage_pos(s, coords) * elem_size, val, elem_size);
}


DENSE_STORAGE* copy_dense_storage(DENSE_STORAGE* rhs, size_t elem_size) {
  DENSE_STORAGE* lhs;
  size_t count = count_dense_storage_elements(rhs), p;
  size_t* shape = malloc(elem_size*rhs->rank);
  if (!shape) return NULL;

  // copy shape array
  for (p = 0; p < rhs->rank; ++p)
    shape[p] = rhs->shape[p];

  //fprintf(stderr, "copy_dense_storage\n");

  lhs = create_dense_storage(elem_size, shape, rhs->rank);

  if (lhs && count) // ensure that allocation worked before copying
    memcpy(lhs->elements, rhs->elements, elem_size * count);

  return lhs;
}


DENSE_STORAGE* create_dense_storage(size_t elem_size, size_t* shape, size_t rank) {
  DENSE_STORAGE* s;
  size_t count;

  if (!(s = malloc(sizeof(DENSE_STORAGE)))) return NULL;
  s->rank       = rank;
  s->shape      = shape;

  //fprintf(stderr, "create_dense_storage: %p\n", s);

  count         = count_dense_storage_elements(s);
  //fprintf(stderr, "count_dense_storage_elements: %d\n", count);

  if (!(s->elements   = malloc(elem_size * count))) {
    free(s->shape);
    free(s);
    s = NULL;
  }

  return s;
}


void delete_dense_storage(DENSE_STORAGE* s) {
  if (s) { // sometimes Ruby passes in NULL storage for some reason (probably on copy construction failure)
    free(s->shape);
    free(s->elements);
    free(s);
  }
}



#endif