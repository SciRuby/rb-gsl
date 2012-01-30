// yale.c - "new yale" storage format for 2D matrices
#ifndef YALE_C
# define YALE_C

#include "nmatrix.h"

YALE_STORAGE* create_yale_storage(int8_t dtype, size_t* shape, size_t rank, size_t init_ndnz) {
  YALE_STORAGE* s;

  if (rank > 2) rb_raise(rb_eNotImpError, "Can only support 2D matrices");

  if (!(s = malloc(sizeof(YALE_STORAGE)))) return NULL;
  s->ndnz  = init_ndnz;
  s->dtype = dtype;
  s->shape = shape;
  s->rank  = rank;

  if (!(s->ija = malloc(sizeof(elem_size) * YALE_CAPACITY(s)))) {
    free(s);
    return NULL;
  } else {
    if (!(s->a = malloc(sizeof(elem_size) * YALE_CAPACITY(s)))) {
      free(s->ija);
      free(s);
      return NULL;
    }
  }

  return s;
}


void* yale_storage_ref(YALE_STORAGE* s, size_t* coords, size_t elem_size) {
  // assume 2D! At some point need to look at how to make an ND yale matrix, maybe?
  if (coords[0] == coords[1]) return YALE_DIAG(s,elem_size,coords[0]);


}

#endif
