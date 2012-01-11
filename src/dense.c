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
  DENSE_STORAGE* lhs = create_dense_storage(elem_size, rhs->shape, rhs->rank);

  if (lhs) // ensure that allocation worked before copying
    memcpy(lhs->elements, rhs->elements, elem_size * count_dense_storage_elements(rhs));

  return lhs;
}


DENSE_STORAGE* create_dense_storage(size_t elem_size, size_t* shape, size_t rank) {
  DENSE_STORAGE* s;
  size_t count, p;

  if (!(s = malloc(sizeof(DENSE_STORAGE)))) return NULL;
  s->rank       = rank;
  s->shape      = shape;
  memcpy(s->shape, shape, sizeof(size_t) * rank);
  count         = count_dense_storage_elements(s);
  if (!(s->elements   = malloc(elem_size * count))) {
    free(s->shape);
    free(s);
  }

  return s;
}


void delete_dense_storage(DENSE_STORAGE* s) {
  free(s->shape);
  free(s->elements);
  free(s);
}


/* int main() {
    size_t shape[] = {3,4,2};
    size_t c0[] = {1,0,0};
    size_t c1[] = {0,0,0};
    void* val = NULL;
    int v2 = 500;
    int init_val = 1;
    DENSE_STORAGE* s = create_dense_storage(sizeof(int), 3, shape, &init_val);

    val = dense_storage_get(s, c1, sizeof(int));
    printf("Got %p: %d\n", val, *((int*)(val)));
    val = dense_storage_get(s, c0, sizeof(int));
    printf("Got %p: %d\n", val, *((int*)(val)));

    dense_storage_set(s, c0, &v2, sizeof(int));
    val = dense_storage_get(s, c0, sizeof(int));
    printf("Got %p: %d\n", val, *((int*)(val)));

    delete_dense_storage(s);

    return 0;
}
*/

#endif