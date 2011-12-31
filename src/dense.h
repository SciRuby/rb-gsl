#ifndef DENSE_H
# define DENSE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>
#endif

typedef struct dense_s {
  size_t* shape;
  size_t  rank;
  void*   elements;
} DENSE_STORAGE;

DENSE_STORAGE* create_dense_storage(size_t elem_size, size_t rank, size_t* shape, void* init_val);
void delete_dense_storage(DENSE_STORAGE* s);

size_t count_dense_storage_elements(DENSE_STORAGE* s);

size_t dense_storage_pos(DENSE_STORAGE* s, size_t* coords);
void* dense_storage_get(DENSE_STORAGE* s, size_t* coords, size_t elem_size);
void dense_storage_set(DENSE_STORAGE* s, size_t* coords, void* val, size_t elem_size);

int main();

#endif