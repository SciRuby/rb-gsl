/////////////////////////////////////////////////////////////////////
// = NMatrix
//
// A linear algebra library for scientific computation in Ruby.
// NMatrix is part of SciRuby.
//
// NMatrix was originally inspired by and derived from NArray, by
// Masahiro Tanaka: http://narray.rubyforge.org
//
// == Copyright Information
//
// SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
// NMatrix is Copyright (c) 2012, Ruby Science Foundation
//
// Please see LICENSE.txt for additional copyright notices.
//
// == Contributing
//
// By contributing source code to SciRuby, you agree to be bound by
// our Contributor Agreement:
//
// * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
//
// == list.h
//
// List-of-lists n-dimensional matrix storage. Uses singly-linked
// lists.

#ifndef LIST_H
#define LIST_H

// Standard Includes

#include <stdlib.h>

// Project Includes

#include "nmatrix.h"
#include "util/sl_list.h"

// Macros

// Types

typedef struct {
	int8_t    dtype;
	size_t    rank;
	size_t*   shape;
	size_t*   offset;
	void*     default_val;
	LIST*     rows;
} LIST_STORAGE;

// Functions

/* Calculate the max number of elements in the list storage structure, based on shape and rank */
inline size_t count_storage_max_elements(const STORAGE* s) {
  return count_dense_storage_elements((DENSE_STORAGE*)s);
}

// Count non-zero elements. See also count_list_storage_nd_elements.
size_t count_list_storage_elements(const LIST_STORAGE* s) {
  return count_list_storage_elements_r(s->rows, s->rank-1);
}

size_t count_list_storage_nd_elements(const LIST_STORAGE* s)
void* list_storage_get(LIST_STORAGE* s, SLICE* slice);
void* list_storage_remove(LIST_STORAGE* s, SLICE* slice);
void* list_storage_insert(LIST_STORAGE* s, SLICE* slice, void* val);
LIST_STORAGE* create_list_storage(int8_t dtype, size_t* shape, size_t rank, void* init_val);
LIST_STORAGE* copy_list_storage(LIST_STORAGE* rhs);
LIST_STORAGE* cast_copy_list_storage(LIST_STORAGE* rhs, int8_t new_dtype);
LIST_STORAGE* scast_copy_list_dense(const DENSE_STORAGE* rhs, int8_t l_dtype);
LIST_STORAGE* scast_copy_list_yale(const YALE_STORAGE* rhs, int8_t l_dtype);
bool list_storage_eqeq(const LIST_STORAGE* left, const LIST_STORAGE* right);
void delete_list_storage(LIST_STORAGE* s);
void mark_list_storage(void* m)

#endif
