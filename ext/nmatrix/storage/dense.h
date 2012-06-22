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
// == dense.h
//
// Dense n-dimensional matrix storage.

#ifndef DENSE_H
#define DENSE_H

// Standard Includes

#include <stdlib.h>

// Project Includes

#include "nmatrix.h"

// Macros

// Types

typedef struct {
	int8_t    dtype;
	size_t    rank;
	size_t*   shape;
	size_t*   offset;
	int       count;
	void*     src;
	void*     elements;
} DENSE_STORAGE;

// Functions

size_t count_dense_storage_elements(const DENSE_STORAGE* s);
bool dense_storage_eqeq(const DENSE_STORAGE* left, const DENSE_STORAGE* right);
bool dense_is_symmetric(const DENSE_STORAGE* mat, int lda, bool hermitian);
size_t dense_storage_pos(DENSE_STORAGE* s, SLICE* slice);
void* dense_storage_get(DENSE_STORAGE* s, SLICE* slice);
void dense_storage_set(DENSE_STORAGE* s, SLICE* slice, void* val);
DENSE_STORAGE* copy_dense_storage(DENSE_STORAGE* rhs);
DENSE_STORAGE* cast_copy_dense_storage(DENSE_STORAGE* rhs, int8_t new_dtype);
DENSE_STORAGE* scast_copy_dense_list(const LIST_STORAGE* rhs, int8_t l_dtype);
DENSE_STORAGE* scast_copy_dense_yale(const YALE_STORAGE* rhs, int8_t l_dtype);
DENSE_STORAGE* create_dense_storage(int8_t dtype, size_t* shape, size_t rank, void* elements, size_t elements_length);
void delete_dense_storage(DENSE_STORAGE* s);
void delete_dense_storage_ref(DENSE_STORAGE* s);
void mark_dense_storage(void* m);

#endif
