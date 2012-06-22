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

/*
 * Lifecycle
 */

DENSE_STORAGE*	dense_storage_create(int8_t dtype, size_t* shape, size_t rank, void* elements, size_t elements_length);
void						dense_storage_delete(DENSE_STORAGE* s);
void						dense_storage_delete_ref(DENSE_STORAGE* s);
void						dense_storage_mark(void* m);

/*
 * Accessors
 */

void*	dense_storage_get(DENSE_STORAGE* s, SLICE* slice);
void	dense_storage_set(DENSE_STORAGE* s, SLICE* slice, void* val);

/*
 * Tests
 */

bool dense_storage_eqeq(const DENSE_STORAGE* left, const DENSE_STORAGE* right);
bool dense_storage_is_symmetric(const DENSE_STORAGE* mat, int lda, bool hermitian);

/*
 * Utility
 */

size_t dense_storage_count_elements(const DENSE_STORAGE* s);
size_t dense_storage_pos(DENSE_STORAGE* s, SLICE* slice);

/*
 * Copying and Casting
 */

DENSE_STORAGE* dense_storage_copy(DENSE_STORAGE* rhs);
DENSE_STORAGE* dense_storage_cast_copy(DENSE_STORAGE* rhs, int8_t new_dtype);
DENSE_STORAGE* dense_storage_from_list(const LIST_STORAGE* rhs, int8_t l_dtype);
DENSE_STORAGE* dense_storage_from_yale(const YALE_STORAGE* rhs, int8_t l_dtype);

#endif
