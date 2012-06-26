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
// == yale.h
//
// "new yale" storage format for 2D matrices (like yale, but with
// the diagonal pulled out for O(1) access).
//
// Specifications:
// * dtype and index dtype must necessarily differ
//      * index dtype is defined by whatever unsigned type can store
//        max(rows,cols)
//      * that means vector ija stores only index dtype, but a stores
//        dtype
// * vectors must be able to grow as necessary
//      * maximum size is rows*cols+1

#ifndef DENSE_H
#define DENSE_H

/*
 * Standard Includes
 */

#include <stdlib.h>

/*
 * Project Includes
 */

#include "nmatrix.h"

/*
 * Macros
 */

#define YALE_GROWTH_CONSTANT 1.5

/*
 * Types
 */

typedef struct {
	int8_t    dtype;
	size_t    rank;
	size_t*   shape;
	size_t*   offset;
	void*     a;
	
	// Strictly non-diagonal non-zero count!
	size_t    ndnz;
	
	size_t    capacity;
	int8_t    index_dtype;
	void*     ija;
} YALE_STORAGE;

/*
 * Functions
 */

///////////////
// Lifecycle //
///////////////

YALE_STORAGE* yale_storage_create(int8_t dtype, size_t* shape, size_t rank, size_t init_capacity);
YALE_STORAGE* yale_storage_create_from_old_yale(int8_t dtype, size_t* shape, char* ia, char* ja, char* a, int8_t from_dtype, int8_t from_index_dtype);
YALE_STORAGE*	yale_storage_create_merged(const YALE_STORAGE* template, const YALE_STORAGE* other);
void					yale_storage_delete(YALE_STORAGE* s);
void					yale_storage_init(YALE_STORAGE* s);
void					yale_storage_mark(void* m);

///////////////
// Accessors //
///////////////

void*	yale_storage_ref(YALE_STORAGE* s, SLICE* slice);
char	yale_storage_set(YALE_STORAGE* s, SLICE* slice, void* v);

///////////
// Tests //
///////////

bool yale_storage_eqeq(const YALE_STORAGE* left, const YALE_STORAGE* right);

/////////////
// Utility //
/////////////

int8_t	yale_storage_index_dtype(YALE_STORAGE* s);
void		yale_storage_print_vectors(YALE_STORAGE* s);

int yale_storage_binary_search(YALE_STORAGE* s, y_size_t left, y_size_t right, y_size_t key);

char yale_storage_set_diagonal(YALE_STORAGE* s, y_size_t i, void* v);

char yale_storage_vector_replace(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n);
char yale_storage_vector_insert_resize(YALE_STORAGE* s, y_size_t current_size, y_size_t pos, y_size_t* j, void* val, y_size_t n, bool struct_only);
char yale_storage_vector_insert(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n, bool struct_only);

/////////////////////////
// Copying and Casting //
/////////////////////////

YALE_STORAGE* yale_storage_cast_copy(YALE_STORAGE* rhs, int8_t new_dtype);
YALE_STORAGE* yale_storage_copy(YALE_STORAGE* rhs);
YALE_STORAGE* yale_storage_from_list(const LIST_STORAGE* rhs, int8_t l_dtype);
YALE_STORAGE* yale_storage_from_dense(const DENSE_STORAGE* rhs, int8_t l_dtype);

#endif
