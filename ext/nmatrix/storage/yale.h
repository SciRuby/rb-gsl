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

// Standard Includes

#include <stdlib.h>

// Project Includes

#include "nmatrix.h"

// Macros

#define YALE_GROWTH_CONSTANT 1.5

// Types

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

// Functions

void print_vectors(YALE_STORAGE* s);
int8_t yale_index_dtype(YALE_STORAGE* s);
bool yale_storage_eqeq(const YALE_STORAGE* left, const YALE_STORAGE* right);
char yale_vector_replace(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n);
char yale_vector_insert_resize(YALE_STORAGE* s, y_size_t current_size, y_size_t pos, y_size_t* j, void* val, y_size_t n, bool struct_only);
char yale_vector_insert(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n, bool struct_only);
void delete_yale_storage(YALE_STORAGE* s);
void mark_yale_storage(void* m);
YALE_STORAGE* copy_yale_storage(YALE_STORAGE* rhs);
YALE_STORAGE* cast_copy_yale_storage(YALE_STORAGE* rhs, int8_t new_dtype);
YALE_STORAGE* scast_copy_yale_list(const LIST_STORAGE* rhs, int8_t l_dtype);
YALE_STORAGE* scast_copy_yale_dense(const DENSE_STORAGE* rhs, int8_t l_dtype);
YALE_STORAGE* create_merged_yale_storage(const YALE_STORAGE* template, const YALE_STORAGE* other);
YALE_STORAGE* create_yale_storage_from_old_yale(int8_t dtype, size_t* shape, char* ia, char* ja, char* a, int8_t from_dtype, int8_t from_index_dtype);
YALE_STORAGE* create_yale_storage(int8_t dtype, size_t* shape, size_t rank, size_t init_capacity);
void init_yale_storage(YALE_STORAGE* s);
char yale_storage_set_diagonal(YALE_STORAGE* s, y_size_t i, void* v);
int yale_storage_binary_search(YALE_STORAGE* s, y_size_t left, y_size_t right, y_size_t key);
char yale_storage_set(YALE_STORAGE* s, SLICE* slice, void* v);
void* yale_storage_ref(YALE_STORAGE* s, SLICE* slice);

#endif
