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

/*
 * Standard Includes
 */

#include <stdlib.h>

/*
 * Project Includes
 */

#include "types.h"

#include "data/data.h"

#include "common.h"

#include "util/sl_list.h"

/*
 * Macros
 */

/*
 * Types
 */

struct LIST_STORAGE : STORAGE {
	// List storage specific elements.
	void* default_val;
	LIST* rows;
};

/*
 * Data
 */
 

/*
 * Functions
 */

////////////////
// Lifecycle //
///////////////

LIST_STORAGE*	list_storage_create(dtype_t dtype, size_t* shape, size_t rank, void* init_val);
void					list_storage_delete(STORAGE* s);
void					list_storage_mark(void*);

///////////////
// Accessors //
///////////////

void* list_storage_ref(STORAGE* s, SLICE* slice);
void* list_storage_get(STORAGE* s, SLICE* slice);
void* list_storage_insert(STORAGE* s, SLICE* slice, void* val);
void* list_storage_remove(STORAGE* s, SLICE* slice);

///////////
// Tests //
///////////

bool list_storage_eqeq(const STORAGE* left, const STORAGE* right);

//////////
// Math //
//////////

STORAGE* list_storage_ew_multiply(const STORAGE* left, const STORAGE* right);
STORAGE* list_storage_matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector);

/////////////
// Utility //
/////////////

size_t list_storage_count_elements_r(const LIST* l, size_t recursions);
size_t list_storage_count_nd_elements(const LIST_STORAGE* s);

/*
 * Count non-zero elements. See also count_list_storage_nd_elements.
 */
inline size_t list_storage_count_elements(const LIST_STORAGE* s) {
  return list_storage_count_elements_r(s->rows, s->rank - 1);
}

/////////////////////////
// Copying and Casting //
/////////////////////////

LIST_STORAGE* list_storage_copy(LIST_STORAGE* rhs);
STORAGE*      list_storage_copy_transposed(const STORAGE* rhs_base);
STORAGE* list_storage_cast_copy(const STORAGE* rhs, dtype_t new_dtype);

#endif // LIST_H
