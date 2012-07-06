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
// == storage.h
//
// Resources for all storage types.

#ifndef STORAGE_H
#define STORAGE_H

/*
 * Standard Includes
 */

#include <stdlib.h>

/*
 * Project Includes
 */

#include "data/data.h"

#include "dense.h"
#include "list.h"
#include "yale.h"

/*
 * Macros
 */

#define NUM_STYPES 3

/*
 * Types
 */

typedef enum {
	DENSE_STORE,
	LIST_STORE,
	YALE_STORE
} stype_t;

typedef struct {
	// Common elements found in all storage types.  Must not be re-arranged.
	dtype_t	dtype;
	size_t	rank;
	size_t*	shape;
	size_t*	offset;
} STORAGE;

// For binary operations involving matrices that need to be casted.
typedef struct {
	STORAGE* left;
	STORAGE* right;
} STORAGE_PAIR;

/*
 * Data
 */

/*
 * Functions
 */

/////////////////////////
// Copying and Casting //
/////////////////////////

DENSE_STORAGE*	dense_storage_from_list(const LIST_STORAGE* rhs, dtype_t l_dtype);
DENSE_STORAGE*	dense_storage_from_yale(const YALE_STORAGE* rhs, dtype_t l_dtype);
LIST_STORAGE*		list_storage_from_dense(const DENSE_STORAGE* rhs, int8_t l_dtype);
LIST_STORAGE*		list_storage_from_yale(const YALE_STORAGE* rhs, int8_t l_dtype);
YALE_STORAGE*		yale_storage_from_list(const LIST_STORAGE* rhs, int8_t l_dtype);
YALE_STORAGE*		yale_storage_from_dense(const DENSE_STORAGE* rhs, int8_t l_dtype);

/* Calculate the max number of elements in the list storage structure, based
 * on shape and rank.
 *
 * FIXME: Massively unsafe.  Should be removed or re-done.
 */
//inline size_t storage_count_max_elements(const STORAGE* s) {
//  return dense_storage_count_elements((DENSE_STORAGE*)s);
//}

#endif
