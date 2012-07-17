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
// This file brings together everything in the storage directory.  It should not
// be included by anything in the storage directory, but should be included by
// files needing to use the storage code.

#ifndef STORAGE_H
#define STORAGE_H

/*
 * Standard Includes
 */

#include <stdlib.h>

/*
 * Project Includes
 */

#include "types.h"

#include "data/data.h"

#include "slice.h"
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

/*
 * Data
 */

extern const char* const STYPE_NAMES[NUM_STYPES];

/*
 * Functions
 */

/////////////////////////
// Copying and Casting //
/////////////////////////

DENSE_STORAGE*	dense_storage_from_list(const LIST_STORAGE* rhs, dtype_t l_dtype);
DENSE_STORAGE*	dense_storage_from_yale(const YALE_STORAGE* rhs, dtype_t l_dtype);
LIST_STORAGE*		list_storage_from_dense(const DENSE_STORAGE* rhs, dtype_t l_dtype);
LIST_STORAGE*		list_storage_from_yale(const YALE_STORAGE* rhs, dtype_t l_dtype);
YALE_STORAGE*		yale_storage_from_list(const LIST_STORAGE* rhs, dtype_t l_dtype);
YALE_STORAGE*		yale_storage_from_dense(const DENSE_STORAGE* rhs, dtype_t l_dtype);

#endif
