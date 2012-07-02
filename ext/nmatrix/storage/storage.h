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

#include "nmatrix.h"

/*
 * Macros
 */

/*
 * Types
 */

typedef enum {
  DENSE,
  LIST,
  YALE
} stype_t;

typedef struct {
	// Common elements found in all storage types.
	
  dtype_t	dtype;
  size_t	rank;
  size_t*	shape;
  size_t*	offset;
  void*		elements;
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

/* Calculate the max number of elements in the list storage structure, based
 * on shape and rank.
 *
 * FIXME: Massively unsafe.  Should be removed or re-done.
 */
inline size_t storage_count_max_elements(const STORAGE* s) {
  return dense_storage_count_elements((DENSE_STORAGE*)s);
}

#endif
