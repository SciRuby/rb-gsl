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
// == common.h
//
// Header file for code common to all storage types.

#ifndef STORAGE_COMMON_H
#define STORAGE_COMMON_H

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

/*
 * Macros
 */

/*
 * Types
 */

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

typedef struct {
	// Coordinate of first element
	size_t*	coords;
	// Lenght of slice
	size_t*	lens;
	// 1 - if all lens eql 1
	uint8_t	is_one_el;
} SLICE;

/*
 * Data
 */

/*
 * Functions
 */

/*
 * Calculate the number of elements in the dense storage structure, based on
 * shape and rank.
 */
inline size_t storage_count_max_elements(size_t rank, const size_t* shape) {
  unsigned int i;
  size_t count = 1;
  
  for (i = rank; i-- > 0;) {
  	count *= shape[i];
  }
  
  return count;
}

#endif
