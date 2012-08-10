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

struct STORAGE {
	// Common elements found in all storage types.  Must not be re-arranged.
	dtype_t	dtype;
	size_t	rank;
	size_t*	shape;
	size_t*	offset;
	
	virtual void empty(void) = 0;
};

// For binary operations involving matrices that need to be casted.
typedef struct {
	STORAGE* left;
	STORAGE* right;
} STORAGE_PAIR;

typedef struct {
	size_t*	coords; // Coordinate of first element
	size_t*	lengths; // Lengths of slice
	uint8_t	is_one_el; // 1 - if all lens eql 1
} SLICE;

/*
 * Data
 */

/*
 * Functions
 */

size_t storage_count_max_elements(const STORAGE* storage);

#endif // STORAGE_COMMON_H
