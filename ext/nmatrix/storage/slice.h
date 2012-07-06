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
// == slice.h
//
// Header file for slice related code.

#ifndef SLICE_H
#define SLICE_H

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

#endif
