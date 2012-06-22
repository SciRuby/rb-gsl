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

#endif
