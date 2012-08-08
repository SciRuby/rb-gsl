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
// == util.cpp
//
// Utility functions and data.

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

#include "types.h"

#include "util.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

/*
 * Forward Declarations
 */


/*
 * Functions
 */

template <typename Type>
Type gcf(Type x, Type y) {
	Type t;
	
	if (x < 0) x = -x;
	if (y < 0) y = -y;
	
	if (x == 0) return y;
	if (y == 0) return x;
	
	while (x > 0) {
		t = x;
		x = y % x;
		y = t;
	}
	
	return y;
}

template <> int16_t
gcf<int16_t>(int16_t, int16_t);
template <> int32_t
gcf<int32_t>(int32_t, int32_t);
template <> int64_t
gcf<int64_t>(int64_t, int64_t);
