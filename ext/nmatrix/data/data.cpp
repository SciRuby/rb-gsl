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
// == data.cpp
//
// Functions and data for dealing the data types.

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

#include "types.h"

#include "data.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

extern "C" {

const char* const DTYPE_NAMES[NUM_DTYPES] = {
	"Byte",
	"Int8",
	"Int16",
	"Int32",
	"Int64",
	"Float32",
	"Float64",
	"Complex64",
	"Complex128",
	"Rational32",
	"Rational64",
	"Rational128",
	"RubyObject"
};

const size_t DTYPE_SIZES[NUM_DTYPES] = {
	sizeof(u_int8_t),
	sizeof(int8_t),
	sizeof(int16_t),
	sizeof(int32_t),
	sizeof(int64_t),
	sizeof(float32_t),
	sizeof(float64_t),
	sizeof(Complex64),
	sizeof(Complex128),
	sizeof(Rational32),
	sizeof(Rational64),
	sizeof(Rational128),
	sizeof(RubyObject)
};

}

/*
 * Forward Declarations
 */

/*
 * Functions
 */
