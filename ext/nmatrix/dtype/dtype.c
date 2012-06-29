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
// == dtype.c
//
// Functions and data for dealing the data types.

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

#include "dtype.h"

/*
 * Macros
 */

/*
 * Global Variables
 */
 
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

/*
 * Forward Declarations
 */

/*
 * Functions
 */
