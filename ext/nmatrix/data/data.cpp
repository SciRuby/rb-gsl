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

#include <ruby.h>

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

const char* const DTYPE_NAMES[NUM_DTYPES] = {
	"byte",
	"int8",
	"int16",
	"int32",
	"int64",
	"float32",
	"float64",
	"complex64",
	"complex128",
	"rational32",
	"rational64",
	"rational128",
	"object"
};

const size_t DTYPE_SIZES[NUM_DTYPES] = {
	sizeof(uint8_t),
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

/*
 * Forward Declarations
 */

/*
 * Functions
 */

/*
 * Converts a RubyObject
 */
void rubyval_to_cval(VALUE val, dtype_t dtype, void* loc) {
	switch (dtype) {
		case BYTE:
			*reinterpret_cast<uint8_t*>(loc)			= static_cast<uint8_t>(RubyObject(val));
			break;
			
		case INT8:
			*reinterpret_cast<int8_t*>(loc)				= static_cast<int8_t>(RubyObject(val));
			break;
			
		case INT16:
			*reinterpret_cast<int16_t*>(loc)			= static_cast<int16_t>(RubyObject(val));
			break;
			
		case INT32:
			*reinterpret_cast<int32_t*>(loc)			= static_cast<int32_t>(RubyObject(val));
			break;
			
		case INT64:
			*reinterpret_cast<int64_t*>(loc)			= static_cast<int64_t>(RubyObject(val));
			break;
			
		case FLOAT32:
			*reinterpret_cast<float32_t*>(loc)		= static_cast<float32_t>(RubyObject(val));
			break;
			
		case FLOAT64:
			*reinterpret_cast<float64_t*>(loc)		= static_cast<float64_t>(RubyObject(val));
			break;
			
		case COMPLEX64:
			*reinterpret_cast<Complex64*>(loc)		= RubyObject(val).to<Complex64>();
			break;
			
		case COMPLEX128:
			*reinterpret_cast<Complex128*>(loc)		= RubyObject(val).to<Complex64>();
			break;
			
		case RATIONAL32:
			*reinterpret_cast<Rational32*>(loc)		= RubyObject(val).to<Rational32>();
			break;
			
		case RATIONAL64:
			*reinterpret_cast<Rational64*>(loc)		= RubyObject(val).to<Rational64>();
			break;
			
		case RATIONAL128:
			*reinterpret_cast<Rational128*>(loc)	= RubyObject(val).to<Rational128>();
			break;
			
		case RUBYOBJ:
			rb_raise(rb_eTypeError, "Attempting a bad conversion from a Ruby value.");
			break;
	}
}

/*
 * Documentation goes here.
 *
 * FIXME: The actual constructors still need to be defined.
 */
RubyObject rubyobj_from_cval(void* val, dtype_t dtype) {
	switch (dtype) {
		case BYTE:
			return RubyObject(static_cast<uint8_t>(*val));
			
		case INT8:
			return RubyObject(static_cast<int8_t>(*val));
			
		case INT16:
			return RubyObject(static_cast<int16_t>(*val));
			
		case INT32:
			return RubyObject(static_cast<int32_t>(*val));
			
		case INT64:
			return RubyObject(static_cast<int64_t>(*val));
			
		case FLOAT32:
			return RubyObject(static_cast<float32_t>(*val));
			
		case FLOAT64:
			return RubyObject(static_cast<float64_t>(*val));
			
		case COMPLEX64:
			return RubyObject(static_cast<Complex64>(*val));
			
		case COMPLEX128:
			return RubyObject(static_cast<Complex128>(*val));
			
		case RATIONAL32:
			return RubyObject(static_cast<Rational32>(*val));
			
		case RATIONAL64:
			return RubyObject(static_cast<Rational64>(*val));
			
		case RATIONAL128:
			return RubyObject(static_cast<Rational128>(*val));
			
		case RUBYOBJ:
			rb_raise(rb_eTypeError, "Attempting a bad conversion from a Ruby value.");
	}
}


/*
 * Allocate and return a piece of data of the correct dtype, converted from a
 * given RubyObject.
 */
void* rubyobj_to_cval(VALUE val, dtype_t dtype) {
  void* ret_val = malloc(DTYPE_SIZES[dtype]);

  rubyobj_to_cval_noalloc(val, dtype, ret_val);

  return ret_val;
}

