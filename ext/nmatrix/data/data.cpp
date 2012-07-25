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

RubyObject rubyobj_from_val(void* val, dtype_t dtype) {
	return RubyObject((VALUE)NULL);
}


// Mostly necessary for copying Ruby Arrays into C arrays, where we want to allocate
// ret_val all at once instead of piecewise.
void rubyobj_to_val_noalloc(RubyObject obj, dtype_t dtype, void* ret_val) {
	
//	switch (dtype) {
//		case BYTE:
//			*(uint8_t*)ret_val	= NUM2UINT(obj);
//			break;
//			
//		case INT8:
//			*(int8_t*)ret_val		= NUM2INT(obj);
//			break;
//		
//		case INT16:
//			*(int16_t*)ret_val		= NUM2INT(obj);
//			break;
//		
//		case INT32:
//			*(int32_t*)ret_val		= NUM2INT(obj);
//			break;
//		
//		case INT64:
//			*(int64_t*)ret_val		= NUM2INT(obj);
//			break;
//		
//		case FLOAT32:
//			*(float32_t*)ret_val = NUM2DBL(obj);
//			break;
//		
//		case FLOAT64:
//			*(float32_t*)ret_val = NUM2DBL(obj);
//			break;
//		
//		case COMPLEX64:
//			*(Complex64*)ret_val
//		case COMPLEX128:
//		case RATIONAL32:
//		case RATIONAL64:
//		case RATIONAL128:
//		case RUBYOBJ:
//			free(ret_val);
//			rb_raise(rb_eException, "Trying to cast a Ruby object to a Ruby object.");
//	}
	
	return ret_val;
}

// Allocate and return a piece of data of the correct dtype, converted from a given RubyObject.
void* rubyobj_to_val(RubyObject obj, dtype_t dtype) {
  void* ret_val = malloc(DTYPE_SIZES[dtype]);

  rubyobj_to_val_noalloc(obj, dtype, ret_val);

  return ret_val;
}
