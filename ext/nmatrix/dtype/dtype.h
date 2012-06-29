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
// == dtype.h
//
// Header file for dealing with data types.

#ifndef DTYPE_H
#define DTYPE_H

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

#ifdef __cplusplus
	// These inlcudes are only needed for C++ programs.
	#include "complex.h"
	#include "object.h"
	#include "rational.h"
#endif

/*
 * Macros
 */

#define NUM_DTYPES 13

/*
 * FIXME: Provide a 2D version of the DTYPE_TEMPLATE_TABLE for operations with
 * different left- and right-hand side data types.
 */ 

#define DTYPE_TEMPLATE_TABLE(fun, ret, ...)					\
	static ret (*ttable[NUM_DTYPES])(__VA_ARGS__) =	{	\
		fun<unsigned char>,															\
		fun<char>,																			\
		fun<short>,																			\
		fun<int>,																				\
		fun<long>,																			\
		fun<float>,																			\
		fun<double>,																		\
		fun<Complex64>,																	\
		fun<Complex128>,																\
		fun<Rational32>,																\
		fun<Rational64>,																\
		fun<Rational128>,																\
		fun<VALUE>																			\
	}

#define EPSILON 1E-10
#define FP_IS_ZERO(n) (-EPSILON < n && n < EPSILON)
#define FP_EQUAL(a, b) (FP_IS_ZERO(a - b))

/*
 * Types
 */

typedef enum {
	BYTE				=  0, // unsigned char
	INT8				=  1, // char
	INT16				=  2, // short
	INT32				=  3, // int
	INT64				=  4, // long
	FLOAT32			=  5, // float
	FLOAT64			=  6, // double
	COMPLEX64		=  7, // Complex64 class
	COMPLEX128	=  8, // Complex128 class
	RATIONAL32	=  9, // Rational32 class
	RATIONAL64	= 10, // Rational64 class
	RATIONAL128	= 11, // Rational128 class
	RUBYOBJ			= 12  // Ruby VALUE type
} dtype;

/*
 * Functions
 */

#endif
