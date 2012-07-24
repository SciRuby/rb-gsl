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
// == data.h
//
// Header file for dealing with data types.

#ifndef DATA_TYPES_H
#define DATA_TYPES_H

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

#include "types.h"

#include "complex.h"
#include "rational.h"
#include "ruby_object.h"

/*
 * Macros
 */

#define NUM_DTYPES 13

/*
 * Defines a static array named ttables that hold function pointers to
 * dtype templated versions of the specified function.
 */
#define DTYPE_TEMPLATE_TABLE(fun, ret, ...)					\
	static ret (*ttable[NUM_DTYPES])(__VA_ARGS__) =	{	\
		fun<uint8_t>,																		\
		fun<int8_t>,																		\
		fun<int16_t>,																		\
		fun<int32_t>,																		\
		fun<int64_t>,																		\
		fun<float32_t>,																	\
		fun<float64_t>,																	\
		fun<Complex64>,																	\
		fun<Complex128>,																\
		fun<Rational32>,																\
		fun<Rational64>,																\
		fun<Rational128>,																\
		fun<RubyObject>																	\
	}

/*
 * Same as DTYPE_TEMPLATE_TABLE but for functions that have two template
 * parameters.
 *
 * The left-hand DType is used as the first index, and the right-hand side is
 * the second index.  Not all left- and right-hand side combinations are valid,
 * and an invalid combination will result in a NULL pointer.
 */
#define LR_DTYPE_TEMPLATE_TABLE(fun, ret, ...)																																																														\
	static ret (*ttable[NUM_DTYPES][NUM_DTYPES])(__VA_ARGS__) = {																																																						\
		{fun<uint8_t, uint8_t>, fun<uint8_t, int8_t>, fun<uint8_t, int16_t>, fun<uint8_t, int32_t>, fun<uint8_t, int64_t>, fun<uint8_t, float32_t>, fun<uint8_t, float64_t>,	\
			fun<uint8_t, Complex64>, fun<uint8_t, Complex128>, fun<uint8_t, Rational32>, fun<uint8_t, Rational64>, fun<uint8_t, Rational128>, NULL},														\
																																																																																					\
		{fun<int8_t, uint8_t>, fun<int8_t, int8_t>, fun<int8_t, int16_t>, fun<int8_t, int32_t>, fun<int8_t, int64_t>, fun<int8_t, float32_t>, fun<int8_t, float64_t>,					\
			fun<int8_t, Complex64>, fun<int8_t, Complex128>, fun<int8_t, Rational32>, fun<int8_t, Rational64>, fun<int8_t, Rational128>, NULL},																	\
																																																																																					\
		{fun<int16_t, uint8_t>, fun<int16_t, int8_t>, fun<int16_t, int16_t>, fun<int16_t, int32_t>, fun<int16_t, int64_t>, fun<int16_t, float32_t>, fun<int16_t, float64_t>,	\
			fun<int16_t, Complex64>, fun<int16_t, Complex128>, fun<int16_t, Rational32>, fun<int16_t, Rational64>, fun<int16_t, Rational128>, NULL},														\
																																																																																					\
		{fun<int32_t, uint8_t>, fun<int32_t, int8_t>, fun<int32_t, int16_t>, fun<int32_t, int32_t>, fun<int32_t, int64_t>, fun<int32_t, float32_t>, fun<int32_t, float64_t>,	\
			fun<int32_t, Complex64>, fun<int32_t, Complex128>, fun<int32_t, Rational32>, fun<int32_t, Rational64>, fun<int32_t, Rational128>, NULL},														\
																																																																																					\
		{fun<int64_t, uint8_t>, fun<int64_t, int8_t>, fun<int64_t, int16_t>, fun<int64_t, int32_t>, fun<int64_t, int64_t>, fun<int64_t, float32_t>, fun<int64_t, float64_t>,	\
			fun<int64_t, Complex64>, fun<int64_t, Complex128>, fun<int64_t, Rational32>, fun<int64_t, Rational64>, fun<int64_t, Rational128>, NULL},														\
																																																																																					\
		{fun<float32_t, uint8_t>, fun<float32_t, int8_t>, fun<float32_t, int16_t>, fun<float32_t, int32_t>, fun<float32_t, int64_t>,																					\
			fun<float32_t, float32_t>, fun<float32_t, float64_t>, fun<float32_t, Complex64>, fun<float32_t, Complex128>,  NULL, NULL, NULL, NULL},															\
																																																																																					\
		{fun<float64_t, uint8_t>, fun<float64_t, int8_t>, fun<float64_t, int16_t>, fun<float64_t, int32_t>, fun<float64_t, int64_t>,																					\
			fun<float64_t, float32_t>, fun<float64_t, float64_t>, fun<float64_t, Complex64>, fun<float64_t, Complex128>, NULL, NULL, NULL, NULL},																\
																																																																																					\
		{fun<Complex64, uint8_t>, fun<Complex64, int8_t>, fun<Complex64, int16_t>, fun<Complex64, int32_t>, fun<Complex64, int64_t>,																					\
			fun<Complex64, float32_t>, fun<Complex64, float64_t>, fun<Complex64, Complex64>, fun<Complex64, Complex128>, NULL, NULL, NULL, NULL},																\
																																																																																					\
		{fun<Complex128, uint8_t>, fun<Complex128, int8_t>, fun<Complex128, int16_t>, fun<Complex128, int32_t>, fun<Complex128, int64_t>,																			\
			fun<Complex128, float32_t>, fun<Complex128, float64_t>, fun<Complex128, Complex64>, fun<Complex128, Complex128>, NULL, NULL, NULL, NULL},														\
																																																																																					\
		{fun<Rational32, uint8_t>, fun<Rational32, int8_t>, fun<Rational32, int16_t>, fun<Rational32, int32_t>, fun<Rational32, int64_t>, NULL, NULL,													\
			NULL, NULL, fun<Rational32, Rational32>, fun<Rational32, Rational64>, fun<Rational32, Rational128>, NULL},																													\
																																																																																					\
		{fun<Rational64, uint8_t>, fun<Rational64, int8_t>, fun<Rational64, int16_t>, fun<Rational64, int32_t>, fun<Rational64, int64_t>, NULL, NULL,													\
			NULL, NULL, fun<Rational64, Rational32>, fun<Rational64, Rational64>, fun<Rational64, Rational128>, NULL},																													\
																																																																																					\
		{fun<Rational128, uint8_t>, fun<Rational128, int8_t>, fun<Rational128, int16_t>, fun<Rational128, int32_t>, fun<Rational128, int64_t>, NULL, NULL,										\
			NULL, NULL, fun<Rational128, Rational32>, fun<Rational128, Rational64>, fun<Rational128, Rational128>, NULL},																												\
																																																																																					\
		{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<RubyObject, RubyObject>}																																	\
	}

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
} dtype_t;

//typedef union {
//  uint8_t b[2];
//  int16_t s;
//} nm_size16_t;

//typedef union {
//  uint8_t b[4];
//  int32_t  i;
//  float    f;
//} nm_size32_t;

//typedef union {
//  uint8_t  b[8];
//  int64_t   q;
//  float     f[2];
//  double    d;
//  Complex64 c;
//} nm_size64_t;

//typedef union {
//  uint8_t   b[16];
//  int64_t    i[2];
//  double     d[2];
//  float      f[4];
//  Complex64  c[2];
//  Complex128 z;
//  Rational32 r[4];
//  Rational64 ra[2];
//  Rational128 rat;
//  VALUE      v[2];
//} nm_size128_t;

/*
 * Data
 */

extern const char* const	DTYPE_NAMES[NUM_DTYPES];
extern const size_t 			DTYPE_SIZES[NUM_DTYPES];

/*
 * Functions
 */

RubyObject	rubyobj_from_val(void* val, dtype_t dtype);
void*				rubyobj_to_val(RubyObject obj, dtype_t dtype);

#endif
