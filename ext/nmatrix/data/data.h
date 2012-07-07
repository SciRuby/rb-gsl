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

//#include <stdio.h>

/*
 * Project Includes
 */

#include "types.h"

#ifdef __cplusplus
	// These inlcudes are only needed for C++ programs.
	#include "complex.h"
	#include "rational.h"
	#include "ruby_object.h"
#endif

/*
 * Macros
 */

#define NUM_DTYPES 13

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
		fun<RubyObject>																	\
	}

#define LR_DTYPE_TEMPLATE_TABLE(fun, ret, ...) \
	static ret (*ttable[NUM_DTYPES][NUM_DTYPES])(__VA__ARGS) = {																																																					\
		{fun<unsigned char, unsigned char>, fun<unsigned char, char>, fun<unsigned char, short>, fun<unsigned char, int>, fun<unsigned char, long>,													\
			fun<unsigned char, float>, fun<unsigned char, double>, fun<unsigned char, Complex64>, fun<unsigned char, Complex128>,																							\
			fun<unsigned char, Rational32>, fun<unsigned char, Rational64>, fun<unsigned char, Rational128>, NULL}, 																				\\ unsigned char	\
																																																																																				\
		{fun<char, unsigned char>, fun<char, char>, fun<char, short>, fun<char, int>, fun<char, long>, fun<char, float>, fun<char, double>,																	\
			fun<char, Complex64>, fun<char, Complex128>, fun<char, Rational32>, fun<char, Rational64>, fun<char, Rational128>, NULL},												\\ char						\
																																																																																				\
		{fun<short, unsigned char>, fun<short, char>, fun<short, short>, fun<short, int>, fun<short, long>, fun<short, float>, fun<short, double>,													\
			fun<short, Complex64>, fun<short, Complex128>, fun<short, Rational32>, fun<short, Rational64>, fun<short, Rational128>, NULL},									\\ short					\
																																																																																				\
		{fun<int, unsigned char>, fun<int, char>, fun<int, short>, fun<int, int>, fun<int, long>, fun<int, float>, fun<int, double>,																				\
			fun<int, Complex64>, fun<int, Complex128>, fun<int, Rational32>, fun<int, Rational64>, fun<int, Rational128>, NULL},														\\ int						\
																																																																																				\
		{fun<long, unsigned char>, fun<long, char>, fun<long, short>, fun<long, int>, fun<long, long>, fun<long, float>, fun<long, double>,																	\
			fun<long, Complex64>, fun<long, Complex128>, fun<long, Rational32>, fun<long, Rational64>, fun<long, Rational128>, NULL},												\\ long						\
																																																																																				\
		{fun<float, unsigned char>, fun<float, char>, fun<float, short>, fun<float, int>, fun<float, long>, fun<float, float>, fun<float, double>,													\
			fun<float, Complex64>, fun<float, Complex128>,  NULL, NULL, NULL, NULL},																																				\\ float					\
																																																																																				\
		{fun<double, unsigned char>, fun<double, char>, fun<double, short>, fun<double, int>, fun<double, long>, fun<double, float>, fun<double, double>,										\
			fun<double, Complex64>, fun<double, Complex128>, NULL, NULL, NULL, NULL},																																				\\ double					\
																																																																																				\
		{fun<Complex64, unsigned char>, fun<Complex64, char>, fun<Complex64, short>, fun<Complex64, int>, fun<Complex64, long>,																							\
			fun<Complex64, float>, fun<Complex64, double>, fun<Complex64, Complex64>, fun<Complex64, Complex128>, NULL, NULL, NULL, NULL},									\\ Complex64			\
																																																																																				\
		{fun<Complex128, unsigned char>, fun<Complex128, char>, fun<Complex128, short>, fun<Complex128, int>, fun<Complex128, long>,																				\
			fun<Complex128, float>, fun<Complex128, double>, fun<Complex128, Complex64>, fun<Complex128, Complex128>, NULL, NULL, NULL, NULL},							\\ Complex128			\
																																																																																				\
		{fun<Rational32, unsigned char>, fun<Rational32, char>, fun<Rational32, short>, fun<Rational32, int>, fun<Rational32, long>, NULL, NULL,														\
			NULL, Null, fun<Rational32, Rational32>, fun<Rational32, Rational64>, fun<Rational32, Rational128>, NULL},																			\\ Rational32			\
																																																																																				\
		{fun<Rational64, unsigned char>, fun<Rational64, char>, fun<Rational64, short>, fun<Rational64, int>, fun<Rational64, long>, NULL, NULL,														\
			NULL, Null, fun<Rational64, Rational32>, fun<Rational64, Rational64>, fun<Rational64, Rational128>, NULL},																			\\ Rational64			\
																																																																																				\
		{fun<Rational128, unsigned char>, fun<Rational128, char>, fun<Rational128, short>, fun<Rational128, int>, fun<Rational128, long>, NULL, NULL,												\
			NULL, Null, fun<Rational128, Rational32>, fun<Rational128, Rational64>, fun<Rational128, Rational128>, NULL},																		\\ Rational128		\
																																																																																				\
		{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<RubyObject, RubyObject>}																							\\ RubyObject			\
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

#ifdef __cplusplus

//typedef union {
//  u_int8_t b[2];
//  int16_t s;
//} nm_size16_t;

//typedef union {
//  u_int8_t b[4];
//  int32_t  i;
//  float    f;
//} nm_size32_t;

//typedef union {
//  u_int8_t  b[8];
//  int64_t   q;
//  float     f[2];
//  double    d;
//  Complex64 c;
//} nm_size64_t;

//typedef union {
//  u_int8_t   b[16];
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

#endif

/*
 * Data
 */

#ifdef __cplusplus
extern "C" {
#endif

extern const char* const	DTYPE_NAMES[NUM_DTYPES];
extern const size_t 			DTYPE_SIZES[NUM_DTYPES];

#ifdef __cplusplus
}
#endif

/*
 * Functions
 */

//inline void* test_function(char* op_name, void* ptr, dtype_t left, dtype_t right) {
//	if (ptr == NULL) {
//		// FIXME: Make this do something useful, like raise a Ruby exception.
//		printf("Operation '%s' is not permitted with data types %s and %s.\n", op_name, DTYPE_NAMES[left], DTYPE_NAMES[right]);
//	}
//	
//	return ptr;
//}

#endif
