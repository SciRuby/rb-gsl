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

#ifndef DATA_H
#define DATA_H

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
#define NUM_ITYPES 4



#define NAMED_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...)					\
	static ret (*(name)[NUM_DTYPES])(__VA_ARGS__) =	{	\
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
	};


/*
 * Defines a static array named ttable that hold function pointers to
 * dtype templated versions of the specified function.
 */
#define DTYPE_TEMPLATE_TABLE(fun, ret, ...)					NAMED_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)

#define NAMED_ITYPE_TEMPLATE_TABLE(name, fun, ret, ...) \
  static ret (*(name)[NUM_ITYPES])(__VA_ARGS__) = { \
    fun<uint8_t>, \
    fun<uint16_t>,  \
    fun<uint32_t>,  \
    fun<uint64_t>  \
  };

#define ITYPE_TEMPLATE_TABLE(fun, ret, ...)   NAMED_ITYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)

#define STYPE_MARK_TABLE(name) \
  static void (*(name)[NUM_STYPES])(void*) = {  \
    dense_storage_mark, \
    list_storage_mark,  \
    yale_storage_mark \
  };

#define STYPE_CAST_COPY_TABLE(name)                                                   \
  static STORAGE* (*(name)[NUM_STYPES][NUM_STYPES])(const STORAGE*, dtype_t) = {      \
    { dense_storage_cast_copy,  dense_storage_from_list,  dense_storage_from_yale },  \
    { list_storage_from_dense,  list_storage_cast_copy,   list_storage_from_yale  },  \
    { yale_storage_from_dense,  yale_storage_from_list,   yale_storage_cast_copy  }   \
  };

// First template argument is the temp variable type (e.g., long int for T = int)
#define BLAS_TEMPLATE_TABLE(table_name, fun, ret, ...) \
  static ret (*(table_name)[NUM_DTYPES])(__VA_ARGS__) = { \
    fun<int16_t,uint8_t>, \
    fun<int16_t, int8_t>, \
    fun<int32_t,int16_t>, \
    fun<int64_t,int32_t>, \
    fun<int64_t,int64_t>, \
    fun<float64_t,float32_t>, \
    fun<float64_t,float64_t>, \
    fun<Complex128,Complex64>,  \
    fun<Complex128,Complex128>, \
    fun<Rational128,Rational32>,  \
    fun<Rational128,Rational64>,  \
    fun<Rational128,Rational128>  \
  }


/*
 * Same as DTYPE_TEMPLATE_TABLE but for functions that have two template
 * parameters.
 *
 * The left-hand DType is used as the first index, and the right-hand side is
 * the second index.  Not all left- and right-hand side combinations are valid,
 * and an invalid combination will result in a NULL pointer.
 */
#define NAMED_LR_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...)																																																														\
	static ret (*(name)[NUM_DTYPES][NUM_DTYPES])(__VA_ARGS__) = {																																																						\
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
	};


#define LR_DTYPE_TEMPLATE_TABLE(fun, ret, ...)    NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)
	
#define NAMED_LRI_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...) \
static ret (*(name)[NUM_DTYPES][NUM_DTYPES][NUM_ITYPES])(__VA_ARGS__) = { \
  {{fun<uint8_t,uint8_t,uint8_t>,fun<uint8_t,uint8_t,uint16_t>,fun<uint8_t,uint8_t,uint32_t>,fun<uint8_t,uint8_t,uint64_t>},{fun<uint8_t,int8_t,uint8_t>,fun<uint8_t,int8_t,uint16_t>,fun<uint8_t,int8_t,uint32_t>,fun<uint8_t,int8_t,uint64_t>},{fun<uint8_t,int16_t,uint8_t>,fun<uint8_t,int16_t,uint16_t>,fun<uint8_t,int16_t,uint32_t>,fun<uint8_t,int16_t,uint64_t>},{fun<uint8_t,int32_t,uint8_t>,fun<uint8_t,int32_t,uint16_t>,fun<uint8_t,int32_t,uint32_t>,fun<uint8_t,int32_t,uint64_t>},{fun<uint8_t,int64_t,uint8_t>,fun<uint8_t,int64_t,uint16_t>,fun<uint8_t,int64_t,uint32_t>,fun<uint8_t,int64_t,uint64_t>},{fun<uint8_t,float32_t,uint8_t>,fun<uint8_t,float32_t,uint16_t>,fun<uint8_t,float32_t,uint32_t>,fun<uint8_t,float32_t,uint64_t>},{fun<uint8_t,float64_t,uint8_t>,fun<uint8_t,float64_t,uint16_t>,fun<uint8_t,float64_t,uint32_t>,fun<uint8_t,float64_t,uint64_t>},{fun<uint8_t,Complex64,uint8_t>,fun<uint8_t,Complex64,uint16_t>,fun<uint8_t,Complex64,uint32_t>,fun<uint8_t,Complex64,uint64_t>},{fun<uint8_t,Complex128,uint8_t>,fun<uint8_t,Complex128,uint16_t>,fun<uint8_t,Complex128,uint32_t>,fun<uint8_t,Complex128,uint64_t>},{fun<uint8_t,Rational32,uint8_t>,fun<uint8_t,Rational32,uint16_t>,fun<uint8_t,Rational32,uint32_t>,fun<uint8_t,Rational32,uint64_t>},{fun<uint8_t,Rational64,uint8_t>,fun<uint8_t,Rational64,uint16_t>,fun<uint8_t,Rational64,uint32_t>,fun<uint8_t,Rational64,uint64_t>},{fun<uint8_t,Rational128,uint8_t>,fun<uint8_t,Rational128,uint16_t>,fun<uint8_t,Rational128,uint32_t>,fun<uint8_t,Rational128,uint64_t>},{fun<uint8_t,RubyObject,uint8_t>,fun<uint8_t,RubyObject,uint16_t>,fun<uint8_t,RubyObject,uint32_t>,fun<uint8_t,RubyObject,uint64_t>}}, \
  {{fun<uint8_t,uint8_t,uint8_t>,fun<uint8_t,uint8_t,uint16_t>,fun<uint8_t,uint8_t,uint32_t>,fun<uint8_t,uint8_t,uint64_t>},{fun<uint8_t,int8_t,uint8_t>,fun<uint8_t,int8_t,uint16_t>,fun<uint8_t,int8_t,uint32_t>,fun<uint8_t,int8_t,uint64_t>},{fun<uint8_t,int16_t,uint8_t>,fun<uint8_t,int16_t,uint16_t>,fun<uint8_t,int16_t,uint32_t>,fun<uint8_t,int16_t,uint64_t>},{fun<uint8_t,int32_t,uint8_t>,fun<uint8_t,int32_t,uint16_t>,fun<uint8_t,int32_t,uint32_t>,fun<uint8_t,int32_t,uint64_t>},{fun<uint8_t,int64_t,uint8_t>,fun<uint8_t,int64_t,uint16_t>,fun<uint8_t,int64_t,uint32_t>,fun<uint8_t,int64_t,uint64_t>},{fun<uint8_t,float32_t,uint8_t>,fun<uint8_t,float32_t,uint16_t>,fun<uint8_t,float32_t,uint32_t>,fun<uint8_t,float32_t,uint64_t>},{fun<uint8_t,float64_t,uint8_t>,fun<uint8_t,float64_t,uint16_t>,fun<uint8_t,float64_t,uint32_t>,fun<uint8_t,float64_t,uint64_t>},{fun<uint8_t,Complex64,uint8_t>,fun<uint8_t,Complex64,uint16_t>,fun<uint8_t,Complex64,uint32_t>,fun<uint8_t,Complex64,uint64_t>},{fun<uint8_t,Complex128,uint8_t>,fun<uint8_t,Complex128,uint16_t>,fun<uint8_t,Complex128,uint32_t>,fun<uint8_t,Complex128,uint64_t>},{fun<uint8_t,Rational32,uint8_t>,fun<uint8_t,Rational32,uint16_t>,fun<uint8_t,Rational32,uint32_t>,fun<uint8_t,Rational32,uint64_t>},{fun<uint8_t,Rational64,uint8_t>,fun<uint8_t,Rational64,uint16_t>,fun<uint8_t,Rational64,uint32_t>,fun<uint8_t,Rational64,uint64_t>},{fun<uint8_t,Rational128,uint8_t>,fun<uint8_t,Rational128,uint16_t>,fun<uint8_t,Rational128,uint32_t>,fun<uint8_t,Rational128,uint64_t>},{fun<uint8_t,RubyObject,uint8_t>,fun<uint8_t,RubyObject,uint16_t>,fun<uint8_t,RubyObject,uint32_t>,fun<uint8_t,RubyObject,uint64_t>}}, \
  {{fun<int8_t,uint8_t,uint8_t>,fun<int8_t,uint8_t,uint16_t>,fun<int8_t,uint8_t,uint32_t>,fun<int8_t,uint8_t,uint64_t>},{fun<int8_t,int8_t,uint8_t>,fun<int8_t,int8_t,uint16_t>,fun<int8_t,int8_t,uint32_t>,fun<int8_t,int8_t,uint64_t>},{fun<int8_t,int16_t,uint8_t>,fun<int8_t,int16_t,uint16_t>,fun<int8_t,int16_t,uint32_t>,fun<int8_t,int16_t,uint64_t>},{fun<int8_t,int32_t,uint8_t>,fun<int8_t,int32_t,uint16_t>,fun<int8_t,int32_t,uint32_t>,fun<int8_t,int32_t,uint64_t>},{fun<int8_t,int64_t,uint8_t>,fun<int8_t,int64_t,uint16_t>,fun<int8_t,int64_t,uint32_t>,fun<int8_t,int64_t,uint64_t>},{fun<int8_t,float32_t,uint8_t>,fun<int8_t,float32_t,uint16_t>,fun<int8_t,float32_t,uint32_t>,fun<int8_t,float32_t,uint64_t>},{fun<int8_t,float64_t,uint8_t>,fun<int8_t,float64_t,uint16_t>,fun<int8_t,float64_t,uint32_t>,fun<int8_t,float64_t,uint64_t>},{fun<int8_t,Complex64,uint8_t>,fun<int8_t,Complex64,uint16_t>,fun<int8_t,Complex64,uint32_t>,fun<int8_t,Complex64,uint64_t>},{fun<int8_t,Complex128,uint8_t>,fun<int8_t,Complex128,uint16_t>,fun<int8_t,Complex128,uint32_t>,fun<int8_t,Complex128,uint64_t>},{fun<int8_t,Rational32,uint8_t>,fun<int8_t,Rational32,uint16_t>,fun<int8_t,Rational32,uint32_t>,fun<int8_t,Rational32,uint64_t>},{fun<int8_t,Rational64,uint8_t>,fun<int8_t,Rational64,uint16_t>,fun<int8_t,Rational64,uint32_t>,fun<int8_t,Rational64,uint64_t>},{fun<int8_t,Rational128,uint8_t>,fun<int8_t,Rational128,uint16_t>,fun<int8_t,Rational128,uint32_t>,fun<int8_t,Rational128,uint64_t>},{fun<int8_t,RubyObject,uint8_t>,fun<int8_t,RubyObject,uint16_t>,fun<int8_t,RubyObject,uint32_t>,fun<int8_t,RubyObject,uint64_t>}}, \
  {{fun<int8_t,uint8_t,uint8_t>,fun<int8_t,uint8_t,uint16_t>,fun<int8_t,uint8_t,uint32_t>,fun<int8_t,uint8_t,uint64_t>},{fun<int8_t,int8_t,uint8_t>,fun<int8_t,int8_t,uint16_t>,fun<int8_t,int8_t,uint32_t>,fun<int8_t,int8_t,uint64_t>},{fun<int8_t,int16_t,uint8_t>,fun<int8_t,int16_t,uint16_t>,fun<int8_t,int16_t,uint32_t>,fun<int8_t,int16_t,uint64_t>},{fun<int8_t,int32_t,uint8_t>,fun<int8_t,int32_t,uint16_t>,fun<int8_t,int32_t,uint32_t>,fun<int8_t,int32_t,uint64_t>},{fun<int8_t,int64_t,uint8_t>,fun<int8_t,int64_t,uint16_t>,fun<int8_t,int64_t,uint32_t>,fun<int8_t,int64_t,uint64_t>},{fun<int8_t,float32_t,uint8_t>,fun<int8_t,float32_t,uint16_t>,fun<int8_t,float32_t,uint32_t>,fun<int8_t,float32_t,uint64_t>},{fun<int8_t,float64_t,uint8_t>,fun<int8_t,float64_t,uint16_t>,fun<int8_t,float64_t,uint32_t>,fun<int8_t,float64_t,uint64_t>},{fun<int8_t,Complex64,uint8_t>,fun<int8_t,Complex64,uint16_t>,fun<int8_t,Complex64,uint32_t>,fun<int8_t,Complex64,uint64_t>},{fun<int8_t,Complex128,uint8_t>,fun<int8_t,Complex128,uint16_t>,fun<int8_t,Complex128,uint32_t>,fun<int8_t,Complex128,uint64_t>},{fun<int8_t,Rational32,uint8_t>,fun<int8_t,Rational32,uint16_t>,fun<int8_t,Rational32,uint32_t>,fun<int8_t,Rational32,uint64_t>},{fun<int8_t,Rational64,uint8_t>,fun<int8_t,Rational64,uint16_t>,fun<int8_t,Rational64,uint32_t>,fun<int8_t,Rational64,uint64_t>},{fun<int8_t,Rational128,uint8_t>,fun<int8_t,Rational128,uint16_t>,fun<int8_t,Rational128,uint32_t>,fun<int8_t,Rational128,uint64_t>},{fun<int8_t,RubyObject,uint8_t>,fun<int8_t,RubyObject,uint16_t>,fun<int8_t,RubyObject,uint32_t>,fun<int8_t,RubyObject,uint64_t>}}, \
  {{fun<int16_t,uint8_t,uint8_t>,fun<int16_t,uint8_t,uint16_t>,fun<int16_t,uint8_t,uint32_t>,fun<int16_t,uint8_t,uint64_t>},{fun<int16_t,int8_t,uint8_t>,fun<int16_t,int8_t,uint16_t>,fun<int16_t,int8_t,uint32_t>,fun<int16_t,int8_t,uint64_t>},{fun<int16_t,int16_t,uint8_t>,fun<int16_t,int16_t,uint16_t>,fun<int16_t,int16_t,uint32_t>,fun<int16_t,int16_t,uint64_t>},{fun<int16_t,int32_t,uint8_t>,fun<int16_t,int32_t,uint16_t>,fun<int16_t,int32_t,uint32_t>,fun<int16_t,int32_t,uint64_t>},{fun<int16_t,int64_t,uint8_t>,fun<int16_t,int64_t,uint16_t>,fun<int16_t,int64_t,uint32_t>,fun<int16_t,int64_t,uint64_t>},{fun<int16_t,float32_t,uint8_t>,fun<int16_t,float32_t,uint16_t>,fun<int16_t,float32_t,uint32_t>,fun<int16_t,float32_t,uint64_t>},{fun<int16_t,float64_t,uint8_t>,fun<int16_t,float64_t,uint16_t>,fun<int16_t,float64_t,uint32_t>,fun<int16_t,float64_t,uint64_t>},{fun<int16_t,Complex64,uint8_t>,fun<int16_t,Complex64,uint16_t>,fun<int16_t,Complex64,uint32_t>,fun<int16_t,Complex64,uint64_t>},{fun<int16_t,Complex128,uint8_t>,fun<int16_t,Complex128,uint16_t>,fun<int16_t,Complex128,uint32_t>,fun<int16_t,Complex128,uint64_t>},{fun<int16_t,Rational32,uint8_t>,fun<int16_t,Rational32,uint16_t>,fun<int16_t,Rational32,uint32_t>,fun<int16_t,Rational32,uint64_t>},{fun<int16_t,Rational64,uint8_t>,fun<int16_t,Rational64,uint16_t>,fun<int16_t,Rational64,uint32_t>,fun<int16_t,Rational64,uint64_t>},{fun<int16_t,Rational128,uint8_t>,fun<int16_t,Rational128,uint16_t>,fun<int16_t,Rational128,uint32_t>,fun<int16_t,Rational128,uint64_t>},{fun<int16_t,RubyObject,uint8_t>,fun<int16_t,RubyObject,uint16_t>,fun<int16_t,RubyObject,uint32_t>,fun<int16_t,RubyObject,uint64_t>}}, \
  {{fun<int16_t,uint8_t,uint8_t>,fun<int16_t,uint8_t,uint16_t>,fun<int16_t,uint8_t,uint32_t>,fun<int16_t,uint8_t,uint64_t>},{fun<int16_t,int8_t,uint8_t>,fun<int16_t,int8_t,uint16_t>,fun<int16_t,int8_t,uint32_t>,fun<int16_t,int8_t,uint64_t>},{fun<int16_t,int16_t,uint8_t>,fun<int16_t,int16_t,uint16_t>,fun<int16_t,int16_t,uint32_t>,fun<int16_t,int16_t,uint64_t>},{fun<int16_t,int32_t,uint8_t>,fun<int16_t,int32_t,uint16_t>,fun<int16_t,int32_t,uint32_t>,fun<int16_t,int32_t,uint64_t>},{fun<int16_t,int64_t,uint8_t>,fun<int16_t,int64_t,uint16_t>,fun<int16_t,int64_t,uint32_t>,fun<int16_t,int64_t,uint64_t>},{fun<int16_t,float32_t,uint8_t>,fun<int16_t,float32_t,uint16_t>,fun<int16_t,float32_t,uint32_t>,fun<int16_t,float32_t,uint64_t>},{fun<int16_t,float64_t,uint8_t>,fun<int16_t,float64_t,uint16_t>,fun<int16_t,float64_t,uint32_t>,fun<int16_t,float64_t,uint64_t>},{fun<int16_t,Complex64,uint8_t>,fun<int16_t,Complex64,uint16_t>,fun<int16_t,Complex64,uint32_t>,fun<int16_t,Complex64,uint64_t>},{fun<int16_t,Complex128,uint8_t>,fun<int16_t,Complex128,uint16_t>,fun<int16_t,Complex128,uint32_t>,fun<int16_t,Complex128,uint64_t>},{fun<int16_t,Rational32,uint8_t>,fun<int16_t,Rational32,uint16_t>,fun<int16_t,Rational32,uint32_t>,fun<int16_t,Rational32,uint64_t>},{fun<int16_t,Rational64,uint8_t>,fun<int16_t,Rational64,uint16_t>,fun<int16_t,Rational64,uint32_t>,fun<int16_t,Rational64,uint64_t>},{fun<int16_t,Rational128,uint8_t>,fun<int16_t,Rational128,uint16_t>,fun<int16_t,Rational128,uint32_t>,fun<int16_t,Rational128,uint64_t>},{fun<int16_t,RubyObject,uint8_t>,fun<int16_t,RubyObject,uint16_t>,fun<int16_t,RubyObject,uint32_t>,fun<int16_t,RubyObject,uint64_t>}}, \
  {{fun<int32_t,uint8_t,uint8_t>,fun<int32_t,uint8_t,uint16_t>,fun<int32_t,uint8_t,uint32_t>,fun<int32_t,uint8_t,uint64_t>},{fun<int32_t,int8_t,uint8_t>,fun<int32_t,int8_t,uint16_t>,fun<int32_t,int8_t,uint32_t>,fun<int32_t,int8_t,uint64_t>},{fun<int32_t,int16_t,uint8_t>,fun<int32_t,int16_t,uint16_t>,fun<int32_t,int16_t,uint32_t>,fun<int32_t,int16_t,uint64_t>},{fun<int32_t,int32_t,uint8_t>,fun<int32_t,int32_t,uint16_t>,fun<int32_t,int32_t,uint32_t>,fun<int32_t,int32_t,uint64_t>},{fun<int32_t,int64_t,uint8_t>,fun<int32_t,int64_t,uint16_t>,fun<int32_t,int64_t,uint32_t>,fun<int32_t,int64_t,uint64_t>},{fun<int32_t,float32_t,uint8_t>,fun<int32_t,float32_t,uint16_t>,fun<int32_t,float32_t,uint32_t>,fun<int32_t,float32_t,uint64_t>},{fun<int32_t,float64_t,uint8_t>,fun<int32_t,float64_t,uint16_t>,fun<int32_t,float64_t,uint32_t>,fun<int32_t,float64_t,uint64_t>},{fun<int32_t,Complex64,uint8_t>,fun<int32_t,Complex64,uint16_t>,fun<int32_t,Complex64,uint32_t>,fun<int32_t,Complex64,uint64_t>},{fun<int32_t,Complex128,uint8_t>,fun<int32_t,Complex128,uint16_t>,fun<int32_t,Complex128,uint32_t>,fun<int32_t,Complex128,uint64_t>},{fun<int32_t,Rational32,uint8_t>,fun<int32_t,Rational32,uint16_t>,fun<int32_t,Rational32,uint32_t>,fun<int32_t,Rational32,uint64_t>},{fun<int32_t,Rational64,uint8_t>,fun<int32_t,Rational64,uint16_t>,fun<int32_t,Rational64,uint32_t>,fun<int32_t,Rational64,uint64_t>},{fun<int32_t,Rational128,uint8_t>,fun<int32_t,Rational128,uint16_t>,fun<int32_t,Rational128,uint32_t>,fun<int32_t,Rational128,uint64_t>},{fun<int32_t,RubyObject,uint8_t>,fun<int32_t,RubyObject,uint16_t>,fun<int32_t,RubyObject,uint32_t>,fun<int32_t,RubyObject,uint64_t>}}, \
  {{fun<int32_t,uint8_t,uint8_t>,fun<int32_t,uint8_t,uint16_t>,fun<int32_t,uint8_t,uint32_t>,fun<int32_t,uint8_t,uint64_t>},{fun<int32_t,int8_t,uint8_t>,fun<int32_t,int8_t,uint16_t>,fun<int32_t,int8_t,uint32_t>,fun<int32_t,int8_t,uint64_t>},{fun<int32_t,int16_t,uint8_t>,fun<int32_t,int16_t,uint16_t>,fun<int32_t,int16_t,uint32_t>,fun<int32_t,int16_t,uint64_t>},{fun<int32_t,int32_t,uint8_t>,fun<int32_t,int32_t,uint16_t>,fun<int32_t,int32_t,uint32_t>,fun<int32_t,int32_t,uint64_t>},{fun<int32_t,int64_t,uint8_t>,fun<int32_t,int64_t,uint16_t>,fun<int32_t,int64_t,uint32_t>,fun<int32_t,int64_t,uint64_t>},{fun<int32_t,float32_t,uint8_t>,fun<int32_t,float32_t,uint16_t>,fun<int32_t,float32_t,uint32_t>,fun<int32_t,float32_t,uint64_t>},{fun<int32_t,float64_t,uint8_t>,fun<int32_t,float64_t,uint16_t>,fun<int32_t,float64_t,uint32_t>,fun<int32_t,float64_t,uint64_t>},{fun<int32_t,Complex64,uint8_t>,fun<int32_t,Complex64,uint16_t>,fun<int32_t,Complex64,uint32_t>,fun<int32_t,Complex64,uint64_t>},{fun<int32_t,Complex128,uint8_t>,fun<int32_t,Complex128,uint16_t>,fun<int32_t,Complex128,uint32_t>,fun<int32_t,Complex128,uint64_t>},{fun<int32_t,Rational32,uint8_t>,fun<int32_t,Rational32,uint16_t>,fun<int32_t,Rational32,uint32_t>,fun<int32_t,Rational32,uint64_t>},{fun<int32_t,Rational64,uint8_t>,fun<int32_t,Rational64,uint16_t>,fun<int32_t,Rational64,uint32_t>,fun<int32_t,Rational64,uint64_t>},{fun<int32_t,Rational128,uint8_t>,fun<int32_t,Rational128,uint16_t>,fun<int32_t,Rational128,uint32_t>,fun<int32_t,Rational128,uint64_t>},{fun<int32_t,RubyObject,uint8_t>,fun<int32_t,RubyObject,uint16_t>,fun<int32_t,RubyObject,uint32_t>,fun<int32_t,RubyObject,uint64_t>}}, \
  {{fun<int64_t,uint8_t,uint8_t>,fun<int64_t,uint8_t,uint16_t>,fun<int64_t,uint8_t,uint32_t>,fun<int64_t,uint8_t,uint64_t>},{fun<int64_t,int8_t,uint8_t>,fun<int64_t,int8_t,uint16_t>,fun<int64_t,int8_t,uint32_t>,fun<int64_t,int8_t,uint64_t>},{fun<int64_t,int16_t,uint8_t>,fun<int64_t,int16_t,uint16_t>,fun<int64_t,int16_t,uint32_t>,fun<int64_t,int16_t,uint64_t>},{fun<int64_t,int32_t,uint8_t>,fun<int64_t,int32_t,uint16_t>,fun<int64_t,int32_t,uint32_t>,fun<int64_t,int32_t,uint64_t>},{fun<int64_t,int64_t,uint8_t>,fun<int64_t,int64_t,uint16_t>,fun<int64_t,int64_t,uint32_t>,fun<int64_t,int64_t,uint64_t>},{fun<int64_t,float32_t,uint8_t>,fun<int64_t,float32_t,uint16_t>,fun<int64_t,float32_t,uint32_t>,fun<int64_t,float32_t,uint64_t>},{fun<int64_t,float64_t,uint8_t>,fun<int64_t,float64_t,uint16_t>,fun<int64_t,float64_t,uint32_t>,fun<int64_t,float64_t,uint64_t>},{fun<int64_t,Complex64,uint8_t>,fun<int64_t,Complex64,uint16_t>,fun<int64_t,Complex64,uint32_t>,fun<int64_t,Complex64,uint64_t>},{fun<int64_t,Complex128,uint8_t>,fun<int64_t,Complex128,uint16_t>,fun<int64_t,Complex128,uint32_t>,fun<int64_t,Complex128,uint64_t>},{fun<int64_t,Rational32,uint8_t>,fun<int64_t,Rational32,uint16_t>,fun<int64_t,Rational32,uint32_t>,fun<int64_t,Rational32,uint64_t>},{fun<int64_t,Rational64,uint8_t>,fun<int64_t,Rational64,uint16_t>,fun<int64_t,Rational64,uint32_t>,fun<int64_t,Rational64,uint64_t>},{fun<int64_t,Rational128,uint8_t>,fun<int64_t,Rational128,uint16_t>,fun<int64_t,Rational128,uint32_t>,fun<int64_t,Rational128,uint64_t>},{fun<int64_t,RubyObject,uint8_t>,fun<int64_t,RubyObject,uint16_t>,fun<int64_t,RubyObject,uint32_t>,fun<int64_t,RubyObject,uint64_t>}}, \
  {{fun<int64_t,uint8_t,uint8_t>,fun<int64_t,uint8_t,uint16_t>,fun<int64_t,uint8_t,uint32_t>,fun<int64_t,uint8_t,uint64_t>},{fun<int64_t,int8_t,uint8_t>,fun<int64_t,int8_t,uint16_t>,fun<int64_t,int8_t,uint32_t>,fun<int64_t,int8_t,uint64_t>},{fun<int64_t,int16_t,uint8_t>,fun<int64_t,int16_t,uint16_t>,fun<int64_t,int16_t,uint32_t>,fun<int64_t,int16_t,uint64_t>},{fun<int64_t,int32_t,uint8_t>,fun<int64_t,int32_t,uint16_t>,fun<int64_t,int32_t,uint32_t>,fun<int64_t,int32_t,uint64_t>},{fun<int64_t,int64_t,uint8_t>,fun<int64_t,int64_t,uint16_t>,fun<int64_t,int64_t,uint32_t>,fun<int64_t,int64_t,uint64_t>},{fun<int64_t,float32_t,uint8_t>,fun<int64_t,float32_t,uint16_t>,fun<int64_t,float32_t,uint32_t>,fun<int64_t,float32_t,uint64_t>},{fun<int64_t,float64_t,uint8_t>,fun<int64_t,float64_t,uint16_t>,fun<int64_t,float64_t,uint32_t>,fun<int64_t,float64_t,uint64_t>},{fun<int64_t,Complex64,uint8_t>,fun<int64_t,Complex64,uint16_t>,fun<int64_t,Complex64,uint32_t>,fun<int64_t,Complex64,uint64_t>},{fun<int64_t,Complex128,uint8_t>,fun<int64_t,Complex128,uint16_t>,fun<int64_t,Complex128,uint32_t>,fun<int64_t,Complex128,uint64_t>},{fun<int64_t,Rational32,uint8_t>,fun<int64_t,Rational32,uint16_t>,fun<int64_t,Rational32,uint32_t>,fun<int64_t,Rational32,uint64_t>},{fun<int64_t,Rational64,uint8_t>,fun<int64_t,Rational64,uint16_t>,fun<int64_t,Rational64,uint32_t>,fun<int64_t,Rational64,uint64_t>},{fun<int64_t,Rational128,uint8_t>,fun<int64_t,Rational128,uint16_t>,fun<int64_t,Rational128,uint32_t>,fun<int64_t,Rational128,uint64_t>},{fun<int64_t,RubyObject,uint8_t>,fun<int64_t,RubyObject,uint16_t>,fun<int64_t,RubyObject,uint32_t>,fun<int64_t,RubyObject,uint64_t>}}, \
  {{fun<float32_t,uint8_t,uint8_t>,fun<float32_t,uint8_t,uint16_t>,fun<float32_t,uint8_t,uint32_t>,fun<float32_t,uint8_t,uint64_t>},{fun<float32_t,int8_t,uint8_t>,fun<float32_t,int8_t,uint16_t>,fun<float32_t,int8_t,uint32_t>,fun<float32_t,int8_t,uint64_t>},{fun<float32_t,int16_t,uint8_t>,fun<float32_t,int16_t,uint16_t>,fun<float32_t,int16_t,uint32_t>,fun<float32_t,int16_t,uint64_t>},{fun<float32_t,int32_t,uint8_t>,fun<float32_t,int32_t,uint16_t>,fun<float32_t,int32_t,uint32_t>,fun<float32_t,int32_t,uint64_t>},{fun<float32_t,int64_t,uint8_t>,fun<float32_t,int64_t,uint16_t>,fun<float32_t,int64_t,uint32_t>,fun<float32_t,int64_t,uint64_t>},{fun<float32_t,float32_t,uint8_t>,fun<float32_t,float32_t,uint16_t>,fun<float32_t,float32_t,uint32_t>,fun<float32_t,float32_t,uint64_t>},{fun<float32_t,float64_t,uint8_t>,fun<float32_t,float64_t,uint16_t>,fun<float32_t,float64_t,uint32_t>,fun<float32_t,float64_t,uint64_t>},{fun<float32_t,Complex64,uint8_t>,fun<float32_t,Complex64,uint16_t>,fun<float32_t,Complex64,uint32_t>,fun<float32_t,Complex64,uint64_t>},{fun<float32_t,Complex128,uint8_t>,fun<float32_t,Complex128,uint16_t>,fun<float32_t,Complex128,uint32_t>,fun<float32_t,Complex128,uint64_t>},{fun<float32_t,Rational32,uint8_t>,fun<float32_t,Rational32,uint16_t>,fun<float32_t,Rational32,uint32_t>,fun<float32_t,Rational32,uint64_t>},{fun<float32_t,Rational64,uint8_t>,fun<float32_t,Rational64,uint16_t>,fun<float32_t,Rational64,uint32_t>,fun<float32_t,Rational64,uint64_t>},{fun<float32_t,Rational128,uint8_t>,fun<float32_t,Rational128,uint16_t>,fun<float32_t,Rational128,uint32_t>,fun<float32_t,Rational128,uint64_t>},{fun<float32_t,RubyObject,uint8_t>,fun<float32_t,RubyObject,uint16_t>,fun<float32_t,RubyObject,uint32_t>,fun<float32_t,RubyObject,uint64_t>}}, \
  {{fun<float32_t,uint8_t,uint8_t>,fun<float32_t,uint8_t,uint16_t>,fun<float32_t,uint8_t,uint32_t>,fun<float32_t,uint8_t,uint64_t>},{fun<float32_t,int8_t,uint8_t>,fun<float32_t,int8_t,uint16_t>,fun<float32_t,int8_t,uint32_t>,fun<float32_t,int8_t,uint64_t>},{fun<float32_t,int16_t,uint8_t>,fun<float32_t,int16_t,uint16_t>,fun<float32_t,int16_t,uint32_t>,fun<float32_t,int16_t,uint64_t>},{fun<float32_t,int32_t,uint8_t>,fun<float32_t,int32_t,uint16_t>,fun<float32_t,int32_t,uint32_t>,fun<float32_t,int32_t,uint64_t>},{fun<float32_t,int64_t,uint8_t>,fun<float32_t,int64_t,uint16_t>,fun<float32_t,int64_t,uint32_t>,fun<float32_t,int64_t,uint64_t>},{fun<float32_t,float32_t,uint8_t>,fun<float32_t,float32_t,uint16_t>,fun<float32_t,float32_t,uint32_t>,fun<float32_t,float32_t,uint64_t>},{fun<float32_t,float64_t,uint8_t>,fun<float32_t,float64_t,uint16_t>,fun<float32_t,float64_t,uint32_t>,fun<float32_t,float64_t,uint64_t>},{fun<float32_t,Complex64,uint8_t>,fun<float32_t,Complex64,uint16_t>,fun<float32_t,Complex64,uint32_t>,fun<float32_t,Complex64,uint64_t>},{fun<float32_t,Complex128,uint8_t>,fun<float32_t,Complex128,uint16_t>,fun<float32_t,Complex128,uint32_t>,fun<float32_t,Complex128,uint64_t>},{fun<float32_t,Rational32,uint8_t>,fun<float32_t,Rational32,uint16_t>,fun<float32_t,Rational32,uint32_t>,fun<float32_t,Rational32,uint64_t>},{fun<float32_t,Rational64,uint8_t>,fun<float32_t,Rational64,uint16_t>,fun<float32_t,Rational64,uint32_t>,fun<float32_t,Rational64,uint64_t>},{fun<float32_t,Rational128,uint8_t>,fun<float32_t,Rational128,uint16_t>,fun<float32_t,Rational128,uint32_t>,fun<float32_t,Rational128,uint64_t>},{fun<float32_t,RubyObject,uint8_t>,fun<float32_t,RubyObject,uint16_t>,fun<float32_t,RubyObject,uint32_t>,fun<float32_t,RubyObject,uint64_t>}}, \
  {{fun<float64_t,uint8_t,uint8_t>,fun<float64_t,uint8_t,uint16_t>,fun<float64_t,uint8_t,uint32_t>,fun<float64_t,uint8_t,uint64_t>},{fun<float64_t,int8_t,uint8_t>,fun<float64_t,int8_t,uint16_t>,fun<float64_t,int8_t,uint32_t>,fun<float64_t,int8_t,uint64_t>},{fun<float64_t,int16_t,uint8_t>,fun<float64_t,int16_t,uint16_t>,fun<float64_t,int16_t,uint32_t>,fun<float64_t,int16_t,uint64_t>},{fun<float64_t,int32_t,uint8_t>,fun<float64_t,int32_t,uint16_t>,fun<float64_t,int32_t,uint32_t>,fun<float64_t,int32_t,uint64_t>},{fun<float64_t,int64_t,uint8_t>,fun<float64_t,int64_t,uint16_t>,fun<float64_t,int64_t,uint32_t>,fun<float64_t,int64_t,uint64_t>},{fun<float64_t,float32_t,uint8_t>,fun<float64_t,float32_t,uint16_t>,fun<float64_t,float32_t,uint32_t>,fun<float64_t,float32_t,uint64_t>},{fun<float64_t,float64_t,uint8_t>,fun<float64_t,float64_t,uint16_t>,fun<float64_t,float64_t,uint32_t>,fun<float64_t,float64_t,uint64_t>},{fun<float64_t,Complex64,uint8_t>,fun<float64_t,Complex64,uint16_t>,fun<float64_t,Complex64,uint32_t>,fun<float64_t,Complex64,uint64_t>},{fun<float64_t,Complex128,uint8_t>,fun<float64_t,Complex128,uint16_t>,fun<float64_t,Complex128,uint32_t>,fun<float64_t,Complex128,uint64_t>},{fun<float64_t,Rational32,uint8_t>,fun<float64_t,Rational32,uint16_t>,fun<float64_t,Rational32,uint32_t>,fun<float64_t,Rational32,uint64_t>},{fun<float64_t,Rational64,uint8_t>,fun<float64_t,Rational64,uint16_t>,fun<float64_t,Rational64,uint32_t>,fun<float64_t,Rational64,uint64_t>},{fun<float64_t,Rational128,uint8_t>,fun<float64_t,Rational128,uint16_t>,fun<float64_t,Rational128,uint32_t>,fun<float64_t,Rational128,uint64_t>},{fun<float64_t,RubyObject,uint8_t>,fun<float64_t,RubyObject,uint16_t>,fun<float64_t,RubyObject,uint32_t>,fun<float64_t,RubyObject,uint64_t>}}, \
  {{fun<float64_t,uint8_t,uint8_t>,fun<float64_t,uint8_t,uint16_t>,fun<float64_t,uint8_t,uint32_t>,fun<float64_t,uint8_t,uint64_t>},{fun<float64_t,int8_t,uint8_t>,fun<float64_t,int8_t,uint16_t>,fun<float64_t,int8_t,uint32_t>,fun<float64_t,int8_t,uint64_t>},{fun<float64_t,int16_t,uint8_t>,fun<float64_t,int16_t,uint16_t>,fun<float64_t,int16_t,uint32_t>,fun<float64_t,int16_t,uint64_t>},{fun<float64_t,int32_t,uint8_t>,fun<float64_t,int32_t,uint16_t>,fun<float64_t,int32_t,uint32_t>,fun<float64_t,int32_t,uint64_t>},{fun<float64_t,int64_t,uint8_t>,fun<float64_t,int64_t,uint16_t>,fun<float64_t,int64_t,uint32_t>,fun<float64_t,int64_t,uint64_t>},{fun<float64_t,float32_t,uint8_t>,fun<float64_t,float32_t,uint16_t>,fun<float64_t,float32_t,uint32_t>,fun<float64_t,float32_t,uint64_t>},{fun<float64_t,float64_t,uint8_t>,fun<float64_t,float64_t,uint16_t>,fun<float64_t,float64_t,uint32_t>,fun<float64_t,float64_t,uint64_t>},{fun<float64_t,Complex64,uint8_t>,fun<float64_t,Complex64,uint16_t>,fun<float64_t,Complex64,uint32_t>,fun<float64_t,Complex64,uint64_t>},{fun<float64_t,Complex128,uint8_t>,fun<float64_t,Complex128,uint16_t>,fun<float64_t,Complex128,uint32_t>,fun<float64_t,Complex128,uint64_t>},{fun<float64_t,Rational32,uint8_t>,fun<float64_t,Rational32,uint16_t>,fun<float64_t,Rational32,uint32_t>,fun<float64_t,Rational32,uint64_t>},{fun<float64_t,Rational64,uint8_t>,fun<float64_t,Rational64,uint16_t>,fun<float64_t,Rational64,uint32_t>,fun<float64_t,Rational64,uint64_t>},{fun<float64_t,Rational128,uint8_t>,fun<float64_t,Rational128,uint16_t>,fun<float64_t,Rational128,uint32_t>,fun<float64_t,Rational128,uint64_t>},{fun<float64_t,RubyObject,uint8_t>,fun<float64_t,RubyObject,uint16_t>,fun<float64_t,RubyObject,uint32_t>,fun<float64_t,RubyObject,uint64_t>}}, \
  {{fun<Complex64,uint8_t,uint8_t>,fun<Complex64,uint8_t,uint16_t>,fun<Complex64,uint8_t,uint32_t>,fun<Complex64,uint8_t,uint64_t>},{fun<Complex64,int8_t,uint8_t>,fun<Complex64,int8_t,uint16_t>,fun<Complex64,int8_t,uint32_t>,fun<Complex64,int8_t,uint64_t>},{fun<Complex64,int16_t,uint8_t>,fun<Complex64,int16_t,uint16_t>,fun<Complex64,int16_t,uint32_t>,fun<Complex64,int16_t,uint64_t>},{fun<Complex64,int32_t,uint8_t>,fun<Complex64,int32_t,uint16_t>,fun<Complex64,int32_t,uint32_t>,fun<Complex64,int32_t,uint64_t>},{fun<Complex64,int64_t,uint8_t>,fun<Complex64,int64_t,uint16_t>,fun<Complex64,int64_t,uint32_t>,fun<Complex64,int64_t,uint64_t>},{fun<Complex64,float32_t,uint8_t>,fun<Complex64,float32_t,uint16_t>,fun<Complex64,float32_t,uint32_t>,fun<Complex64,float32_t,uint64_t>},{fun<Complex64,float64_t,uint8_t>,fun<Complex64,float64_t,uint16_t>,fun<Complex64,float64_t,uint32_t>,fun<Complex64,float64_t,uint64_t>},{fun<Complex64,Complex64,uint8_t>,fun<Complex64,Complex64,uint16_t>,fun<Complex64,Complex64,uint32_t>,fun<Complex64,Complex64,uint64_t>},{fun<Complex64,Complex128,uint8_t>,fun<Complex64,Complex128,uint16_t>,fun<Complex64,Complex128,uint32_t>,fun<Complex64,Complex128,uint64_t>},{fun<Complex64,Rational32,uint8_t>,fun<Complex64,Rational32,uint16_t>,fun<Complex64,Rational32,uint32_t>,fun<Complex64,Rational32,uint64_t>},{fun<Complex64,Rational64,uint8_t>,fun<Complex64,Rational64,uint16_t>,fun<Complex64,Rational64,uint32_t>,fun<Complex64,Rational64,uint64_t>},{fun<Complex64,Rational128,uint8_t>,fun<Complex64,Rational128,uint16_t>,fun<Complex64,Rational128,uint32_t>,fun<Complex64,Rational128,uint64_t>},{fun<Complex64,RubyObject,uint8_t>,fun<Complex64,RubyObject,uint16_t>,fun<Complex64,RubyObject,uint32_t>,fun<Complex64,RubyObject,uint64_t>}}, \
  {{fun<Complex64,uint8_t,uint8_t>,fun<Complex64,uint8_t,uint16_t>,fun<Complex64,uint8_t,uint32_t>,fun<Complex64,uint8_t,uint64_t>},{fun<Complex64,int8_t,uint8_t>,fun<Complex64,int8_t,uint16_t>,fun<Complex64,int8_t,uint32_t>,fun<Complex64,int8_t,uint64_t>},{fun<Complex64,int16_t,uint8_t>,fun<Complex64,int16_t,uint16_t>,fun<Complex64,int16_t,uint32_t>,fun<Complex64,int16_t,uint64_t>},{fun<Complex64,int32_t,uint8_t>,fun<Complex64,int32_t,uint16_t>,fun<Complex64,int32_t,uint32_t>,fun<Complex64,int32_t,uint64_t>},{fun<Complex64,int64_t,uint8_t>,fun<Complex64,int64_t,uint16_t>,fun<Complex64,int64_t,uint32_t>,fun<Complex64,int64_t,uint64_t>},{fun<Complex64,float32_t,uint8_t>,fun<Complex64,float32_t,uint16_t>,fun<Complex64,float32_t,uint32_t>,fun<Complex64,float32_t,uint64_t>},{fun<Complex64,float64_t,uint8_t>,fun<Complex64,float64_t,uint16_t>,fun<Complex64,float64_t,uint32_t>,fun<Complex64,float64_t,uint64_t>},{fun<Complex64,Complex64,uint8_t>,fun<Complex64,Complex64,uint16_t>,fun<Complex64,Complex64,uint32_t>,fun<Complex64,Complex64,uint64_t>},{fun<Complex64,Complex128,uint8_t>,fun<Complex64,Complex128,uint16_t>,fun<Complex64,Complex128,uint32_t>,fun<Complex64,Complex128,uint64_t>},{fun<Complex64,Rational32,uint8_t>,fun<Complex64,Rational32,uint16_t>,fun<Complex64,Rational32,uint32_t>,fun<Complex64,Rational32,uint64_t>},{fun<Complex64,Rational64,uint8_t>,fun<Complex64,Rational64,uint16_t>,fun<Complex64,Rational64,uint32_t>,fun<Complex64,Rational64,uint64_t>},{fun<Complex64,Rational128,uint8_t>,fun<Complex64,Rational128,uint16_t>,fun<Complex64,Rational128,uint32_t>,fun<Complex64,Rational128,uint64_t>},{fun<Complex64,RubyObject,uint8_t>,fun<Complex64,RubyObject,uint16_t>,fun<Complex64,RubyObject,uint32_t>,fun<Complex64,RubyObject,uint64_t>}}, \
  {{fun<Complex128,uint8_t,uint8_t>,fun<Complex128,uint8_t,uint16_t>,fun<Complex128,uint8_t,uint32_t>,fun<Complex128,uint8_t,uint64_t>},{fun<Complex128,int8_t,uint8_t>,fun<Complex128,int8_t,uint16_t>,fun<Complex128,int8_t,uint32_t>,fun<Complex128,int8_t,uint64_t>},{fun<Complex128,int16_t,uint8_t>,fun<Complex128,int16_t,uint16_t>,fun<Complex128,int16_t,uint32_t>,fun<Complex128,int16_t,uint64_t>},{fun<Complex128,int32_t,uint8_t>,fun<Complex128,int32_t,uint16_t>,fun<Complex128,int32_t,uint32_t>,fun<Complex128,int32_t,uint64_t>},{fun<Complex128,int64_t,uint8_t>,fun<Complex128,int64_t,uint16_t>,fun<Complex128,int64_t,uint32_t>,fun<Complex128,int64_t,uint64_t>},{fun<Complex128,float32_t,uint8_t>,fun<Complex128,float32_t,uint16_t>,fun<Complex128,float32_t,uint32_t>,fun<Complex128,float32_t,uint64_t>},{fun<Complex128,float64_t,uint8_t>,fun<Complex128,float64_t,uint16_t>,fun<Complex128,float64_t,uint32_t>,fun<Complex128,float64_t,uint64_t>},{fun<Complex128,Complex64,uint8_t>,fun<Complex128,Complex64,uint16_t>,fun<Complex128,Complex64,uint32_t>,fun<Complex128,Complex64,uint64_t>},{fun<Complex128,Complex128,uint8_t>,fun<Complex128,Complex128,uint16_t>,fun<Complex128,Complex128,uint32_t>,fun<Complex128,Complex128,uint64_t>},{fun<Complex128,Rational32,uint8_t>,fun<Complex128,Rational32,uint16_t>,fun<Complex128,Rational32,uint32_t>,fun<Complex128,Rational32,uint64_t>},{fun<Complex128,Rational64,uint8_t>,fun<Complex128,Rational64,uint16_t>,fun<Complex128,Rational64,uint32_t>,fun<Complex128,Rational64,uint64_t>},{fun<Complex128,Rational128,uint8_t>,fun<Complex128,Rational128,uint16_t>,fun<Complex128,Rational128,uint32_t>,fun<Complex128,Rational128,uint64_t>},{fun<Complex128,RubyObject,uint8_t>,fun<Complex128,RubyObject,uint16_t>,fun<Complex128,RubyObject,uint32_t>,fun<Complex128,RubyObject,uint64_t>}}, \
  {{fun<Complex128,uint8_t,uint8_t>,fun<Complex128,uint8_t,uint16_t>,fun<Complex128,uint8_t,uint32_t>,fun<Complex128,uint8_t,uint64_t>},{fun<Complex128,int8_t,uint8_t>,fun<Complex128,int8_t,uint16_t>,fun<Complex128,int8_t,uint32_t>,fun<Complex128,int8_t,uint64_t>},{fun<Complex128,int16_t,uint8_t>,fun<Complex128,int16_t,uint16_t>,fun<Complex128,int16_t,uint32_t>,fun<Complex128,int16_t,uint64_t>},{fun<Complex128,int32_t,uint8_t>,fun<Complex128,int32_t,uint16_t>,fun<Complex128,int32_t,uint32_t>,fun<Complex128,int32_t,uint64_t>},{fun<Complex128,int64_t,uint8_t>,fun<Complex128,int64_t,uint16_t>,fun<Complex128,int64_t,uint32_t>,fun<Complex128,int64_t,uint64_t>},{fun<Complex128,float32_t,uint8_t>,fun<Complex128,float32_t,uint16_t>,fun<Complex128,float32_t,uint32_t>,fun<Complex128,float32_t,uint64_t>},{fun<Complex128,float64_t,uint8_t>,fun<Complex128,float64_t,uint16_t>,fun<Complex128,float64_t,uint32_t>,fun<Complex128,float64_t,uint64_t>},{fun<Complex128,Complex64,uint8_t>,fun<Complex128,Complex64,uint16_t>,fun<Complex128,Complex64,uint32_t>,fun<Complex128,Complex64,uint64_t>},{fun<Complex128,Complex128,uint8_t>,fun<Complex128,Complex128,uint16_t>,fun<Complex128,Complex128,uint32_t>,fun<Complex128,Complex128,uint64_t>},{fun<Complex128,Rational32,uint8_t>,fun<Complex128,Rational32,uint16_t>,fun<Complex128,Rational32,uint32_t>,fun<Complex128,Rational32,uint64_t>},{fun<Complex128,Rational64,uint8_t>,fun<Complex128,Rational64,uint16_t>,fun<Complex128,Rational64,uint32_t>,fun<Complex128,Rational64,uint64_t>},{fun<Complex128,Rational128,uint8_t>,fun<Complex128,Rational128,uint16_t>,fun<Complex128,Rational128,uint32_t>,fun<Complex128,Rational128,uint64_t>},{fun<Complex128,RubyObject,uint8_t>,fun<Complex128,RubyObject,uint16_t>,fun<Complex128,RubyObject,uint32_t>,fun<Complex128,RubyObject,uint64_t>}}, \
  {{fun<Rational32,uint8_t,uint8_t>,fun<Rational32,uint8_t,uint16_t>,fun<Rational32,uint8_t,uint32_t>,fun<Rational32,uint8_t,uint64_t>},{fun<Rational32,int8_t,uint8_t>,fun<Rational32,int8_t,uint16_t>,fun<Rational32,int8_t,uint32_t>,fun<Rational32,int8_t,uint64_t>},{fun<Rational32,int16_t,uint8_t>,fun<Rational32,int16_t,uint16_t>,fun<Rational32,int16_t,uint32_t>,fun<Rational32,int16_t,uint64_t>},{fun<Rational32,int32_t,uint8_t>,fun<Rational32,int32_t,uint16_t>,fun<Rational32,int32_t,uint32_t>,fun<Rational32,int32_t,uint64_t>},{fun<Rational32,int64_t,uint8_t>,fun<Rational32,int64_t,uint16_t>,fun<Rational32,int64_t,uint32_t>,fun<Rational32,int64_t,uint64_t>},{fun<Rational32,float32_t,uint8_t>,fun<Rational32,float32_t,uint16_t>,fun<Rational32,float32_t,uint32_t>,fun<Rational32,float32_t,uint64_t>},{fun<Rational32,float64_t,uint8_t>,fun<Rational32,float64_t,uint16_t>,fun<Rational32,float64_t,uint32_t>,fun<Rational32,float64_t,uint64_t>},{fun<Rational32,Complex64,uint8_t>,fun<Rational32,Complex64,uint16_t>,fun<Rational32,Complex64,uint32_t>,fun<Rational32,Complex64,uint64_t>},{fun<Rational32,Complex128,uint8_t>,fun<Rational32,Complex128,uint16_t>,fun<Rational32,Complex128,uint32_t>,fun<Rational32,Complex128,uint64_t>},{fun<Rational32,Rational32,uint8_t>,fun<Rational32,Rational32,uint16_t>,fun<Rational32,Rational32,uint32_t>,fun<Rational32,Rational32,uint64_t>},{fun<Rational32,Rational64,uint8_t>,fun<Rational32,Rational64,uint16_t>,fun<Rational32,Rational64,uint32_t>,fun<Rational32,Rational64,uint64_t>},{fun<Rational32,Rational128,uint8_t>,fun<Rational32,Rational128,uint16_t>,fun<Rational32,Rational128,uint32_t>,fun<Rational32,Rational128,uint64_t>},{fun<Rational32,RubyObject,uint8_t>,fun<Rational32,RubyObject,uint16_t>,fun<Rational32,RubyObject,uint32_t>,fun<Rational32,RubyObject,uint64_t>}}, \
  {{fun<Rational32,uint8_t,uint8_t>,fun<Rational32,uint8_t,uint16_t>,fun<Rational32,uint8_t,uint32_t>,fun<Rational32,uint8_t,uint64_t>},{fun<Rational32,int8_t,uint8_t>,fun<Rational32,int8_t,uint16_t>,fun<Rational32,int8_t,uint32_t>,fun<Rational32,int8_t,uint64_t>},{fun<Rational32,int16_t,uint8_t>,fun<Rational32,int16_t,uint16_t>,fun<Rational32,int16_t,uint32_t>,fun<Rational32,int16_t,uint64_t>},{fun<Rational32,int32_t,uint8_t>,fun<Rational32,int32_t,uint16_t>,fun<Rational32,int32_t,uint32_t>,fun<Rational32,int32_t,uint64_t>},{fun<Rational32,int64_t,uint8_t>,fun<Rational32,int64_t,uint16_t>,fun<Rational32,int64_t,uint32_t>,fun<Rational32,int64_t,uint64_t>},{fun<Rational32,float32_t,uint8_t>,fun<Rational32,float32_t,uint16_t>,fun<Rational32,float32_t,uint32_t>,fun<Rational32,float32_t,uint64_t>},{fun<Rational32,float64_t,uint8_t>,fun<Rational32,float64_t,uint16_t>,fun<Rational32,float64_t,uint32_t>,fun<Rational32,float64_t,uint64_t>},{fun<Rational32,Complex64,uint8_t>,fun<Rational32,Complex64,uint16_t>,fun<Rational32,Complex64,uint32_t>,fun<Rational32,Complex64,uint64_t>},{fun<Rational32,Complex128,uint8_t>,fun<Rational32,Complex128,uint16_t>,fun<Rational32,Complex128,uint32_t>,fun<Rational32,Complex128,uint64_t>},{fun<Rational32,Rational32,uint8_t>,fun<Rational32,Rational32,uint16_t>,fun<Rational32,Rational32,uint32_t>,fun<Rational32,Rational32,uint64_t>},{fun<Rational32,Rational64,uint8_t>,fun<Rational32,Rational64,uint16_t>,fun<Rational32,Rational64,uint32_t>,fun<Rational32,Rational64,uint64_t>},{fun<Rational32,Rational128,uint8_t>,fun<Rational32,Rational128,uint16_t>,fun<Rational32,Rational128,uint32_t>,fun<Rational32,Rational128,uint64_t>},{fun<Rational32,RubyObject,uint8_t>,fun<Rational32,RubyObject,uint16_t>,fun<Rational32,RubyObject,uint32_t>,fun<Rational32,RubyObject,uint64_t>}}, \
  {{fun<Rational64,uint8_t,uint8_t>,fun<Rational64,uint8_t,uint16_t>,fun<Rational64,uint8_t,uint32_t>,fun<Rational64,uint8_t,uint64_t>},{fun<Rational64,int8_t,uint8_t>,fun<Rational64,int8_t,uint16_t>,fun<Rational64,int8_t,uint32_t>,fun<Rational64,int8_t,uint64_t>},{fun<Rational64,int16_t,uint8_t>,fun<Rational64,int16_t,uint16_t>,fun<Rational64,int16_t,uint32_t>,fun<Rational64,int16_t,uint64_t>},{fun<Rational64,int32_t,uint8_t>,fun<Rational64,int32_t,uint16_t>,fun<Rational64,int32_t,uint32_t>,fun<Rational64,int32_t,uint64_t>},{fun<Rational64,int64_t,uint8_t>,fun<Rational64,int64_t,uint16_t>,fun<Rational64,int64_t,uint32_t>,fun<Rational64,int64_t,uint64_t>},{fun<Rational64,float32_t,uint8_t>,fun<Rational64,float32_t,uint16_t>,fun<Rational64,float32_t,uint32_t>,fun<Rational64,float32_t,uint64_t>},{fun<Rational64,float64_t,uint8_t>,fun<Rational64,float64_t,uint16_t>,fun<Rational64,float64_t,uint32_t>,fun<Rational64,float64_t,uint64_t>},{fun<Rational64,Complex64,uint8_t>,fun<Rational64,Complex64,uint16_t>,fun<Rational64,Complex64,uint32_t>,fun<Rational64,Complex64,uint64_t>},{fun<Rational64,Complex128,uint8_t>,fun<Rational64,Complex128,uint16_t>,fun<Rational64,Complex128,uint32_t>,fun<Rational64,Complex128,uint64_t>},{fun<Rational64,Rational32,uint8_t>,fun<Rational64,Rational32,uint16_t>,fun<Rational64,Rational32,uint32_t>,fun<Rational64,Rational32,uint64_t>},{fun<Rational64,Rational64,uint8_t>,fun<Rational64,Rational64,uint16_t>,fun<Rational64,Rational64,uint32_t>,fun<Rational64,Rational64,uint64_t>},{fun<Rational64,Rational128,uint8_t>,fun<Rational64,Rational128,uint16_t>,fun<Rational64,Rational128,uint32_t>,fun<Rational64,Rational128,uint64_t>},{fun<Rational64,RubyObject,uint8_t>,fun<Rational64,RubyObject,uint16_t>,fun<Rational64,RubyObject,uint32_t>,fun<Rational64,RubyObject,uint64_t>}}, \
  {{fun<Rational64,uint8_t,uint8_t>,fun<Rational64,uint8_t,uint16_t>,fun<Rational64,uint8_t,uint32_t>,fun<Rational64,uint8_t,uint64_t>},{fun<Rational64,int8_t,uint8_t>,fun<Rational64,int8_t,uint16_t>,fun<Rational64,int8_t,uint32_t>,fun<Rational64,int8_t,uint64_t>},{fun<Rational64,int16_t,uint8_t>,fun<Rational64,int16_t,uint16_t>,fun<Rational64,int16_t,uint32_t>,fun<Rational64,int16_t,uint64_t>},{fun<Rational64,int32_t,uint8_t>,fun<Rational64,int32_t,uint16_t>,fun<Rational64,int32_t,uint32_t>,fun<Rational64,int32_t,uint64_t>},{fun<Rational64,int64_t,uint8_t>,fun<Rational64,int64_t,uint16_t>,fun<Rational64,int64_t,uint32_t>,fun<Rational64,int64_t,uint64_t>},{fun<Rational64,float32_t,uint8_t>,fun<Rational64,float32_t,uint16_t>,fun<Rational64,float32_t,uint32_t>,fun<Rational64,float32_t,uint64_t>},{fun<Rational64,float64_t,uint8_t>,fun<Rational64,float64_t,uint16_t>,fun<Rational64,float64_t,uint32_t>,fun<Rational64,float64_t,uint64_t>},{fun<Rational64,Complex64,uint8_t>,fun<Rational64,Complex64,uint16_t>,fun<Rational64,Complex64,uint32_t>,fun<Rational64,Complex64,uint64_t>},{fun<Rational64,Complex128,uint8_t>,fun<Rational64,Complex128,uint16_t>,fun<Rational64,Complex128,uint32_t>,fun<Rational64,Complex128,uint64_t>},{fun<Rational64,Rational32,uint8_t>,fun<Rational64,Rational32,uint16_t>,fun<Rational64,Rational32,uint32_t>,fun<Rational64,Rational32,uint64_t>},{fun<Rational64,Rational64,uint8_t>,fun<Rational64,Rational64,uint16_t>,fun<Rational64,Rational64,uint32_t>,fun<Rational64,Rational64,uint64_t>},{fun<Rational64,Rational128,uint8_t>,fun<Rational64,Rational128,uint16_t>,fun<Rational64,Rational128,uint32_t>,fun<Rational64,Rational128,uint64_t>},{fun<Rational64,RubyObject,uint8_t>,fun<Rational64,RubyObject,uint16_t>,fun<Rational64,RubyObject,uint32_t>,fun<Rational64,RubyObject,uint64_t>}}, \
  {{fun<Rational128,uint8_t,uint8_t>,fun<Rational128,uint8_t,uint16_t>,fun<Rational128,uint8_t,uint32_t>,fun<Rational128,uint8_t,uint64_t>},{fun<Rational128,int8_t,uint8_t>,fun<Rational128,int8_t,uint16_t>,fun<Rational128,int8_t,uint32_t>,fun<Rational128,int8_t,uint64_t>},{fun<Rational128,int16_t,uint8_t>,fun<Rational128,int16_t,uint16_t>,fun<Rational128,int16_t,uint32_t>,fun<Rational128,int16_t,uint64_t>},{fun<Rational128,int32_t,uint8_t>,fun<Rational128,int32_t,uint16_t>,fun<Rational128,int32_t,uint32_t>,fun<Rational128,int32_t,uint64_t>},{fun<Rational128,int64_t,uint8_t>,fun<Rational128,int64_t,uint16_t>,fun<Rational128,int64_t,uint32_t>,fun<Rational128,int64_t,uint64_t>},{fun<Rational128,float32_t,uint8_t>,fun<Rational128,float32_t,uint16_t>,fun<Rational128,float32_t,uint32_t>,fun<Rational128,float32_t,uint64_t>},{fun<Rational128,float64_t,uint8_t>,fun<Rational128,float64_t,uint16_t>,fun<Rational128,float64_t,uint32_t>,fun<Rational128,float64_t,uint64_t>},{fun<Rational128,Complex64,uint8_t>,fun<Rational128,Complex64,uint16_t>,fun<Rational128,Complex64,uint32_t>,fun<Rational128,Complex64,uint64_t>},{fun<Rational128,Complex128,uint8_t>,fun<Rational128,Complex128,uint16_t>,fun<Rational128,Complex128,uint32_t>,fun<Rational128,Complex128,uint64_t>},{fun<Rational128,Rational32,uint8_t>,fun<Rational128,Rational32,uint16_t>,fun<Rational128,Rational32,uint32_t>,fun<Rational128,Rational32,uint64_t>},{fun<Rational128,Rational64,uint8_t>,fun<Rational128,Rational64,uint16_t>,fun<Rational128,Rational64,uint32_t>,fun<Rational128,Rational64,uint64_t>},{fun<Rational128,Rational128,uint8_t>,fun<Rational128,Rational128,uint16_t>,fun<Rational128,Rational128,uint32_t>,fun<Rational128,Rational128,uint64_t>},{fun<Rational128,RubyObject,uint8_t>,fun<Rational128,RubyObject,uint16_t>,fun<Rational128,RubyObject,uint32_t>,fun<Rational128,RubyObject,uint64_t>}}, \
  {{fun<Rational128,uint8_t,uint8_t>,fun<Rational128,uint8_t,uint16_t>,fun<Rational128,uint8_t,uint32_t>,fun<Rational128,uint8_t,uint64_t>},{fun<Rational128,int8_t,uint8_t>,fun<Rational128,int8_t,uint16_t>,fun<Rational128,int8_t,uint32_t>,fun<Rational128,int8_t,uint64_t>},{fun<Rational128,int16_t,uint8_t>,fun<Rational128,int16_t,uint16_t>,fun<Rational128,int16_t,uint32_t>,fun<Rational128,int16_t,uint64_t>},{fun<Rational128,int32_t,uint8_t>,fun<Rational128,int32_t,uint16_t>,fun<Rational128,int32_t,uint32_t>,fun<Rational128,int32_t,uint64_t>},{fun<Rational128,int64_t,uint8_t>,fun<Rational128,int64_t,uint16_t>,fun<Rational128,int64_t,uint32_t>,fun<Rational128,int64_t,uint64_t>},{fun<Rational128,float32_t,uint8_t>,fun<Rational128,float32_t,uint16_t>,fun<Rational128,float32_t,uint32_t>,fun<Rational128,float32_t,uint64_t>},{fun<Rational128,float64_t,uint8_t>,fun<Rational128,float64_t,uint16_t>,fun<Rational128,float64_t,uint32_t>,fun<Rational128,float64_t,uint64_t>},{fun<Rational128,Complex64,uint8_t>,fun<Rational128,Complex64,uint16_t>,fun<Rational128,Complex64,uint32_t>,fun<Rational128,Complex64,uint64_t>},{fun<Rational128,Complex128,uint8_t>,fun<Rational128,Complex128,uint16_t>,fun<Rational128,Complex128,uint32_t>,fun<Rational128,Complex128,uint64_t>},{fun<Rational128,Rational32,uint8_t>,fun<Rational128,Rational32,uint16_t>,fun<Rational128,Rational32,uint32_t>,fun<Rational128,Rational32,uint64_t>},{fun<Rational128,Rational64,uint8_t>,fun<Rational128,Rational64,uint16_t>,fun<Rational128,Rational64,uint32_t>,fun<Rational128,Rational64,uint64_t>},{fun<Rational128,Rational128,uint8_t>,fun<Rational128,Rational128,uint16_t>,fun<Rational128,Rational128,uint32_t>,fun<Rational128,Rational128,uint64_t>},{fun<Rational128,RubyObject,uint8_t>,fun<Rational128,RubyObject,uint16_t>,fun<Rational128,RubyObject,uint32_t>,fun<Rational128,RubyObject,uint64_t>}}, \
  {{fun<RubyObject,uint8_t,uint8_t>,fun<RubyObject,uint8_t,uint16_t>,fun<RubyObject,uint8_t,uint32_t>,fun<RubyObject,uint8_t,uint64_t>},{fun<RubyObject,int8_t,uint8_t>,fun<RubyObject,int8_t,uint16_t>,fun<RubyObject,int8_t,uint32_t>,fun<RubyObject,int8_t,uint64_t>},{fun<RubyObject,int16_t,uint8_t>,fun<RubyObject,int16_t,uint16_t>,fun<RubyObject,int16_t,uint32_t>,fun<RubyObject,int16_t,uint64_t>},{fun<RubyObject,int32_t,uint8_t>,fun<RubyObject,int32_t,uint16_t>,fun<RubyObject,int32_t,uint32_t>,fun<RubyObject,int32_t,uint64_t>},{fun<RubyObject,int64_t,uint8_t>,fun<RubyObject,int64_t,uint16_t>,fun<RubyObject,int64_t,uint32_t>,fun<RubyObject,int64_t,uint64_t>},{fun<RubyObject,float32_t,uint8_t>,fun<RubyObject,float32_t,uint16_t>,fun<RubyObject,float32_t,uint32_t>,fun<RubyObject,float32_t,uint64_t>},{fun<RubyObject,float64_t,uint8_t>,fun<RubyObject,float64_t,uint16_t>,fun<RubyObject,float64_t,uint32_t>,fun<RubyObject,float64_t,uint64_t>},{fun<RubyObject,Complex64,uint8_t>,fun<RubyObject,Complex64,uint16_t>,fun<RubyObject,Complex64,uint32_t>,fun<RubyObject,Complex64,uint64_t>},{fun<RubyObject,Complex128,uint8_t>,fun<RubyObject,Complex128,uint16_t>,fun<RubyObject,Complex128,uint32_t>,fun<RubyObject,Complex128,uint64_t>},{fun<RubyObject,Rational32,uint8_t>,fun<RubyObject,Rational32,uint16_t>,fun<RubyObject,Rational32,uint32_t>,fun<RubyObject,Rational32,uint64_t>},{fun<RubyObject,Rational64,uint8_t>,fun<RubyObject,Rational64,uint16_t>,fun<RubyObject,Rational64,uint32_t>,fun<RubyObject,Rational64,uint64_t>},{fun<RubyObject,Rational128,uint8_t>,fun<RubyObject,Rational128,uint16_t>,fun<RubyObject,Rational128,uint32_t>,fun<RubyObject,Rational128,uint64_t>},{fun<RubyObject,RubyObject,uint8_t>,fun<RubyObject,RubyObject,uint16_t>,fun<RubyObject,RubyObject,uint32_t>,fun<RubyObject,RubyObject,uint64_t>}} \
};

#define LRI_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...)    NAMED_LRI_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)


#define NAMED_LI_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...)\
static ret (*(name)[NUM_DTYPES][NUM_ITYPES])(__VA_ARGS__) = {\
  { fun<uint8_t,uint8_t>,fun<uint8_t,uint16_t>,fun<uint8_t,uint32_t>,fun<uint8_t,uint64_t>},\
  { fun<int8_t,uint8_t>,fun<int8_t,uint16_t>,fun<int8_t,uint32_t>,fun<int8_t,uint64_t>},\
  { fun<int16_t,uint8_t>,fun<int16_t,uint16_t>,fun<int16_t,uint32_t>,fun<int16_t,uint64_t>},\
  { fun<int32_t,uint8_t>,fun<int32_t,uint16_t>,fun<int32_t,uint32_t>,fun<int32_t,uint64_t>},\
  { fun<int64_t,uint8_t>,fun<int64_t,uint16_t>,fun<int64_t,uint32_t>,fun<int64_t,uint64_t>},\
  { fun<float32_t,uint8_t>,fun<float32_t,uint16_t>,fun<float32_t,uint32_t>,fun<float32_t,uint64_t>},\
  { fun<float64_t,uint8_t>,fun<float64_t,uint16_t>,fun<float64_t,uint32_t>,fun<float64_t,uint64_t>},\
  { fun<Complex64,uint8_t>,fun<Complex64,uint16_t>,fun<Complex64,uint32_t>,fun<Complex64,uint64_t>},\
  { fun<Complex128,uint8_t>,fun<Complex128,uint16_t>,fun<Complex128,uint32_t>,fun<Complex128,uint64_t>},\
  { fun<Rational32,uint8_t>,fun<Rational32,uint16_t>,fun<Rational32,uint32_t>,fun<Rational32,uint64_t>},\
  { fun<Rational64,uint8_t>,fun<Rational64,uint16_t>,fun<Rational64,uint32_t>,fun<Rational64,uint64_t>},\
  { fun<Rational128,uint8_t>,fun<Rational128,uint16_t>,fun<Rational128,uint32_t>,fun<Rational128,uint64_t>},\
  { fun<RubyObject,uint8_t>,fun<RubyObject,uint16_t>,fun<RubyObject,uint32_t>,fun<RubyObject,uint64_t>} \
};

#define LI_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...)    NAMED_LI_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)



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

/*
 * Index Types
 */

typedef enum {
  UINT8 = 0,
  UINT16 = 1,
  UINT32 = 2,
  UINT64 = 3
} itype_t;

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

// regular data types
extern const char* const	DTYPE_NAMES[NUM_DTYPES];
extern const size_t 			DTYPE_SIZES[NUM_DTYPES];

// index data types
extern const char* const  ITYPE_NAMES[NUM_ITYPES];
extern const size_t 			ITYPE_SIZES[NUM_ITYPES];

extern const dtype_t      Upcast[NUM_DTYPES][NUM_DTYPES];

/*
 * Functions
 */

void*				rubyobj_to_cval(VALUE val, dtype_t dtype);
void				rubyval_to_cval(VALUE val, dtype_t dtype, void* loc);
RubyObject	rubyobj_from_cval(void* val, dtype_t dtype);

#endif // DATA_H
