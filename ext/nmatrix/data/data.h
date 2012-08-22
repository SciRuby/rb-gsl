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

#include "nmatrix.h"

#include "types.h"

#include "complex.h"
#include "rational.h"
#include "ruby_object.h"

namespace nm {
  /*
   * Constants
   */

  const size_t NUM_DTYPES = 13;
  const size_t NUM_ITYPES = 4;
  const size_t NUM_EWOPS = 5;

  enum ewop_t {
    EW_ADD,
    EW_SUB,
    EW_MUL,
    EW_DIV,
    EW_MOD
  };
  
} // end of namespace nm

/*
 * Macros
 */

extern "C" {

#define STYPE_MARK_TABLE(name)									\
  static void (*(name)[nm::NUM_STYPES])(void*) = {	\
    nm_dense_storage_mark,											\
    nm_list_storage_mark,												\
    nm_yale_storage_mark												\
  };

#define STYPE_CAST_COPY_TABLE(name)                                                   \
  static STORAGE* (*(name)[nm::NUM_STYPES][nm::NUM_STYPES])(const STORAGE*, dtype_t) = {      \
    { nm_dense_storage_cast_copy,  nm_dense_storage_from_list,  nm_dense_storage_from_yale },  \
    { nm_list_storage_from_dense,  nm_list_storage_cast_copy,   nm_list_storage_from_yale  },  \
    { nm_yale_storage_from_dense,  nm_yale_storage_from_list,   nm_yale_storage_cast_copy  }   \
  };

/*
 * Defines a static array that hold function pointers to dtype templated
 * versions of the specified function.
 */
#define DTYPE_TEMPLATE_TABLE(fun, ret, ...) NAMED_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)

#define NAMED_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...) \
	static ret (*(name)[nm::NUM_DTYPES])(__VA_ARGS__) =	{			\
		fun<uint8_t>,																				\
		fun<int8_t>,																				\
		fun<int16_t>,																				\
		fun<int32_t>,																				\
		fun<int64_t>,																				\
		fun<float32_t>,																			\
		fun<float64_t>,																			\
		fun<nm::Complex64>,																  \
		fun<nm::Complex128>,																\
		fun<nm::Rational32>,																\
		fun<nm::Rational64>,																\
		fun<nm::Rational128>,																\
		fun<nm::RubyObject>																	\
	};

/*
 * Defines a static array that holds function pointers to itype templated
 * versions of the specified function.
 */
#define ITYPE_TEMPLATE_TABLE(fun, ret, ...) NAMED_ITYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)

#define NAMED_ITYPE_TEMPLATE_TABLE(name, fun, ret, ...)	\
  static ret (*(name)[nm::NUM_ITYPES])(__VA_ARGS__) = {			\
    fun<uint8_t>,																				\
    fun<uint16_t>,																			\
    fun<uint32_t>,																			\
    fun<uint64_t>																				\
  };

#define ITYPE_TEMPLATE_TABLE(fun, ret, ...) NAMED_ITYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)



/*
 * Same as DTYPE_TEMPLATE_TABLE but for functions that have two template
 * parameters.
 *
 * The left-hand DType is used as the first index, and the right-hand side is
 * the second index.  Not all left- and right-hand side combinations are valid,
 * and an invalid combination will result in a NULL pointer.
 */
#define LR_DTYPE_TEMPLATE_TABLE(fun, ret, ...) NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)

#define NAMED_LR_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...)																																																								\
	static ret (*(name)[nm::NUM_DTYPES][nm::NUM_DTYPES])(__VA_ARGS__) = {																																																						\
		{fun<uint8_t, uint8_t>, fun<uint8_t, int8_t>, fun<uint8_t, int16_t>, fun<uint8_t, int32_t>, fun<uint8_t, int64_t>, fun<uint8_t, float32_t>, fun<uint8_t, float64_t>,	\
			fun<uint8_t, nm::Complex64>, fun<uint8_t, nm::Complex128>, fun<uint8_t, nm::Rational32>, fun<uint8_t, nm::Rational64>,																							\
			fun<uint8_t, nm::Rational128>, NULL},																																																																\
																																																																																					\
		{fun<int8_t, uint8_t>, fun<int8_t, int8_t>, fun<int8_t, int16_t>, fun<int8_t, int32_t>, fun<int8_t, int64_t>, fun<int8_t, float32_t>, fun<int8_t, float64_t>,					\
			fun<int8_t, nm::Complex64>, fun<int8_t, nm::Complex128>, fun<int8_t, nm::Rational32>, fun<int8_t, nm::Rational64>, fun<int8_t, nm::Rational128>, NULL},							\
																																																																																					\
		{fun<int16_t, uint8_t>, fun<int16_t, int8_t>, fun<int16_t, int16_t>, fun<int16_t, int32_t>, fun<int16_t, int64_t>, fun<int16_t, float32_t>, fun<int16_t, float64_t>,	\
			fun<int16_t, nm::Complex64>, fun<int16_t, nm::Complex128>, fun<int16_t, nm::Rational32>, fun<int16_t, nm::Rational64>, fun<int16_t, nm::Rational128>, NULL},				\
																																																																																					\
		{fun<int32_t, uint8_t>, fun<int32_t, int8_t>, fun<int32_t, int16_t>, fun<int32_t, int32_t>, fun<int32_t, int64_t>, fun<int32_t, float32_t>, fun<int32_t, float64_t>,	\
			fun<int32_t, nm::Complex64>, fun<int32_t, nm::Complex128>, fun<int32_t, nm::Rational32>, fun<int32_t, nm::Rational64>, fun<int32_t, nm::Rational128>, NULL},				\
																																																																																					\
		{fun<int64_t, uint8_t>, fun<int64_t, int8_t>, fun<int64_t, int16_t>, fun<int64_t, int32_t>, fun<int64_t, int64_t>, fun<int64_t, float32_t>, fun<int64_t, float64_t>,	\
			fun<int64_t, nm::Complex64>, fun<int64_t, nm::Complex128>, fun<int64_t, nm::Rational32>, fun<int64_t, nm::Rational64>, fun<int64_t, nm::Rational128>, NULL},				\
																																																																																					\
		{fun<float32_t, uint8_t>, fun<float32_t, int8_t>, fun<float32_t, int16_t>, fun<float32_t, int32_t>, fun<float32_t, int64_t>,																					\
			fun<float32_t, float32_t>, fun<float32_t, float64_t>, fun<float32_t, nm::Complex64>, fun<float32_t, nm::Complex128>,  fun<float32_t, nm::Rational32>,								\
			fun<float32_t, nm::Rational64>, fun<float32_t, nm::Rational128>, NULL},																																															\
                                                                                                                                                                          \
		{fun<float64_t, uint8_t>, fun<float64_t, int8_t>, fun<float64_t, int16_t>, fun<float64_t, int32_t>, fun<float64_t, int64_t>,																					\
			fun<float64_t, float32_t>, fun<float64_t, float64_t>, fun<float64_t, nm::Complex64>, fun<float64_t, nm::Complex128>, fun<float64_t, nm::Rational32>,                \
			fun<float64_t, nm::Rational64>, fun<float64_t, nm::Rational128>, NULL},																																															\
                                                                                                                                                                          \
		{fun<nm::Complex64, uint8_t>, fun<nm::Complex64, int8_t>, fun<nm::Complex64, int16_t>, fun<nm::Complex64, int32_t>, fun<nm::Complex64, int64_t>,											\
			fun<nm::Complex64, float32_t>, fun<nm::Complex64, float64_t>, fun<nm::Complex64, nm::Complex64>, fun<nm::Complex64, nm::Complex128>,																\
			fun<nm::Complex64, nm::Rational32>, fun<nm::Complex64, nm::Rational64>, fun<nm::Complex64, nm::Rational128>, NULL},																									\
																																																																																					\
		{fun<nm::Complex128, uint8_t>, fun<nm::Complex128, int8_t>, fun<nm::Complex128, int16_t>, fun<nm::Complex128, int32_t>, fun<nm::Complex128, int64_t>,									\
			fun<nm::Complex128, float32_t>, fun<nm::Complex128, float64_t>, fun<nm::Complex128, nm::Complex64>, fun<nm::Complex128, nm::Complex128>,														\
			fun<nm::Complex128, nm::Rational32>, fun<nm::Complex128, nm::Rational64>, fun<nm::Complex128, nm::Rational128>, NULL},																							\
																																																																																					\
		{fun<nm::Rational32, uint8_t>, fun<nm::Rational32, int8_t>, fun<nm::Rational32, int16_t>, fun<nm::Rational32, int32_t>, fun<nm::Rational32, int64_t>, NULL, NULL,			\
			NULL, NULL, fun<nm::Rational32, nm::Rational32>, fun<nm::Rational32, nm::Rational64>, fun<nm::Rational32, nm::Rational128>, NULL},																	\
																																																																																					\
		{fun<nm::Rational64, uint8_t>, fun<nm::Rational64, int8_t>, fun<nm::Rational64, int16_t>, fun<nm::Rational64, int32_t>, fun<nm::Rational64, int64_t>, NULL, NULL,			\
			NULL, NULL, fun<nm::Rational64, nm::Rational32>, fun<nm::Rational64, nm::Rational64>, fun<nm::Rational64, nm::Rational128>, NULL},																	\
																																																																																					\
		{fun<nm::Rational128, uint8_t>, fun<nm::Rational128, int8_t>, fun<nm::Rational128, int16_t>, fun<nm::Rational128, int32_t>, fun<nm::Rational128, int64_t>, NULL,			\
			NULL, NULL, NULL, fun<nm::Rational128, nm::Rational32>, fun<nm::Rational128, nm::Rational64>, fun<nm::Rational128, nm::Rational128>, NULL},													\
																																																																																					\
		{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::RubyObject, nm::RubyObject>}																													\
	};

/*
 * Defines a static array that holds function pointers to operation, and left-
 * and right-side dtype templated version sof the specified function.
 */
#define OP_LR_DTYPE_TEMPLATE_TABLE(fun, ret, ...) NAMED_OP_LR_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)

#define NAMED_OP_LR_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...) 																																																							\
	static ret (*(name)[nm::NUM_EWOPS][nm::NUM_DTYPES][nm::NUM_DTYPES])(__VA_ARGS__) = {																																																	\
		{																																																																																				\
			{fun<nm::EW_ADD, uint8_t, uint8_t>, fun<nm::EW_ADD, uint8_t, int8_t>, fun<nm::EW_ADD, uint8_t, int16_t>, fun<nm::EW_ADD, uint8_t, int32_t>, fun<nm::EW_ADD, uint8_t, int64_t>,						\
				fun<nm::EW_ADD, uint8_t, float32_t>, fun<nm::EW_ADD, uint8_t, float64_t>, fun<nm::EW_ADD, uint8_t, nm::Complex64>, fun<nm::EW_ADD, uint8_t, nm::Complex128>,												\
				fun<nm::EW_ADD, uint8_t, nm::Rational32>, fun<nm::EW_ADD, uint8_t, nm::Rational64>, fun<nm::EW_ADD, uint8_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_ADD, int8_t, uint8_t>, fun<nm::EW_ADD, int8_t, int8_t>, fun<nm::EW_ADD, int8_t, int16_t>, fun<nm::EW_ADD, int8_t, int32_t>, fun<nm::EW_ADD, int8_t, int64_t>,									\
				fun<nm::EW_ADD, int8_t, float32_t>, fun<nm::EW_ADD, int8_t, float64_t>, fun<nm::EW_ADD, int8_t, nm::Complex64>, fun<nm::EW_ADD, int8_t, nm::Complex128>,														\
				fun<nm::EW_ADD, int8_t, nm::Rational32>, fun<nm::EW_ADD, int8_t, nm::Rational64>, fun<nm::EW_ADD, int8_t, nm::Rational128>, NULL},																							\
																																																																																						\
			{fun<nm::EW_ADD, int16_t, uint8_t>, fun<nm::EW_ADD, int16_t, int8_t>, fun<nm::EW_ADD, int16_t, int16_t>, fun<nm::EW_ADD, int16_t, int32_t>, fun<nm::EW_ADD, int16_t, int64_t>,						\
				fun<nm::EW_ADD, int16_t, float32_t>, fun<nm::EW_ADD, int16_t, float64_t>, fun<nm::EW_ADD, int16_t, nm::Complex64>, fun<nm::EW_ADD, int16_t, nm::Complex128>,												\
				fun<nm::EW_ADD, int16_t, nm::Rational32>, fun<nm::EW_ADD, int16_t, nm::Rational64>, fun<nm::EW_ADD, int16_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_ADD, int32_t, uint8_t>, fun<nm::EW_ADD, int32_t, int8_t>, fun<nm::EW_ADD, int32_t, int16_t>, fun<nm::EW_ADD, int32_t, int32_t>, fun<nm::EW_ADD, int32_t, int64_t>,						\
				fun<nm::EW_ADD, int32_t, float32_t>, fun<nm::EW_ADD, int32_t, float64_t>, fun<nm::EW_ADD, int32_t, nm::Complex64>, fun<nm::EW_ADD, int32_t, nm::Complex128>,												\
				fun<nm::EW_ADD, int32_t, nm::Rational32>, fun<nm::EW_ADD, int32_t, nm::Rational64>, fun<nm::EW_ADD, int32_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_ADD, int64_t, uint8_t>, fun<nm::EW_ADD, int64_t, int8_t>, fun<nm::EW_ADD, int64_t, int16_t>, fun<nm::EW_ADD, int64_t, int32_t>, fun<nm::EW_ADD, int64_t, int64_t>,						\
				fun<nm::EW_ADD, int64_t, float32_t>, fun<nm::EW_ADD, int64_t, float64_t>, fun<nm::EW_ADD, int64_t, nm::Complex64>, fun<nm::EW_ADD, int64_t, nm::Complex128>,												\
				fun<nm::EW_ADD, int64_t, nm::Rational32>, fun<nm::EW_ADD, int64_t, nm::Rational64>, fun<nm::EW_ADD, int64_t, nm::Rational128>, NULL}, 																					\
																																																																																						\
			{fun<nm::EW_ADD, float32_t, uint8_t>, fun<nm::EW_ADD, float32_t, int8_t>, fun<nm::EW_ADD, float32_t, int16_t>, fun<nm::EW_ADD, float32_t, int32_t>, fun<nm::EW_ADD, float32_t, int64_t>,	\
				fun<nm::EW_ADD, float32_t, float32_t>, fun<nm::EW_ADD, float32_t, float64_t>, fun<nm::EW_ADD, float32_t, nm::Complex64>, fun<nm::EW_ADD, float32_t, nm::Complex128>,								\
				fun<nm::EW_ADD, float32_t, nm::Rational32>, fun<nm::EW_ADD, float32_t, nm::Rational64>, fun<nm::EW_ADD, float32_t, nm::Rational128>, NULL},																			\
																																																																																						\
			{fun<nm::EW_ADD, float64_t, uint8_t>, fun<nm::EW_ADD, float64_t, int8_t>, fun<nm::EW_ADD, float64_t, int16_t>, fun<nm::EW_ADD, float64_t, int32_t>, fun<nm::EW_ADD, float64_t, int64_t>,	\
				fun<nm::EW_ADD, float64_t, float32_t>, fun<nm::EW_ADD, float64_t, float64_t>, fun<nm::EW_ADD, float64_t, nm::Complex64>, fun<nm::EW_ADD, float64_t, nm::Complex128>,								\
				fun<nm::EW_ADD, float64_t, nm::Rational32>, fun<nm::EW_ADD, float64_t, nm::Rational64>, fun<nm::EW_ADD, float64_t, nm::Rational128>, NULL},																			\
																																																																																						\
			{fun<nm::EW_ADD, nm::Complex64, uint8_t>, fun<nm::EW_ADD, nm::Complex64, int8_t>, fun<nm::EW_ADD, nm::Complex64, int16_t>, fun<nm::EW_ADD, nm::Complex64, int32_t>,										\
				fun<nm::EW_ADD, nm::Complex64, int64_t>, fun<nm::EW_ADD, nm::Complex64, float32_t>, fun<nm::EW_ADD, nm::Complex64, float64_t>, fun<nm::EW_ADD, nm::Complex64, nm::Complex64>,				\
				fun<nm::EW_ADD, nm::Complex64, nm::Complex128>, fun<nm::EW_ADD, nm::Complex64, nm::Rational32>, fun<nm::EW_ADD, nm::Complex64, nm::Rational64>,																	\
				fun<nm::EW_ADD, nm::Complex64, nm::Rational128>, NULL},																																																									\
																																																																																						\
			{fun<nm::EW_ADD, nm::Complex128, uint8_t>, fun<nm::EW_ADD, nm::Complex128, int8_t>, fun<nm::EW_ADD, nm::Complex128, int16_t>, fun<nm::EW_ADD, nm::Complex128, int32_t>,								\
				fun<nm::EW_ADD, nm::Complex128, int64_t>, fun<nm::EW_ADD, nm::Complex128, float32_t>, fun<nm::EW_ADD, nm::Complex128, float64_t>, fun<nm::EW_ADD, nm::Complex128, nm::Complex64>,		\
				fun<nm::EW_ADD, nm::Complex128, nm::Complex128>,	fun<nm::EW_ADD, nm::Complex128, nm::Rational32>, fun<nm::EW_ADD, nm::Complex128, nm::Rational64>,															\
				fun<nm::EW_ADD, nm::Complex128, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_ADD, nm::Rational32, uint8_t>, fun<nm::EW_ADD, nm::Rational32, int8_t>, fun<nm::EW_ADD, nm::Rational32, int16_t>, fun<nm::EW_ADD, nm::Rational32, int32_t>,								\
				fun<nm::EW_ADD, nm::Rational32, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_ADD, nm::Rational32, nm::Rational32>, fun<nm::EW_ADD, nm::Rational32, nm::Rational64>,							\
				fun<nm::EW_ADD, nm::Rational32, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_ADD, nm::Rational64, uint8_t>, fun<nm::EW_ADD, nm::Rational64, int8_t>, fun<nm::EW_ADD, nm::Rational64, int16_t>, fun<nm::EW_ADD, nm::Rational64, int32_t>,								\
				fun<nm::EW_ADD, nm::Rational64, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_ADD, nm::Rational64, nm::Rational32>, fun<nm::EW_ADD, nm::Rational64, nm::Rational64>,							\
				fun<nm::EW_ADD, nm::Rational64, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_ADD, nm::Rational128, uint8_t>, fun<nm::EW_ADD, nm::Rational128, int8_t>, fun<nm::EW_ADD, nm::Rational128, int16_t>, fun<nm::EW_ADD, nm::Rational128, int32_t>,						\
				fun<nm::EW_ADD, nm::Rational128, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_ADD, nm::Rational128, nm::Rational32>, fun<nm::EW_ADD, nm::Rational128, nm::Rational64>,					\
				fun<nm::EW_ADD, nm::Rational128, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_ADD, nm::RubyObject, nm::RubyObject>}																									\
		},																																																																																			\
																																																																																						\
		{																																																																																				\
			{fun<nm::EW_SUB, uint8_t, uint8_t>, fun<nm::EW_SUB, uint8_t, int8_t>, fun<nm::EW_SUB, uint8_t, int16_t>, fun<nm::EW_SUB, uint8_t, int32_t>, fun<nm::EW_SUB, uint8_t, int64_t>,						\
				fun<nm::EW_SUB, uint8_t, float32_t>, fun<nm::EW_SUB, uint8_t, float64_t>, fun<nm::EW_SUB, uint8_t, nm::Complex64>, fun<nm::EW_SUB, uint8_t, nm::Complex128>,												\
				fun<nm::EW_SUB, uint8_t, nm::Rational32>, fun<nm::EW_SUB, uint8_t, nm::Rational64>, fun<nm::EW_SUB, uint8_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_SUB, int8_t, uint8_t>, fun<nm::EW_SUB, int8_t, int8_t>, fun<nm::EW_SUB, int8_t, int16_t>, fun<nm::EW_SUB, int8_t, int32_t>, fun<nm::EW_SUB, int8_t, int64_t>,									\
				fun<nm::EW_SUB, int8_t, float32_t>, fun<nm::EW_SUB, int8_t, float64_t>, fun<nm::EW_SUB, int8_t, nm::Complex64>, fun<nm::EW_SUB, int8_t, nm::Complex128>,														\
				fun<nm::EW_SUB, int8_t, nm::Rational32>, fun<nm::EW_SUB, int8_t, nm::Rational64>, fun<nm::EW_SUB, int8_t, nm::Rational128>, NULL},																							\
																																																																																						\
			{fun<nm::EW_SUB, int16_t, uint8_t>, fun<nm::EW_SUB, int16_t, int8_t>, fun<nm::EW_SUB, int16_t, int16_t>, fun<nm::EW_SUB, int16_t, int32_t>, fun<nm::EW_SUB, int16_t, int64_t>,						\
				fun<nm::EW_SUB, int16_t, float32_t>, fun<nm::EW_SUB, int16_t, float64_t>, fun<nm::EW_SUB, int16_t, nm::Complex64>, fun<nm::EW_SUB, int16_t, nm::Complex128>,												\
				fun<nm::EW_SUB, int16_t, nm::Rational32>, fun<nm::EW_SUB, int16_t, nm::Rational64>, fun<nm::EW_SUB, int16_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_SUB, int32_t, uint8_t>, fun<nm::EW_SUB, int32_t, int8_t>, fun<nm::EW_SUB, int32_t, int16_t>, fun<nm::EW_SUB, int32_t, int32_t>, fun<nm::EW_SUB, int32_t, int64_t>,						\
				fun<nm::EW_SUB, int32_t, float32_t>, fun<nm::EW_SUB, int32_t, float64_t>, fun<nm::EW_SUB, int32_t, nm::Complex64>, fun<nm::EW_SUB, int32_t, nm::Complex128>,												\
				fun<nm::EW_SUB, int32_t, nm::Rational32>, fun<nm::EW_SUB, int32_t, nm::Rational64>, fun<nm::EW_SUB, int32_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_SUB, int64_t, uint8_t>, fun<nm::EW_SUB, int64_t, int8_t>, fun<nm::EW_SUB, int64_t, int16_t>, fun<nm::EW_SUB, int64_t, int32_t>, fun<nm::EW_SUB, int64_t, int64_t>,						\
				fun<nm::EW_SUB, int64_t, float32_t>, fun<nm::EW_SUB, int64_t, float64_t>, fun<nm::EW_SUB, int64_t, nm::Complex64>, fun<nm::EW_SUB, int64_t, nm::Complex128>,												\
				fun<nm::EW_SUB, int64_t, nm::Rational32>, fun<nm::EW_SUB, int64_t, nm::Rational64>, fun<nm::EW_SUB, int64_t, nm::Rational128>, NULL}, 																					\
																																																																																						\
			{fun<nm::EW_SUB, float32_t, uint8_t>, fun<nm::EW_SUB, float32_t, int8_t>, fun<nm::EW_SUB, float32_t, int16_t>, fun<nm::EW_SUB, float32_t, int32_t>, fun<nm::EW_SUB, float32_t, int64_t>,	\
				fun<nm::EW_SUB, float32_t, float32_t>, fun<nm::EW_SUB, float32_t, float64_t>, fun<nm::EW_SUB, float32_t, nm::Complex64>, fun<nm::EW_SUB, float32_t, nm::Complex128>,								\
				fun<nm::EW_SUB, float32_t, nm::Rational32>, fun<nm::EW_SUB, float32_t, nm::Rational64>, fun<nm::EW_SUB, float32_t, nm::Rational128>, NULL},																			\
																																																																																						\
			{fun<nm::EW_SUB, float64_t, uint8_t>, fun<nm::EW_SUB, float64_t, int8_t>, fun<nm::EW_SUB, float64_t, int16_t>, fun<nm::EW_SUB, float64_t, int32_t>, fun<nm::EW_SUB, float64_t, int64_t>,	\
				fun<nm::EW_SUB, float64_t, float32_t>, fun<nm::EW_SUB, float64_t, float64_t>, fun<nm::EW_SUB, float64_t, nm::Complex64>, fun<nm::EW_SUB, float64_t, nm::Complex128>,								\
				fun<nm::EW_SUB, float64_t, nm::Rational32>, fun<nm::EW_SUB, float64_t, nm::Rational64>, fun<nm::EW_SUB, float64_t, nm::Rational128>, NULL},																			\
																																																																																						\
			{fun<nm::EW_SUB, nm::Complex64, uint8_t>, fun<nm::EW_SUB, nm::Complex64, int8_t>, fun<nm::EW_SUB, nm::Complex64, int16_t>, fun<nm::EW_SUB, nm::Complex64, int32_t>,										\
				fun<nm::EW_SUB, nm::Complex64, int64_t>, fun<nm::EW_SUB, nm::Complex64, float32_t>, fun<nm::EW_SUB, nm::Complex64, float64_t>, fun<nm::EW_SUB, nm::Complex64, nm::Complex64>,				\
				fun<nm::EW_SUB, nm::Complex64, nm::Complex128>, fun<nm::EW_SUB, nm::Complex64, nm::Rational32>, fun<nm::EW_SUB, nm::Complex64, nm::Rational64>,																	\
				fun<nm::EW_SUB, nm::Complex64, nm::Rational128>, NULL},																																																									\
																																																																																						\
			{fun<nm::EW_SUB, nm::Complex128, uint8_t>, fun<nm::EW_SUB, nm::Complex128, int8_t>, fun<nm::EW_SUB, nm::Complex128, int16_t>, fun<nm::EW_SUB, nm::Complex128, int32_t>,								\
				fun<nm::EW_SUB, nm::Complex128, int64_t>, fun<nm::EW_SUB, nm::Complex128, float32_t>, fun<nm::EW_SUB, nm::Complex128, float64_t>, fun<nm::EW_SUB, nm::Complex128, nm::Complex64>,		\
				fun<nm::EW_SUB, nm::Complex128, nm::Complex128>,	fun<nm::EW_SUB, nm::Complex128, nm::Rational32>, fun<nm::EW_SUB, nm::Complex128, nm::Rational64>,															\
				fun<nm::EW_SUB, nm::Complex128, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_SUB, nm::Rational32, uint8_t>, fun<nm::EW_SUB, nm::Rational32, int8_t>, fun<nm::EW_SUB, nm::Rational32, int16_t>, fun<nm::EW_SUB, nm::Rational32, int32_t>,								\
				fun<nm::EW_SUB, nm::Rational32, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_SUB, nm::Rational32, nm::Rational32>, fun<nm::EW_SUB, nm::Rational32, nm::Rational64>,							\
				fun<nm::EW_SUB, nm::Rational32, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_SUB, nm::Rational64, uint8_t>, fun<nm::EW_SUB, nm::Rational64, int8_t>, fun<nm::EW_SUB, nm::Rational64, int16_t>, fun<nm::EW_SUB, nm::Rational64, int32_t>,								\
				fun<nm::EW_SUB, nm::Rational64, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_SUB, nm::Rational64, nm::Rational32>, fun<nm::EW_SUB, nm::Rational64, nm::Rational64>,							\
				fun<nm::EW_SUB, nm::Rational64, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_SUB, nm::Rational128, uint8_t>, fun<nm::EW_SUB, nm::Rational128, int8_t>, fun<nm::EW_SUB, nm::Rational128, int16_t>, fun<nm::EW_SUB, nm::Rational128, int32_t>,						\
				fun<nm::EW_SUB, nm::Rational128, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_SUB, nm::Rational128, nm::Rational32>, fun<nm::EW_SUB, nm::Rational128, nm::Rational64>,					\
				fun<nm::EW_SUB, nm::Rational128, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_SUB, nm::RubyObject, nm::RubyObject>}																									\
		},																																																																																			\
																																																																																						\
		{																																																																																				\
			{fun<nm::EW_MUL, uint8_t, uint8_t>, fun<nm::EW_MUL, uint8_t, int8_t>, fun<nm::EW_MUL, uint8_t, int16_t>, fun<nm::EW_MUL, uint8_t, int32_t>, fun<nm::EW_MUL, uint8_t, int64_t>,						\
				fun<nm::EW_MUL, uint8_t, float32_t>, fun<nm::EW_MUL, uint8_t, float64_t>, fun<nm::EW_MUL, uint8_t, nm::Complex64>, fun<nm::EW_MUL, uint8_t, nm::Complex128>,												\
				fun<nm::EW_MUL, uint8_t, nm::Rational32>, fun<nm::EW_MUL, uint8_t, nm::Rational64>, fun<nm::EW_MUL, uint8_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_MUL, int8_t, uint8_t>, fun<nm::EW_MUL, int8_t, int8_t>, fun<nm::EW_MUL, int8_t, int16_t>, fun<nm::EW_MUL, int8_t, int32_t>, fun<nm::EW_MUL, int8_t, int64_t>,									\
				fun<nm::EW_MUL, int8_t, float32_t>, fun<nm::EW_MUL, int8_t, float64_t>, fun<nm::EW_MUL, int8_t, nm::Complex64>, fun<nm::EW_MUL, int8_t, nm::Complex128>,														\
				fun<nm::EW_MUL, int8_t, nm::Rational32>, fun<nm::EW_MUL, int8_t, nm::Rational64>, fun<nm::EW_MUL, int8_t, nm::Rational128>, NULL},																							\
																																																																																						\
			{fun<nm::EW_MUL, int16_t, uint8_t>, fun<nm::EW_MUL, int16_t, int8_t>, fun<nm::EW_MUL, int16_t, int16_t>, fun<nm::EW_MUL, int16_t, int32_t>, fun<nm::EW_MUL, int16_t, int64_t>,						\
				fun<nm::EW_MUL, int16_t, float32_t>, fun<nm::EW_MUL, int16_t, float64_t>, fun<nm::EW_MUL, int16_t, nm::Complex64>, fun<nm::EW_MUL, int16_t, nm::Complex128>,												\
				fun<nm::EW_MUL, int16_t, nm::Rational32>, fun<nm::EW_MUL, int16_t, nm::Rational64>, fun<nm::EW_MUL, int16_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_MUL, int32_t, uint8_t>, fun<nm::EW_MUL, int32_t, int8_t>, fun<nm::EW_MUL, int32_t, int16_t>, fun<nm::EW_MUL, int32_t, int32_t>, fun<nm::EW_MUL, int32_t, int64_t>,						\
				fun<nm::EW_MUL, int32_t, float32_t>, fun<nm::EW_MUL, int32_t, float64_t>, fun<nm::EW_MUL, int32_t, nm::Complex64>, fun<nm::EW_MUL, int32_t, nm::Complex128>,												\
				fun<nm::EW_MUL, int32_t, nm::Rational32>, fun<nm::EW_MUL, int32_t, nm::Rational64>, fun<nm::EW_MUL, int32_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_MUL, int64_t, uint8_t>, fun<nm::EW_MUL, int64_t, int8_t>, fun<nm::EW_MUL, int64_t, int16_t>, fun<nm::EW_MUL, int64_t, int32_t>, fun<nm::EW_MUL, int64_t, int64_t>,						\
				fun<nm::EW_MUL, int64_t, float32_t>, fun<nm::EW_MUL, int64_t, float64_t>, fun<nm::EW_MUL, int64_t, nm::Complex64>, fun<nm::EW_MUL, int64_t, nm::Complex128>,												\
				fun<nm::EW_MUL, int64_t, nm::Rational32>, fun<nm::EW_MUL, int64_t, nm::Rational64>, fun<nm::EW_MUL, int64_t, nm::Rational128>, NULL}, 																					\
																																																																																						\
			{fun<nm::EW_MUL, float32_t, uint8_t>, fun<nm::EW_MUL, float32_t, int8_t>, fun<nm::EW_MUL, float32_t, int16_t>, fun<nm::EW_MUL, float32_t, int32_t>, fun<nm::EW_MUL, float32_t, int64_t>,	\
				fun<nm::EW_MUL, float32_t, float32_t>, fun<nm::EW_MUL, float32_t, float64_t>, fun<nm::EW_MUL, float32_t, nm::Complex64>, fun<nm::EW_MUL, float32_t, nm::Complex128>,								\
				fun<nm::EW_MUL, float32_t, nm::Rational32>, fun<nm::EW_MUL, float32_t, nm::Rational64>, fun<nm::EW_MUL, float32_t, nm::Rational128>, NULL},																			\
																																																																																						\
			{fun<nm::EW_MUL, float64_t, uint8_t>, fun<nm::EW_MUL, float64_t, int8_t>, fun<nm::EW_MUL, float64_t, int16_t>, fun<nm::EW_MUL, float64_t, int32_t>, fun<nm::EW_MUL, float64_t, int64_t>,	\
				fun<nm::EW_MUL, float64_t, float32_t>, fun<nm::EW_MUL, float64_t, float64_t>, fun<nm::EW_MUL, float64_t, nm::Complex64>, fun<nm::EW_MUL, float64_t, nm::Complex128>,								\
				fun<nm::EW_MUL, float64_t, nm::Rational32>, fun<nm::EW_MUL, float64_t, nm::Rational64>, fun<nm::EW_MUL, float64_t, nm::Rational128>, NULL},																			\
																																																																																						\
			{fun<nm::EW_MUL, nm::Complex64, uint8_t>, fun<nm::EW_MUL, nm::Complex64, int8_t>, fun<nm::EW_MUL, nm::Complex64, int16_t>, fun<nm::EW_MUL, nm::Complex64, int32_t>,										\
				fun<nm::EW_MUL, nm::Complex64, int64_t>, fun<nm::EW_MUL, nm::Complex64, float32_t>, fun<nm::EW_MUL, nm::Complex64, float64_t>, fun<nm::EW_MUL, nm::Complex64, nm::Complex64>,				\
				fun<nm::EW_MUL, nm::Complex64, nm::Complex128>, fun<nm::EW_MUL, nm::Complex64, nm::Rational32>, fun<nm::EW_MUL, nm::Complex64, nm::Rational64>,																	\
				fun<nm::EW_MUL, nm::Complex64, nm::Rational128>, NULL},																																																									\
																																																																																						\
			{fun<nm::EW_MUL, nm::Complex128, uint8_t>, fun<nm::EW_MUL, nm::Complex128, int8_t>, fun<nm::EW_MUL, nm::Complex128, int16_t>, fun<nm::EW_MUL, nm::Complex128, int32_t>,								\
				fun<nm::EW_MUL, nm::Complex128, int64_t>, fun<nm::EW_MUL, nm::Complex128, float32_t>, fun<nm::EW_MUL, nm::Complex128, float64_t>, fun<nm::EW_MUL, nm::Complex128, nm::Complex64>,		\
				fun<nm::EW_MUL, nm::Complex128, nm::Complex128>,	fun<nm::EW_MUL, nm::Complex128, nm::Rational32>, fun<nm::EW_MUL, nm::Complex128, nm::Rational64>,															\
				fun<nm::EW_MUL, nm::Complex128, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_MUL, nm::Rational32, uint8_t>, fun<nm::EW_MUL, nm::Rational32, int8_t>, fun<nm::EW_MUL, nm::Rational32, int16_t>, fun<nm::EW_MUL, nm::Rational32, int32_t>,								\
				fun<nm::EW_MUL, nm::Rational32, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_MUL, nm::Rational32, nm::Rational32>, fun<nm::EW_MUL, nm::Rational32, nm::Rational64>,							\
				fun<nm::EW_MUL, nm::Rational32, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_MUL, nm::Rational64, uint8_t>, fun<nm::EW_MUL, nm::Rational64, int8_t>, fun<nm::EW_MUL, nm::Rational64, int16_t>, fun<nm::EW_MUL, nm::Rational64, int32_t>,								\
				fun<nm::EW_MUL, nm::Rational64, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_MUL, nm::Rational64, nm::Rational32>, fun<nm::EW_MUL, nm::Rational64, nm::Rational64>,							\
				fun<nm::EW_MUL, nm::Rational64, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_MUL, nm::Rational128, uint8_t>, fun<nm::EW_MUL, nm::Rational128, int8_t>, fun<nm::EW_MUL, nm::Rational128, int16_t>, fun<nm::EW_MUL, nm::Rational128, int32_t>,						\
				fun<nm::EW_MUL, nm::Rational128, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_MUL, nm::Rational128, nm::Rational32>, fun<nm::EW_MUL, nm::Rational128, nm::Rational64>,					\
				fun<nm::EW_MUL, nm::Rational128, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_MUL, nm::RubyObject, nm::RubyObject>}																									\
		},																																																																																			\
																																																																																						\
		{																																																																																				\
			{fun<nm::EW_DIV, uint8_t, uint8_t>, fun<nm::EW_DIV, uint8_t, int8_t>, fun<nm::EW_DIV, uint8_t, int16_t>, fun<nm::EW_DIV, uint8_t, int32_t>, fun<nm::EW_DIV, uint8_t, int64_t>,						\
				fun<nm::EW_DIV, uint8_t, float32_t>, fun<nm::EW_DIV, uint8_t, float64_t>, fun<nm::EW_DIV, uint8_t, nm::Complex64>, fun<nm::EW_DIV, uint8_t, nm::Complex128>,												\
				fun<nm::EW_DIV, uint8_t, nm::Rational32>, fun<nm::EW_DIV, uint8_t, nm::Rational64>, fun<nm::EW_DIV, uint8_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_DIV, int8_t, uint8_t>, fun<nm::EW_DIV, int8_t, int8_t>, fun<nm::EW_DIV, int8_t, int16_t>, fun<nm::EW_DIV, int8_t, int32_t>, fun<nm::EW_DIV, int8_t, int64_t>,									\
				fun<nm::EW_DIV, int8_t, float32_t>, fun<nm::EW_DIV, int8_t, float64_t>, fun<nm::EW_DIV, int8_t, nm::Complex64>, fun<nm::EW_DIV, int8_t, nm::Complex128>,														\
				fun<nm::EW_DIV, int8_t, nm::Rational32>, fun<nm::EW_DIV, int8_t, nm::Rational64>, fun<nm::EW_DIV, int8_t, nm::Rational128>, NULL},																							\
																																																																																						\
			{fun<nm::EW_DIV, int16_t, uint8_t>, fun<nm::EW_DIV, int16_t, int8_t>, fun<nm::EW_DIV, int16_t, int16_t>, fun<nm::EW_DIV, int16_t, int32_t>, fun<nm::EW_DIV, int16_t, int64_t>,						\
				fun<nm::EW_DIV, int16_t, float32_t>, fun<nm::EW_DIV, int16_t, float64_t>, fun<nm::EW_DIV, int16_t, nm::Complex64>, fun<nm::EW_DIV, int16_t, nm::Complex128>,												\
				fun<nm::EW_DIV, int16_t, nm::Rational32>, fun<nm::EW_DIV, int16_t, nm::Rational64>, fun<nm::EW_DIV, int16_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_DIV, int32_t, uint8_t>, fun<nm::EW_DIV, int32_t, int8_t>, fun<nm::EW_DIV, int32_t, int16_t>, fun<nm::EW_DIV, int32_t, int32_t>, fun<nm::EW_DIV, int32_t, int64_t>,						\
				fun<nm::EW_DIV, int32_t, float32_t>, fun<nm::EW_DIV, int32_t, float64_t>, fun<nm::EW_DIV, int32_t, nm::Complex64>, fun<nm::EW_DIV, int32_t, nm::Complex128>,												\
				fun<nm::EW_DIV, int32_t, nm::Rational32>, fun<nm::EW_DIV, int32_t, nm::Rational64>, fun<nm::EW_DIV, int32_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_DIV, int64_t, uint8_t>, fun<nm::EW_DIV, int64_t, int8_t>, fun<nm::EW_DIV, int64_t, int16_t>, fun<nm::EW_DIV, int64_t, int32_t>, fun<nm::EW_DIV, int64_t, int64_t>,						\
				fun<nm::EW_DIV, int64_t, float32_t>, fun<nm::EW_DIV, int64_t, float64_t>, fun<nm::EW_DIV, int64_t, nm::Complex64>, fun<nm::EW_DIV, int64_t, nm::Complex128>,												\
				fun<nm::EW_DIV, int64_t, nm::Rational32>, fun<nm::EW_DIV, int64_t, nm::Rational64>, fun<nm::EW_DIV, int64_t, nm::Rational128>, NULL}, 																					\
																																																																																						\
			{fun<nm::EW_DIV, float32_t, uint8_t>, fun<nm::EW_DIV, float32_t, int8_t>, fun<nm::EW_DIV, float32_t, int16_t>, fun<nm::EW_DIV, float32_t, int32_t>, fun<nm::EW_DIV, float32_t, int64_t>,	\
				fun<nm::EW_DIV, float32_t, float32_t>, fun<nm::EW_DIV, float32_t, float64_t>, fun<nm::EW_DIV, float32_t, nm::Complex64>, fun<nm::EW_DIV, float32_t, nm::Complex128>,								\
				fun<nm::EW_DIV, float32_t, nm::Rational32>, fun<nm::EW_DIV, float32_t, nm::Rational64>, fun<nm::EW_DIV, float32_t, nm::Rational128>, NULL},																			\
																																																																																						\
			{fun<nm::EW_DIV, float64_t, uint8_t>, fun<nm::EW_DIV, float64_t, int8_t>, fun<nm::EW_DIV, float64_t, int16_t>, fun<nm::EW_DIV, float64_t, int32_t>, fun<nm::EW_DIV, float64_t, int64_t>,	\
				fun<nm::EW_DIV, float64_t, float32_t>, fun<nm::EW_DIV, float64_t, float64_t>, fun<nm::EW_DIV, float64_t, nm::Complex64>, fun<nm::EW_DIV, float64_t, nm::Complex128>,								\
				fun<nm::EW_DIV, float64_t, nm::Rational32>, fun<nm::EW_DIV, float64_t, nm::Rational64>, fun<nm::EW_DIV, float64_t, nm::Rational128>, NULL},																			\
																																																																																						\
			{fun<nm::EW_DIV, nm::Complex64, uint8_t>, fun<nm::EW_DIV, nm::Complex64, int8_t>, fun<nm::EW_DIV, nm::Complex64, int16_t>, fun<nm::EW_DIV, nm::Complex64, int32_t>,										\
				fun<nm::EW_DIV, nm::Complex64, int64_t>, fun<nm::EW_DIV, nm::Complex64, float32_t>, fun<nm::EW_DIV, nm::Complex64, float64_t>, fun<nm::EW_DIV, nm::Complex64, nm::Complex64>,				\
				fun<nm::EW_DIV, nm::Complex64, nm::Complex128>, fun<nm::EW_DIV, nm::Complex64, nm::Rational32>, fun<nm::EW_DIV, nm::Complex64, nm::Rational64>,																	\
				fun<nm::EW_DIV, nm::Complex64, nm::Rational128>, NULL},																																																									\
																																																																																						\
			{fun<nm::EW_DIV, nm::Complex128, uint8_t>, fun<nm::EW_DIV, nm::Complex128, int8_t>, fun<nm::EW_DIV, nm::Complex128, int16_t>, fun<nm::EW_DIV, nm::Complex128, int32_t>,								\
				fun<nm::EW_DIV, nm::Complex128, int64_t>, fun<nm::EW_DIV, nm::Complex128, float32_t>, fun<nm::EW_DIV, nm::Complex128, float64_t>, fun<nm::EW_DIV, nm::Complex128, nm::Complex64>,		\
				fun<nm::EW_DIV, nm::Complex128, nm::Complex128>,	fun<nm::EW_DIV, nm::Complex128, nm::Rational32>, fun<nm::EW_DIV, nm::Complex128, nm::Rational64>,															\
				fun<nm::EW_DIV, nm::Complex128, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_DIV, nm::Rational32, uint8_t>, fun<nm::EW_DIV, nm::Rational32, int8_t>, fun<nm::EW_DIV, nm::Rational32, int16_t>, fun<nm::EW_DIV, nm::Rational32, int32_t>,								\
				fun<nm::EW_DIV, nm::Rational32, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_DIV, nm::Rational32, nm::Rational32>, fun<nm::EW_DIV, nm::Rational32, nm::Rational64>,							\
				fun<nm::EW_DIV, nm::Rational32, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_DIV, nm::Rational64, uint8_t>, fun<nm::EW_DIV, nm::Rational64, int8_t>, fun<nm::EW_DIV, nm::Rational64, int16_t>, fun<nm::EW_DIV, nm::Rational64, int32_t>,								\
				fun<nm::EW_DIV, nm::Rational64, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_DIV, nm::Rational64, nm::Rational32>, fun<nm::EW_DIV, nm::Rational64, nm::Rational64>,							\
				fun<nm::EW_DIV, nm::Rational64, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_DIV, nm::Rational128, uint8_t>, fun<nm::EW_DIV, nm::Rational128, int8_t>, fun<nm::EW_DIV, nm::Rational128, int16_t>, fun<nm::EW_DIV, nm::Rational128, int32_t>,						\
				fun<nm::EW_DIV, nm::Rational128, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_DIV, nm::Rational128, nm::Rational32>, fun<nm::EW_DIV, nm::Rational128, nm::Rational64>,					\
				fun<nm::EW_DIV, nm::Rational128, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_DIV, nm::RubyObject, nm::RubyObject>}																									\
		},																																																																																			\
																																																																																						\
		{																																																																																				\
			{fun<nm::EW_MOD, uint8_t, uint8_t>, fun<nm::EW_MOD, uint8_t, int8_t>, fun<nm::EW_MOD, uint8_t, int16_t>, fun<nm::EW_MOD, uint8_t, int32_t>, fun<nm::EW_MOD, uint8_t, int64_t>,						\
				fun<nm::EW_MOD, uint8_t, float32_t>, fun<nm::EW_MOD, uint8_t, float64_t>, fun<nm::EW_MOD, uint8_t, nm::Complex64>, fun<nm::EW_MOD, uint8_t, nm::Complex128>,												\
				fun<nm::EW_MOD, uint8_t, nm::Rational32>, fun<nm::EW_MOD, uint8_t, nm::Rational64>, fun<nm::EW_MOD, uint8_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_MOD, int8_t, uint8_t>, fun<nm::EW_MOD, int8_t, int8_t>, fun<nm::EW_MOD, int8_t, int16_t>, fun<nm::EW_MOD, int8_t, int32_t>, fun<nm::EW_MOD, int8_t, int64_t>,									\
				fun<nm::EW_MOD, int8_t, float32_t>, fun<nm::EW_MOD, int8_t, float64_t>, fun<nm::EW_MOD, int8_t, nm::Complex64>, fun<nm::EW_MOD, int8_t, nm::Complex128>,														\
				fun<nm::EW_MOD, int8_t, nm::Rational32>, fun<nm::EW_MOD, int8_t, nm::Rational64>, fun<nm::EW_MOD, int8_t, nm::Rational128>, NULL},																							\
																																																																																						\
			{fun<nm::EW_MOD, int16_t, uint8_t>, fun<nm::EW_MOD, int16_t, int8_t>, fun<nm::EW_MOD, int16_t, int16_t>, fun<nm::EW_MOD, int16_t, int32_t>, fun<nm::EW_MOD, int16_t, int64_t>,						\
				fun<nm::EW_MOD, int16_t, float32_t>, fun<nm::EW_MOD, int16_t, float64_t>, fun<nm::EW_MOD, int16_t, nm::Complex64>, fun<nm::EW_MOD, int16_t, nm::Complex128>,												\
				fun<nm::EW_MOD, int16_t, nm::Rational32>, fun<nm::EW_MOD, int16_t, nm::Rational64>, fun<nm::EW_MOD, int16_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_MOD, int32_t, uint8_t>, fun<nm::EW_MOD, int32_t, int8_t>, fun<nm::EW_MOD, int32_t, int16_t>, fun<nm::EW_MOD, int32_t, int32_t>, fun<nm::EW_MOD, int32_t, int64_t>,						\
				fun<nm::EW_MOD, int32_t, float32_t>, fun<nm::EW_MOD, int32_t, float64_t>, fun<nm::EW_MOD, int32_t, nm::Complex64>, fun<nm::EW_MOD, int32_t, nm::Complex128>,												\
				fun<nm::EW_MOD, int32_t, nm::Rational32>, fun<nm::EW_MOD, int32_t, nm::Rational64>, fun<nm::EW_MOD, int32_t, nm::Rational128>, NULL},																						\
																																																																																						\
			{fun<nm::EW_MOD, int64_t, uint8_t>, fun<nm::EW_MOD, int64_t, int8_t>, fun<nm::EW_MOD, int64_t, int16_t>, fun<nm::EW_MOD, int64_t, int32_t>, fun<nm::EW_MOD, int64_t, int64_t>,						\
				fun<nm::EW_MOD, int64_t, float32_t>, fun<nm::EW_MOD, int64_t, float64_t>, fun<nm::EW_MOD, int64_t, nm::Complex64>, fun<nm::EW_MOD, int64_t, nm::Complex128>,												\
				fun<nm::EW_MOD, int64_t, nm::Rational32>, fun<nm::EW_MOD, int64_t, nm::Rational64>, fun<nm::EW_MOD, int64_t, nm::Rational128>, NULL}, 																					\
																																																																																						\
			{fun<nm::EW_MOD, float32_t, uint8_t>, fun<nm::EW_MOD, float32_t, int8_t>, fun<nm::EW_MOD, float32_t, int16_t>, fun<nm::EW_MOD, float32_t, int32_t>, fun<nm::EW_MOD, float32_t, int64_t>,	\
				fun<nm::EW_MOD, float32_t, float32_t>, fun<nm::EW_MOD, float32_t, float64_t>, fun<nm::EW_MOD, float32_t, nm::Complex64>, fun<nm::EW_MOD, float32_t, nm::Complex128>,								\
				fun<nm::EW_MOD, float32_t, nm::Rational32>, fun<nm::EW_MOD, float32_t, nm::Rational64>, fun<nm::EW_MOD, float32_t, nm::Rational128>, NULL},																			\
																																																																																						\
			{fun<nm::EW_MOD, float64_t, uint8_t>, fun<nm::EW_MOD, float64_t, int8_t>, fun<nm::EW_MOD, float64_t, int16_t>, fun<nm::EW_MOD, float64_t, int32_t>, fun<nm::EW_MOD, float64_t, int64_t>,	\
				fun<nm::EW_MOD, float64_t, float32_t>, fun<nm::EW_MOD, float64_t, float64_t>, fun<nm::EW_MOD, float64_t, nm::Complex64>, fun<nm::EW_MOD, float64_t, nm::Complex128>,								\
				fun<nm::EW_MOD, float64_t, nm::Rational32>, fun<nm::EW_MOD, float64_t, nm::Rational64>, fun<nm::EW_MOD, float64_t, nm::Rational128>, NULL},																			\
																																																																																						\
			{fun<nm::EW_MOD, nm::Complex64, uint8_t>, fun<nm::EW_MOD, nm::Complex64, int8_t>, fun<nm::EW_MOD, nm::Complex64, int16_t>, fun<nm::EW_MOD, nm::Complex64, int32_t>,										\
				fun<nm::EW_MOD, nm::Complex64, int64_t>, fun<nm::EW_MOD, nm::Complex64, float32_t>, fun<nm::EW_MOD, nm::Complex64, float64_t>, fun<nm::EW_MOD, nm::Complex64, nm::Complex64>,				\
				fun<nm::EW_MOD, nm::Complex64, nm::Complex128>, fun<nm::EW_MOD, nm::Complex64, nm::Rational32>, fun<nm::EW_MOD, nm::Complex64, nm::Rational64>,																	\
				fun<nm::EW_MOD, nm::Complex64, nm::Rational128>, NULL},																																																									\
																																																																																						\
			{fun<nm::EW_MOD, nm::Complex128, uint8_t>, fun<nm::EW_MOD, nm::Complex128, int8_t>, fun<nm::EW_MOD, nm::Complex128, int16_t>, fun<nm::EW_MOD, nm::Complex128, int32_t>,								\
				fun<nm::EW_MOD, nm::Complex128, int64_t>, fun<nm::EW_MOD, nm::Complex128, float32_t>, fun<nm::EW_MOD, nm::Complex128, float64_t>, fun<nm::EW_MOD, nm::Complex128, nm::Complex64>,		\
				fun<nm::EW_MOD, nm::Complex128, nm::Complex128>,	fun<nm::EW_MOD, nm::Complex128, nm::Rational32>, fun<nm::EW_MOD, nm::Complex128, nm::Rational64>,															\
				fun<nm::EW_MOD, nm::Complex128, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_MOD, nm::Rational32, uint8_t>, fun<nm::EW_MOD, nm::Rational32, int8_t>, fun<nm::EW_MOD, nm::Rational32, int16_t>, fun<nm::EW_MOD, nm::Rational32, int32_t>,								\
				fun<nm::EW_MOD, nm::Rational32, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_MOD, nm::Rational32, nm::Rational32>, fun<nm::EW_MOD, nm::Rational32, nm::Rational64>,							\
				fun<nm::EW_MOD, nm::Rational32, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_MOD, nm::Rational64, uint8_t>, fun<nm::EW_MOD, nm::Rational64, int8_t>, fun<nm::EW_MOD, nm::Rational64, int16_t>, fun<nm::EW_MOD, nm::Rational64, int32_t>,								\
				fun<nm::EW_MOD, nm::Rational64, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_MOD, nm::Rational64, nm::Rational32>, fun<nm::EW_MOD, nm::Rational64, nm::Rational64>,							\
				fun<nm::EW_MOD, nm::Rational64, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{fun<nm::EW_MOD, nm::Rational128, uint8_t>, fun<nm::EW_MOD, nm::Rational128, int8_t>, fun<nm::EW_MOD, nm::Rational128, int16_t>, fun<nm::EW_MOD, nm::Rational128, int32_t>,						\
				fun<nm::EW_MOD, nm::Rational128, int64_t>, NULL, NULL, NULL, NULL, fun<nm::EW_MOD, nm::Rational128, nm::Rational32>, fun<nm::EW_MOD, nm::Rational128, nm::Rational64>,					\
				fun<nm::EW_MOD, nm::Rational128, nm::Rational128>, NULL},																																																								\
																																																																																						\
			{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, fun<nm::EW_MOD, nm::RubyObject, nm::RubyObject>}																									\
		}																																																																																				\
	};

/*
 * Defines a static array that holds function pointers to left dtype, right
 * dtype, and itype templated versions of the specified function.
 */
#define LRI_DTYPE_TEMPLATE_TABLE(fun, ret, ...) NAMED_LRI_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)

#define NAMED_LRI_DTYPE_TEMPLATE_TABLE(name,  fun,  ret,  ...) \
static ret (*(name)[nm::NUM_DTYPES][nm::NUM_DTYPES][nm::NUM_ITYPES])(__VA_ARGS__) = { \
  { \
    {fun<uint8_t, uint8_t, uint8_t>, fun<uint8_t, uint8_t, uint16_t>, fun<uint8_t, uint8_t, uint32_t>, fun<uint8_t, uint8_t, uint64_t> },  \
    {fun<uint8_t, int8_t, uint8_t>, fun<uint8_t, int8_t, uint16_t>, fun<uint8_t, int8_t, uint32_t>, fun<uint8_t, int8_t, uint64_t> },  \
    {fun<uint8_t, int16_t, uint8_t>, fun<uint8_t, int16_t, uint16_t>, fun<uint8_t, int16_t, uint32_t>, fun<uint8_t, int16_t, uint64_t> },  \
    {fun<uint8_t, int32_t, uint8_t>, fun<uint8_t, int32_t, uint16_t>, fun<uint8_t, int32_t, uint32_t>, fun<uint8_t, int32_t, uint64_t> },  \
    {fun<uint8_t, int64_t, uint8_t>, fun<uint8_t, int64_t, uint16_t>, fun<uint8_t, int64_t, uint32_t>, fun<uint8_t, int64_t, uint64_t> },  \
    {fun<uint8_t, float32_t, uint8_t>, fun<uint8_t, float32_t, uint16_t>, fun<uint8_t, float32_t, uint32_t>, fun<uint8_t, float32_t, uint64_t> },  \
    {fun<uint8_t, float64_t, uint8_t>, fun<uint8_t, float64_t, uint16_t>, fun<uint8_t, float64_t, uint32_t>, fun<uint8_t, float64_t, uint64_t> },  \
    {fun<uint8_t, nm::Complex64, uint8_t>, fun<uint8_t, nm::Complex64, uint16_t>, fun<uint8_t, nm::Complex64, uint32_t>, fun<uint8_t, nm::Complex64, uint64_t> },  \
    {fun<uint8_t, nm::Complex128, uint8_t>, fun<uint8_t, nm::Complex128, uint16_t>, fun<uint8_t, nm::Complex128, uint32_t>, fun<uint8_t, nm::Complex128, uint64_t> },  \
    {fun<uint8_t, nm::Rational32, uint8_t>, fun<uint8_t, nm::Rational32, uint16_t>, fun<uint8_t, nm::Rational32, uint32_t>, fun<uint8_t, nm::Rational32, uint64_t> },  \
    {fun<uint8_t, nm::Rational64, uint8_t>, fun<uint8_t, nm::Rational64, uint16_t>, fun<uint8_t, nm::Rational64, uint32_t>, fun<uint8_t, nm::Rational64, uint64_t> },  \
    {fun<uint8_t, nm::Rational128, uint8_t>, fun<uint8_t, nm::Rational128, uint16_t>, fun<uint8_t, nm::Rational128, uint32_t>, fun<uint8_t, nm::Rational128, uint64_t> },  \
    {fun<uint8_t, nm::RubyObject, uint8_t>, fun<uint8_t, nm::RubyObject, uint16_t>, fun<uint8_t, nm::RubyObject, uint32_t>, fun<uint8_t, nm::RubyObject, uint64_t> } \
  },  \
  {   \
    {fun<int8_t, uint8_t, uint8_t>, fun<int8_t, uint8_t, uint16_t>, fun<int8_t, uint8_t, uint32_t>, fun<int8_t, uint8_t, uint64_t> },  \
    {fun<int8_t, int8_t, uint8_t>, fun<int8_t, int8_t, uint16_t>, fun<int8_t, int8_t, uint32_t>, fun<int8_t, int8_t, uint64_t> },  \
    {fun<int8_t, int16_t, uint8_t>, fun<int8_t, int16_t, uint16_t>, fun<int8_t, int16_t, uint32_t>, fun<int8_t, int16_t, uint64_t> },  \
    {fun<int8_t, int32_t, uint8_t>, fun<int8_t, int32_t, uint16_t>, fun<int8_t, int32_t, uint32_t>, fun<int8_t, int32_t, uint64_t> },  \
    {fun<int8_t, int64_t, uint8_t>, fun<int8_t, int64_t, uint16_t>, fun<int8_t, int64_t, uint32_t>, fun<int8_t, int64_t, uint64_t> },  \
    {fun<int8_t, float32_t, uint8_t>, fun<int8_t, float32_t, uint16_t>, fun<int8_t, float32_t, uint32_t>, fun<int8_t, float32_t, uint64_t> },  \
    {fun<int8_t, float64_t, uint8_t>, fun<int8_t, float64_t, uint16_t>, fun<int8_t, float64_t, uint32_t>, fun<int8_t, float64_t, uint64_t> },  \
    {fun<int8_t, nm::Complex64, uint8_t>, fun<int8_t, nm::Complex64, uint16_t>, fun<int8_t, nm::Complex64, uint32_t>, fun<int8_t, nm::Complex64, uint64_t> },  \
    {fun<int8_t, nm::Complex128, uint8_t>, fun<int8_t, nm::Complex128, uint16_t>, fun<int8_t, nm::Complex128, uint32_t>, fun<int8_t, nm::Complex128, uint64_t> },  \
    {fun<int8_t, nm::Rational32, uint8_t>, fun<int8_t, nm::Rational32, uint16_t>, fun<int8_t, nm::Rational32, uint32_t>, fun<int8_t, nm::Rational32, uint64_t> },  \
    {fun<int8_t, nm::Rational64, uint8_t>, fun<int8_t, nm::Rational64, uint16_t>, fun<int8_t, nm::Rational64, uint32_t>, fun<int8_t, nm::Rational64, uint64_t> },  \
    {fun<int8_t, nm::Rational128, uint8_t>, fun<int8_t, nm::Rational128, uint16_t>, fun<int8_t, nm::Rational128, uint32_t>, fun<int8_t, nm::Rational128, uint64_t> },  \
    {fun<int8_t, nm::RubyObject, uint8_t>, fun<int8_t, nm::RubyObject, uint16_t>, fun<int8_t, nm::RubyObject, uint32_t>, fun<int8_t, nm::RubyObject, uint64_t> } \
  },  \
  {{fun<int16_t, uint8_t, uint8_t>, fun<int16_t, uint8_t, uint16_t>, fun<int16_t, uint8_t, uint32_t>, fun<int16_t, uint8_t, uint64_t> }, {fun<int16_t, int8_t, uint8_t>, fun<int16_t, int8_t, uint16_t>, fun<int16_t, int8_t, uint32_t>, fun<int16_t, int8_t, uint64_t> }, {fun<int16_t, int16_t, uint8_t>, fun<int16_t, int16_t, uint16_t>, fun<int16_t, int16_t, uint32_t>, fun<int16_t, int16_t, uint64_t> }, {fun<int16_t, int32_t, uint8_t>, fun<int16_t, int32_t, uint16_t>, fun<int16_t, int32_t, uint32_t>, fun<int16_t, int32_t, uint64_t> }, {fun<int16_t, int64_t, uint8_t>, fun<int16_t, int64_t, uint16_t>, fun<int16_t, int64_t, uint32_t>, fun<int16_t, int64_t, uint64_t> }, {fun<int16_t, float32_t, uint8_t>, fun<int16_t, float32_t, uint16_t>, fun<int16_t, float32_t, uint32_t>, fun<int16_t, float32_t, uint64_t> }, {fun<int16_t, float64_t, uint8_t>, fun<int16_t, float64_t, uint16_t>, fun<int16_t, float64_t, uint32_t>, fun<int16_t, float64_t, uint64_t> }, {fun<int16_t, nm::Complex64, uint8_t>, fun<int16_t, nm::Complex64, uint16_t>, fun<int16_t, nm::Complex64, uint32_t>, fun<int16_t, nm::Complex64, uint64_t> }, {fun<int16_t, nm::Complex128, uint8_t>, fun<int16_t, nm::Complex128, uint16_t>, fun<int16_t, nm::Complex128, uint32_t>, fun<int16_t, nm::Complex128, uint64_t> }, {fun<int16_t, nm::Rational32, uint8_t>, fun<int16_t, nm::Rational32, uint16_t>, fun<int16_t, nm::Rational32, uint32_t>, fun<int16_t, nm::Rational32, uint64_t> }, {fun<int16_t, nm::Rational64, uint8_t>, fun<int16_t, nm::Rational64, uint16_t>, fun<int16_t, nm::Rational64, uint32_t>, fun<int16_t, nm::Rational64, uint64_t> }, {fun<int16_t, nm::Rational128, uint8_t>, fun<int16_t, nm::Rational128, uint16_t>, fun<int16_t, nm::Rational128, uint32_t>, fun<int16_t, nm::Rational128, uint64_t> }, {fun<int16_t, nm::RubyObject, uint8_t>, fun<int16_t, nm::RubyObject, uint16_t>, fun<int16_t, nm::RubyObject, uint32_t>, fun<int16_t, nm::RubyObject, uint64_t> }},  \
  {{fun<int32_t, uint8_t, uint8_t>, fun<int32_t, uint8_t, uint16_t>, fun<int32_t, uint8_t, uint32_t>, fun<int32_t, uint8_t, uint64_t> }, {fun<int32_t, int8_t, uint8_t>, fun<int32_t, int8_t, uint16_t>, fun<int32_t, int8_t, uint32_t>, fun<int32_t, int8_t, uint64_t> }, {fun<int32_t, int16_t, uint8_t>, fun<int32_t, int16_t, uint16_t>, fun<int32_t, int16_t, uint32_t>, fun<int32_t, int16_t, uint64_t> }, {fun<int32_t, int32_t, uint8_t>, fun<int32_t, int32_t, uint16_t>, fun<int32_t, int32_t, uint32_t>, fun<int32_t, int32_t, uint64_t> }, {fun<int32_t, int64_t, uint8_t>, fun<int32_t, int64_t, uint16_t>, fun<int32_t, int64_t, uint32_t>, fun<int32_t, int64_t, uint64_t> }, {fun<int32_t, float32_t, uint8_t>, fun<int32_t, float32_t, uint16_t>, fun<int32_t, float32_t, uint32_t>, fun<int32_t, float32_t, uint64_t> }, {fun<int32_t, float64_t, uint8_t>, fun<int32_t, float64_t, uint16_t>, fun<int32_t, float64_t, uint32_t>, fun<int32_t, float64_t, uint64_t> }, {fun<int32_t, nm::Complex64, uint8_t>, fun<int32_t, nm::Complex64, uint16_t>, fun<int32_t, nm::Complex64, uint32_t>, fun<int32_t, nm::Complex64, uint64_t> }, {fun<int32_t, nm::Complex128, uint8_t>, fun<int32_t, nm::Complex128, uint16_t>, fun<int32_t, nm::Complex128, uint32_t>, fun<int32_t, nm::Complex128, uint64_t> }, {fun<int32_t, nm::Rational32, uint8_t>, fun<int32_t, nm::Rational32, uint16_t>, fun<int32_t, nm::Rational32, uint32_t>, fun<int32_t, nm::Rational32, uint64_t> }, {fun<int32_t, nm::Rational64, uint8_t>, fun<int32_t, nm::Rational64, uint16_t>, fun<int32_t, nm::Rational64, uint32_t>, fun<int32_t, nm::Rational64, uint64_t> }, {fun<int32_t, nm::Rational128, uint8_t>, fun<int32_t, nm::Rational128, uint16_t>, fun<int32_t, nm::Rational128, uint32_t>, fun<int32_t, nm::Rational128, uint64_t> }, {fun<int32_t, nm::RubyObject, uint8_t>, fun<int32_t, nm::RubyObject, uint16_t>, fun<int32_t, nm::RubyObject, uint32_t>, fun<int32_t, nm::RubyObject, uint64_t> }},  \
  {{fun<int64_t, uint8_t, uint8_t>, fun<int64_t, uint8_t, uint16_t>, fun<int64_t, uint8_t, uint32_t>, fun<int64_t, uint8_t, uint64_t> }, {fun<int64_t, int8_t, uint8_t>, fun<int64_t, int8_t, uint16_t>, fun<int64_t, int8_t, uint32_t>, fun<int64_t, int8_t, uint64_t> }, {fun<int64_t, int16_t, uint8_t>, fun<int64_t, int16_t, uint16_t>, fun<int64_t, int16_t, uint32_t>, fun<int64_t, int16_t, uint64_t> }, {fun<int64_t, int32_t, uint8_t>, fun<int64_t, int32_t, uint16_t>, fun<int64_t, int32_t, uint32_t>, fun<int64_t, int32_t, uint64_t> }, {fun<int64_t, int64_t, uint8_t>, fun<int64_t, int64_t, uint16_t>, fun<int64_t, int64_t, uint32_t>, fun<int64_t, int64_t, uint64_t> }, {fun<int64_t, float32_t, uint8_t>, fun<int64_t, float32_t, uint16_t>, fun<int64_t, float32_t, uint32_t>, fun<int64_t, float32_t, uint64_t> }, {fun<int64_t, float64_t, uint8_t>, fun<int64_t, float64_t, uint16_t>, fun<int64_t, float64_t, uint32_t>, fun<int64_t, float64_t, uint64_t> }, {fun<int64_t, nm::Complex64, uint8_t>, fun<int64_t, nm::Complex64, uint16_t>, fun<int64_t, nm::Complex64, uint32_t>, fun<int64_t, nm::Complex64, uint64_t> }, {fun<int64_t, nm::Complex128, uint8_t>, fun<int64_t, nm::Complex128, uint16_t>, fun<int64_t, nm::Complex128, uint32_t>, fun<int64_t, nm::Complex128, uint64_t> }, {fun<int64_t, nm::Rational32, uint8_t>, fun<int64_t, nm::Rational32, uint16_t>, fun<int64_t, nm::Rational32, uint32_t>, fun<int64_t, nm::Rational32, uint64_t> }, {fun<int64_t, nm::Rational64, uint8_t>, fun<int64_t, nm::Rational64, uint16_t>, fun<int64_t, nm::Rational64, uint32_t>, fun<int64_t, nm::Rational64, uint64_t> }, {fun<int64_t, nm::Rational128, uint8_t>, fun<int64_t, nm::Rational128, uint16_t>, fun<int64_t, nm::Rational128, uint32_t>, fun<int64_t, nm::Rational128, uint64_t> }, {fun<int64_t, nm::RubyObject, uint8_t>, fun<int64_t, nm::RubyObject, uint16_t>, fun<int64_t, nm::RubyObject, uint32_t>, fun<int64_t, nm::RubyObject, uint64_t> }},  \
  {{fun<float32_t, uint8_t, uint8_t>, fun<float32_t, uint8_t, uint16_t>, fun<float32_t, uint8_t, uint32_t>, fun<float32_t, uint8_t, uint64_t> }, {fun<float32_t, int8_t, uint8_t>, fun<float32_t, int8_t, uint16_t>, fun<float32_t, int8_t, uint32_t>, fun<float32_t, int8_t, uint64_t> }, {fun<float32_t, int16_t, uint8_t>, fun<float32_t, int16_t, uint16_t>, fun<float32_t, int16_t, uint32_t>, fun<float32_t, int16_t, uint64_t> }, {fun<float32_t, int32_t, uint8_t>, fun<float32_t, int32_t, uint16_t>, fun<float32_t, int32_t, uint32_t>, fun<float32_t, int32_t, uint64_t> }, {fun<float32_t, int64_t, uint8_t>, fun<float32_t, int64_t, uint16_t>, fun<float32_t, int64_t, uint32_t>, fun<float32_t, int64_t, uint64_t> }, {fun<float32_t, float32_t, uint8_t>, fun<float32_t, float32_t, uint16_t>, fun<float32_t, float32_t, uint32_t>, fun<float32_t, float32_t, uint64_t> }, {fun<float32_t, float64_t, uint8_t>, fun<float32_t, float64_t, uint16_t>, fun<float32_t, float64_t, uint32_t>, fun<float32_t, float64_t, uint64_t> }, {fun<float32_t, nm::Complex64, uint8_t>, fun<float32_t, nm::Complex64, uint16_t>, fun<float32_t, nm::Complex64, uint32_t>, fun<float32_t, nm::Complex64, uint64_t> }, {fun<float32_t, nm::Complex128, uint8_t>, fun<float32_t, nm::Complex128, uint16_t>, fun<float32_t, nm::Complex128, uint32_t>, fun<float32_t, nm::Complex128, uint64_t> }, {fun<float32_t, nm::Rational32, uint8_t>, fun<float32_t, nm::Rational32, uint16_t>, fun<float32_t, nm::Rational32, uint32_t>, fun<float32_t, nm::Rational32, uint64_t> }, {fun<float32_t, nm::Rational64, uint8_t>, fun<float32_t, nm::Rational64, uint16_t>, fun<float32_t, nm::Rational64, uint32_t>, fun<float32_t, nm::Rational64, uint64_t> }, {fun<float32_t, nm::Rational128, uint8_t>, fun<float32_t, nm::Rational128, uint16_t>, fun<float32_t, nm::Rational128, uint32_t>, fun<float32_t, nm::Rational128, uint64_t> }, {fun<float32_t, nm::RubyObject, uint8_t>, fun<float32_t, nm::RubyObject, uint16_t>, fun<float32_t, nm::RubyObject, uint32_t>, fun<float32_t, nm::RubyObject, uint64_t> }},  \
  {{fun<float64_t, uint8_t, uint8_t>, fun<float64_t, uint8_t, uint16_t>, fun<float64_t, uint8_t, uint32_t>, fun<float64_t, uint8_t, uint64_t> }, {fun<float64_t, int8_t, uint8_t>, fun<float64_t, int8_t, uint16_t>, fun<float64_t, int8_t, uint32_t>, fun<float64_t, int8_t, uint64_t> }, {fun<float64_t, int16_t, uint8_t>, fun<float64_t, int16_t, uint16_t>, fun<float64_t, int16_t, uint32_t>, fun<float64_t, int16_t, uint64_t> }, {fun<float64_t, int32_t, uint8_t>, fun<float64_t, int32_t, uint16_t>, fun<float64_t, int32_t, uint32_t>, fun<float64_t, int32_t, uint64_t> }, {fun<float64_t, int64_t, uint8_t>, fun<float64_t, int64_t, uint16_t>, fun<float64_t, int64_t, uint32_t>, fun<float64_t, int64_t, uint64_t> }, {fun<float64_t, float32_t, uint8_t>, fun<float64_t, float32_t, uint16_t>, fun<float64_t, float32_t, uint32_t>, fun<float64_t, float32_t, uint64_t> }, {fun<float64_t, float64_t, uint8_t>, fun<float64_t, float64_t, uint16_t>, fun<float64_t, float64_t, uint32_t>, fun<float64_t, float64_t, uint64_t> }, {fun<float64_t, nm::Complex64, uint8_t>, fun<float64_t, nm::Complex64, uint16_t>, fun<float64_t, nm::Complex64, uint32_t>, fun<float64_t, nm::Complex64, uint64_t> }, {fun<float64_t, nm::Complex128, uint8_t>, fun<float64_t, nm::Complex128, uint16_t>, fun<float64_t, nm::Complex128, uint32_t>, fun<float64_t, nm::Complex128, uint64_t> }, {fun<float64_t, nm::Rational32, uint8_t>, fun<float64_t, nm::Rational32, uint16_t>, fun<float64_t, nm::Rational32, uint32_t>, fun<float64_t, nm::Rational32, uint64_t> }, {fun<float64_t, nm::Rational64, uint8_t>, fun<float64_t, nm::Rational64, uint16_t>, fun<float64_t, nm::Rational64, uint32_t>, fun<float64_t, nm::Rational64, uint64_t> }, {fun<float64_t, nm::Rational128, uint8_t>, fun<float64_t, nm::Rational128, uint16_t>, fun<float64_t, nm::Rational128, uint32_t>, fun<float64_t, nm::Rational128, uint64_t> }, {fun<float64_t, nm::RubyObject, uint8_t>, fun<float64_t, nm::RubyObject, uint16_t>, fun<float64_t, nm::RubyObject, uint32_t>, fun<float64_t, nm::RubyObject, uint64_t> }},  \
  {{fun<nm::Complex64, uint8_t, uint8_t>, fun<nm::Complex64, uint8_t, uint16_t>, fun<nm::Complex64, uint8_t, uint32_t>, fun<nm::Complex64, uint8_t, uint64_t> }, {fun<nm::Complex64, int8_t, uint8_t>, fun<nm::Complex64, int8_t, uint16_t>, fun<nm::Complex64, int8_t, uint32_t>, fun<nm::Complex64, int8_t, uint64_t> }, {fun<nm::Complex64, int16_t, uint8_t>, fun<nm::Complex64, int16_t, uint16_t>, fun<nm::Complex64, int16_t, uint32_t>, fun<nm::Complex64, int16_t, uint64_t> }, {fun<nm::Complex64, int32_t, uint8_t>, fun<nm::Complex64, int32_t, uint16_t>, fun<nm::Complex64, int32_t, uint32_t>, fun<nm::Complex64, int32_t, uint64_t> }, {fun<nm::Complex64, int64_t, uint8_t>, fun<nm::Complex64, int64_t, uint16_t>, fun<nm::Complex64, int64_t, uint32_t>, fun<nm::Complex64, int64_t, uint64_t> }, {fun<nm::Complex64, float32_t, uint8_t>, fun<nm::Complex64, float32_t, uint16_t>, fun<nm::Complex64, float32_t, uint32_t>, fun<nm::Complex64, float32_t, uint64_t> }, {fun<nm::Complex64, float64_t, uint8_t>, fun<nm::Complex64, float64_t, uint16_t>, fun<nm::Complex64, float64_t, uint32_t>, fun<nm::Complex64, float64_t, uint64_t> }, {fun<nm::Complex64, nm::Complex64, uint8_t>, fun<nm::Complex64, nm::Complex64, uint16_t>, fun<nm::Complex64, nm::Complex64, uint32_t>, fun<nm::Complex64, nm::Complex64, uint64_t> }, {fun<nm::Complex64, nm::Complex128, uint8_t>, fun<nm::Complex64, nm::Complex128, uint16_t>, fun<nm::Complex64, nm::Complex128, uint32_t>, fun<nm::Complex64, nm::Complex128, uint64_t> }, {fun<nm::Complex64, nm::Rational32, uint8_t>, fun<nm::Complex64, nm::Rational32, uint16_t>, fun<nm::Complex64, nm::Rational32, uint32_t>, fun<nm::Complex64, nm::Rational32, uint64_t> }, {fun<nm::Complex64, nm::Rational64, uint8_t>, fun<nm::Complex64, nm::Rational64, uint16_t>, fun<nm::Complex64, nm::Rational64, uint32_t>, fun<nm::Complex64, nm::Rational64, uint64_t> }, {fun<nm::Complex64, nm::Rational128, uint8_t>, fun<nm::Complex64, nm::Rational128, uint16_t>, fun<nm::Complex64, nm::Rational128, uint32_t>, fun<nm::Complex64, nm::Rational128, uint64_t> }, {fun<nm::Complex64, nm::RubyObject, uint8_t>, fun<nm::Complex64, nm::RubyObject, uint16_t>, fun<nm::Complex64, nm::RubyObject, uint32_t>, fun<nm::Complex64, nm::RubyObject, uint64_t> }},  \
  {{fun<nm::Complex128, uint8_t, uint8_t>, fun<nm::Complex128, uint8_t, uint16_t>, fun<nm::Complex128, uint8_t, uint32_t>, fun<nm::Complex128, uint8_t, uint64_t> }, {fun<nm::Complex128, int8_t, uint8_t>, fun<nm::Complex128, int8_t, uint16_t>, fun<nm::Complex128, int8_t, uint32_t>, fun<nm::Complex128, int8_t, uint64_t> }, {fun<nm::Complex128, int16_t, uint8_t>, fun<nm::Complex128, int16_t, uint16_t>, fun<nm::Complex128, int16_t, uint32_t>, fun<nm::Complex128, int16_t, uint64_t> }, {fun<nm::Complex128, int32_t, uint8_t>, fun<nm::Complex128, int32_t, uint16_t>, fun<nm::Complex128, int32_t, uint32_t>, fun<nm::Complex128, int32_t, uint64_t> }, {fun<nm::Complex128, int64_t, uint8_t>, fun<nm::Complex128, int64_t, uint16_t>, fun<nm::Complex128, int64_t, uint32_t>, fun<nm::Complex128, int64_t, uint64_t> }, {fun<nm::Complex128, float32_t, uint8_t>, fun<nm::Complex128, float32_t, uint16_t>, fun<nm::Complex128, float32_t, uint32_t>, fun<nm::Complex128, float32_t, uint64_t> }, {fun<nm::Complex128, float64_t, uint8_t>, fun<nm::Complex128, float64_t, uint16_t>, fun<nm::Complex128, float64_t, uint32_t>, fun<nm::Complex128, float64_t, uint64_t> }, {fun<nm::Complex128, nm::Complex64, uint8_t>, fun<nm::Complex128, nm::Complex64, uint16_t>, fun<nm::Complex128, nm::Complex64, uint32_t>, fun<nm::Complex128, nm::Complex64, uint64_t> }, {fun<nm::Complex128, nm::Complex128, uint8_t>, fun<nm::Complex128, nm::Complex128, uint16_t>, fun<nm::Complex128, nm::Complex128, uint32_t>, fun<nm::Complex128, nm::Complex128, uint64_t> }, {fun<nm::Complex128, nm::Rational32, uint8_t>, fun<nm::Complex128, nm::Rational32, uint16_t>, fun<nm::Complex128, nm::Rational32, uint32_t>, fun<nm::Complex128, nm::Rational32, uint64_t> }, {fun<nm::Complex128, nm::Rational64, uint8_t>, fun<nm::Complex128, nm::Rational64, uint16_t>, fun<nm::Complex128, nm::Rational64, uint32_t>, fun<nm::Complex128, nm::Rational64, uint64_t> }, {fun<nm::Complex128, nm::Rational128, uint8_t>, fun<nm::Complex128, nm::Rational128, uint16_t>, fun<nm::Complex128, nm::Rational128, uint32_t>, fun<nm::Complex128, nm::Rational128, uint64_t> }, {fun<nm::Complex128, nm::RubyObject, uint8_t>, fun<nm::Complex128, nm::RubyObject, uint16_t>, fun<nm::Complex128, nm::RubyObject, uint32_t>, fun<nm::Complex128, nm::RubyObject, uint64_t> }},  \
  {{fun<nm::Rational32, uint8_t, uint8_t>, fun<nm::Rational32, uint8_t, uint16_t>, fun<nm::Rational32, uint8_t, uint32_t>, fun<nm::Rational32, uint8_t, uint64_t> }, {fun<nm::Rational32, int8_t, uint8_t>, fun<nm::Rational32, int8_t, uint16_t>, fun<nm::Rational32, int8_t, uint32_t>, fun<nm::Rational32, int8_t, uint64_t> }, {fun<nm::Rational32, int16_t, uint8_t>, fun<nm::Rational32, int16_t, uint16_t>, fun<nm::Rational32, int16_t, uint32_t>, fun<nm::Rational32, int16_t, uint64_t> }, {fun<nm::Rational32, int32_t, uint8_t>, fun<nm::Rational32, int32_t, uint16_t>, fun<nm::Rational32, int32_t, uint32_t>, fun<nm::Rational32, int32_t, uint64_t> }, {fun<nm::Rational32, int64_t, uint8_t>, fun<nm::Rational32, int64_t, uint16_t>, fun<nm::Rational32, int64_t, uint32_t>, fun<nm::Rational32, int64_t, uint64_t> }, {NULL, NULL, NULL, NULL }, {NULL, NULL, NULL, NULL }, {NULL, NULL, NULL, NULL }, {NULL, NULL, NULL, NULL }, {fun<nm::Rational32, nm::Rational32, uint8_t>, fun<nm::Rational32, nm::Rational32, uint16_t>, fun<nm::Rational32, nm::Rational32, uint32_t>, fun<nm::Rational32, nm::Rational32, uint64_t> }, {fun<nm::Rational32, nm::Rational64, uint8_t>, fun<nm::Rational32, nm::Rational64, uint16_t>, fun<nm::Rational32, nm::Rational64, uint32_t>, fun<nm::Rational32, nm::Rational64, uint64_t> }, {fun<nm::Rational32, nm::Rational128, uint8_t>, fun<nm::Rational32, nm::Rational128, uint16_t>, fun<nm::Rational32, nm::Rational128, uint32_t>, fun<nm::Rational32, nm::Rational128, uint64_t> }, {fun<nm::Rational32, nm::RubyObject, uint8_t>, fun<nm::Rational32, nm::RubyObject, uint16_t>, fun<nm::Rational32, nm::RubyObject, uint32_t>, fun<nm::Rational32, nm::RubyObject, uint64_t>}},  \
  {{fun<nm::Rational64, uint8_t, uint8_t>, fun<nm::Rational64, uint8_t, uint16_t>, fun<nm::Rational64, uint8_t, uint32_t>, fun<nm::Rational64, uint8_t, uint64_t> }, {fun<nm::Rational64, int8_t, uint8_t>, fun<nm::Rational64, int8_t, uint16_t>, fun<nm::Rational64, int8_t, uint32_t>, fun<nm::Rational64, int8_t, uint64_t> }, {fun<nm::Rational64, int16_t, uint8_t>, fun<nm::Rational64, int16_t, uint16_t>, fun<nm::Rational64, int16_t, uint32_t>, fun<nm::Rational64, int16_t, uint64_t> }, {fun<nm::Rational64, int32_t, uint8_t>, fun<nm::Rational64, int32_t, uint16_t>, fun<nm::Rational64, int32_t, uint32_t>, fun<nm::Rational64, int32_t, uint64_t> }, {fun<nm::Rational64, int64_t, uint8_t>, fun<nm::Rational64, int64_t, uint16_t>, fun<nm::Rational64, int64_t, uint32_t>, fun<nm::Rational64, int64_t, uint64_t> }, {NULL, NULL, NULL, NULL }, {NULL, NULL, NULL, NULL }, {NULL, NULL, NULL, NULL }, {NULL, NULL, NULL, NULL }, {fun<nm::Rational64, nm::Rational32, uint8_t>, fun<nm::Rational64, nm::Rational32, uint16_t>, fun<nm::Rational64, nm::Rational32, uint32_t>, fun<nm::Rational64, nm::Rational32, uint64_t> }, {fun<nm::Rational64, nm::Rational64, uint8_t>, fun<nm::Rational64, nm::Rational64, uint16_t>, fun<nm::Rational64, nm::Rational64, uint32_t>, fun<nm::Rational64, nm::Rational64, uint64_t> }, {fun<nm::Rational64, nm::Rational128, uint8_t>, fun<nm::Rational64, nm::Rational128, uint16_t>, fun<nm::Rational64, nm::Rational128, uint32_t>, fun<nm::Rational64, nm::Rational128, uint64_t> }, {fun<nm::Rational64, nm::RubyObject, uint8_t>, fun<nm::Rational64, nm::RubyObject, uint16_t>, fun<nm::Rational64, nm::RubyObject, uint32_t>, fun<nm::Rational64, nm::RubyObject, uint64_t>}},  \
  {{fun<nm::Rational128, uint8_t, uint8_t>, fun<nm::Rational128, uint8_t, uint16_t>, fun<nm::Rational128, uint8_t, uint32_t>, fun<nm::Rational128, uint8_t, uint64_t> }, {fun<nm::Rational128, int8_t, uint8_t>, fun<nm::Rational128, int8_t, uint16_t>, fun<nm::Rational128, int8_t, uint32_t>, fun<nm::Rational128, int8_t, uint64_t> }, {fun<nm::Rational128, int16_t, uint8_t>, fun<nm::Rational128, int16_t, uint16_t>, fun<nm::Rational128, int16_t, uint32_t>, fun<nm::Rational128, int16_t, uint64_t> }, {fun<nm::Rational128, int32_t, uint8_t>, fun<nm::Rational128, int32_t, uint16_t>, fun<nm::Rational128, int32_t, uint32_t>, fun<nm::Rational128, int32_t, uint64_t> }, {fun<nm::Rational128, int64_t, uint8_t>, fun<nm::Rational128, int64_t, uint16_t>, fun<nm::Rational128, int64_t, uint32_t>, fun<nm::Rational128, int64_t, uint64_t> }, {fun<nm::Rational128, float32_t, uint8_t>, fun<nm::Rational128, float32_t, uint16_t>, fun<nm::Rational128, float32_t, uint32_t>, fun<nm::Rational128, float32_t, uint64_t> }, {NULL, NULL, NULL, NULL }, {NULL, NULL, NULL, NULL }, {NULL, NULL, NULL, NULL }, {fun<nm::Rational128, nm::Rational32, uint8_t>, fun<nm::Rational128, nm::Rational32, uint16_t>, fun<nm::Rational128, nm::Rational32, uint32_t>, fun<nm::Rational128, nm::Rational32, uint64_t> }, {fun<nm::Rational128, nm::Rational64, uint8_t>, fun<nm::Rational128, nm::Rational64, uint16_t>, fun<nm::Rational128, nm::Rational64, uint32_t>, fun<nm::Rational128, nm::Rational64, uint64_t> }, {fun<nm::Rational128, nm::Rational128, uint8_t>, fun<nm::Rational128, nm::Rational128, uint16_t>, fun<nm::Rational128, nm::Rational128, uint32_t>, fun<nm::Rational128, nm::Rational128, uint64_t> }, {fun<nm::Rational128, nm::RubyObject, uint8_t>, fun<nm::Rational128, nm::RubyObject, uint16_t>, fun<nm::Rational128, nm::RubyObject, uint32_t>, fun<nm::Rational128, nm::RubyObject, uint64_t>}},  \
  {{fun<nm::RubyObject, uint8_t, uint8_t>, fun<nm::RubyObject, uint8_t, uint16_t>, fun<nm::RubyObject, uint8_t, uint32_t>, fun<nm::RubyObject, uint8_t, uint64_t> }, {fun<nm::RubyObject, int8_t, uint8_t>, fun<nm::RubyObject, int8_t, uint16_t>, fun<nm::RubyObject, int8_t, uint32_t>, fun<nm::RubyObject, int8_t, uint64_t> }, {fun<nm::RubyObject, int16_t, uint8_t>, fun<nm::RubyObject, int16_t, uint16_t>, fun<nm::RubyObject, int16_t, uint32_t>, fun<nm::RubyObject, int16_t, uint64_t> }, {fun<nm::RubyObject, int32_t, uint8_t>, fun<nm::RubyObject, int32_t, uint16_t>, fun<nm::RubyObject, int32_t, uint32_t>, fun<nm::RubyObject, int32_t, uint64_t> }, {fun<nm::RubyObject, int64_t, uint8_t>, fun<nm::RubyObject, int64_t, uint16_t>, fun<nm::RubyObject, int64_t, uint32_t>, fun<nm::RubyObject, int64_t, uint64_t> }, {fun<nm::RubyObject, float32_t, uint8_t>, fun<nm::RubyObject, float32_t, uint16_t>, fun<nm::RubyObject, float32_t, uint32_t>, fun<nm::RubyObject, float32_t, uint64_t> }, {fun<nm::RubyObject, float64_t, uint8_t>, fun<nm::RubyObject, float64_t, uint16_t>, fun<nm::RubyObject, float64_t, uint32_t>, fun<nm::RubyObject, float64_t, uint64_t> }, {fun<nm::RubyObject, nm::Complex64, uint8_t>, fun<nm::RubyObject, nm::Complex64, uint16_t>, fun<nm::RubyObject, nm::Complex64, uint32_t>, fun<nm::RubyObject, nm::Complex64, uint64_t> }, {fun<nm::RubyObject, nm::Complex128, uint8_t>, fun<nm::RubyObject, nm::Complex128, uint16_t>, fun<nm::RubyObject, nm::Complex128, uint32_t>, fun<nm::RubyObject, nm::Complex128, uint64_t> }, {fun<nm::RubyObject, nm::Rational32, uint8_t>, fun<nm::RubyObject, nm::Rational32, uint16_t>, fun<nm::RubyObject, nm::Rational32, uint32_t>, fun<nm::RubyObject, nm::Rational32, uint64_t> }, {fun<nm::RubyObject, nm::Rational64, uint8_t>, fun<nm::RubyObject, nm::Rational64, uint16_t>, fun<nm::RubyObject, nm::Rational64, uint32_t>, fun<nm::RubyObject, nm::Rational64, uint64_t> }, {fun<nm::RubyObject, nm::Rational128, uint8_t>, fun<nm::RubyObject, nm::Rational128, uint16_t>, fun<nm::RubyObject, nm::Rational128, uint32_t>, fun<nm::RubyObject, nm::Rational128, uint64_t> }, {fun<nm::RubyObject, nm::RubyObject, uint8_t>, fun<nm::RubyObject, nm::RubyObject, uint16_t>, fun<nm::RubyObject, nm::RubyObject, uint32_t>, fun<nm::RubyObject, nm::RubyObject, uint64_t>}} \
};

/*
 * Defines a static array that holds function pointers to left dtype and itype
 * templated versions of the specified function.
 */
#define LI_DTYPE_TEMPLATE_TABLE(fun, ret, ...) NAMED_LI_DTYPE_TEMPLATE_TABLE(ttable, fun, ret, __VA_ARGS__)

#define NAMED_LI_DTYPE_TEMPLATE_TABLE(name, fun, ret, ...) 																																			\
	static ret (*(name)[nm::NUM_DTYPES][nm::NUM_ITYPES])(__VA_ARGS__) = {																																	\
		{ fun<uint8_t,uint8_t>,fun<uint8_t,uint16_t>,fun<uint8_t,uint32_t>,fun<uint8_t,uint64_t> },																	\
		{ fun<int8_t,uint8_t>,fun<int8_t,uint16_t>,fun<int8_t,uint32_t>,fun<int8_t,uint64_t> },																			\
		{ fun<int16_t,uint8_t>,fun<int16_t,uint16_t>,fun<int16_t,uint32_t>,fun<int16_t,uint64_t> },																	\
		{ fun<int32_t,uint8_t>,fun<int32_t,uint16_t>,fun<int32_t,uint32_t>,fun<int32_t,uint64_t> },																	\
		{ fun<int64_t,uint8_t>,fun<int64_t,uint16_t>,fun<int64_t,uint32_t>,fun<int64_t,uint64_t> },																	\
		{ fun<float32_t,uint8_t>,fun<float32_t,uint16_t>,fun<float32_t,uint32_t>,fun<float32_t,uint64_t> },													\
		{ fun<float64_t,uint8_t>,fun<float64_t,uint16_t>,fun<float64_t,uint32_t>,fun<float64_t,uint64_t> },													\
		{ fun<nm::Complex64,uint8_t>,fun<nm::Complex64,uint16_t>,fun<nm::Complex64,uint32_t>,fun<nm::Complex64,uint64_t> },					\
		{ fun<nm::Complex128,uint8_t>,fun<nm::Complex128,uint16_t>,fun<nm::Complex128,uint32_t>,fun<nm::Complex128,uint64_t> },			\
		{ fun<nm::Rational32,uint8_t>,fun<nm::Rational32,uint16_t>,fun<nm::Rational32,uint32_t>,fun<nm::Rational32,uint64_t> },			\
		{ fun<nm::Rational64,uint8_t>,fun<nm::Rational64,uint16_t>,fun<nm::Rational64,uint32_t>,fun<nm::Rational64,uint64_t> },			\
		{ fun<nm::Rational128,uint8_t>,fun<nm::Rational128,uint16_t>,fun<nm::Rational128,uint32_t>,fun<nm::Rational128,uint64_t> },	\
		{ fun<nm::RubyObject,uint8_t>,fun<nm::RubyObject,uint16_t>,fun<nm::RubyObject,uint32_t>,fun<nm::RubyObject,uint64_t>}				\
	};

} // end of extern "C" block


extern "C" {


/*
 * Data
 */

// regular data types
extern const char* const	DTYPE_NAMES[nm::NUM_DTYPES];
extern const size_t 			DTYPE_SIZES[nm::NUM_DTYPES];

// index data types
extern const char* const  ITYPE_NAMES[nm::NUM_ITYPES];
extern const size_t 			ITYPE_SIZES[nm::NUM_ITYPES];

extern const dtype_t      Upcast[nm::NUM_DTYPES][nm::NUM_DTYPES];

/*
 * Functions
 */


void*	    			rubyobj_to_cval(VALUE val, dtype_t dtype);
void  		  		rubyval_to_cval(VALUE val, dtype_t dtype, void* loc);
nm::RubyObject	rubyobj_from_cval(void* val, dtype_t dtype);
nm::RubyObject  rubyobj_from_cval_by_itype(void* val, itype_t itype);

} // end of extern "C" block

#endif // DATA_H
