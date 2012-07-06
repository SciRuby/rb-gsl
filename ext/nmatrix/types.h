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
// == types.h
//
// Resources for all storage types.

#ifndef NMATRIX_TYPES_H
#define NMATRIX_TYPES_H

// This include is needed here so we can figure out what else to include.
#include "nmatrix_config.h"

/*
 * Standard Includes
 */

#ifdef HAVE_STDINT_H
	#include <stdint.h>
#endif

/*
 * Project Includes
 */

/*
 * Macros
 */

/*
 * Types
 */

///////////////////////////////////////////////////////////////////////////
// Please modify these types if your system has any different type sizes //
///////////////////////////////////////////////////////////////////////////

// NM_BYTE : unsigned 8-bit integer
#ifndef HAVE_U_INT8_T
	#ifdef HAVE_UINT8_T
		typedef uint8_t						u_int8_t;
	#else
		typedef unsigned char			u_int8_t;
	#endif
#endif

#ifndef HAVE_INT8_T
	typedef char								int8_t;
#endif

#ifndef HAVE_INT16_T
	#if SIZEOF_SHORT == 2
		typedef short							int16_t;
	#else
		---->> Please define int16_t manually because sizeof(short) != 2. <<----
	# endif
#endif

#ifndef HAVE_INT32_T
	#if SIZEOF_LONG == 4
		typedef long							int32_t;
	#else
		#if SIZEOF_INT == 4
			typedef int							int32_t;
		#else
			---->> Please define int32_t manually because sizeof(long) != 4. <<----
		#endif
	#endif
#endif

// unsigned 32-bit integer
#ifndef HAVE_U_INT32_T
	#ifdef HAVE_UINT32_T
		typedef uint32_t					u_int32_t;
	#else
		#if SIZEOF_LONG == 4
			typedef unsigned long		u_int32_t;
		#else
			#if SIZEOF_INT == 4
				typedef unsigned int	u_int32_t;
			#else
				---->> Please define u_int32_t manually because sizeof(long) != 4. <<----
			#endif
		#endif
	#endif
#endif

#ifndef HAVE_INT64_T
	#if SIZEOF_QUAD == 8
		typedef quad							int64_t;
	#else
		#if SIZEOF_LONG == 8
			typedef long						int64_t;
		#else
			---->> Please define int64_t manually because sizeof(quad) != 8. <<----
		#endif
	#endif
#endif

// unsigned 64-bit integer
#ifndef HAVE_U_INT64_T
	#ifdef HAVE_UINT64_T
		typedef uint64_t					u_int64_t;
	#else
		#if SIZEOF_QUAD == 8
			typedef unsigned quad		u_int64_t;
		#else
			#if SIZEOF_LONG == 8
				typedef unsigned long	u_int64_t;
			#else
				---->> Please define u_int64_t manually because sizeof(quad) != 8. <<----
			#endif
		#endif
	#endif
#endif

// If you modify this, make sure to modify the definition of y_size_t and Y_SIZE_T!
#ifndef HAVE_SIZE_T
	typedef u_int64_t						size_t;
	#define NM_SIZE_T NM_INT64
#else
	#if SIZEOF_SIZE_T == 8
		#define NM_SIZE_T NM_INT64
	#else
		#if SIZEOF_SIZE_T == 4
			#define NM_SIZE_T NM_INT32
		#else
			---->> Please define size_t and y_size_t manually because sizeof(size_t) is neither 8 nor 4. <<----
		#endif
	#endif
#endif

typedef float		float32;
typedef double	float64;

/*
 * For when we need to return array indices. This must never be larger than
 * size_t.
 */
typedef uint32_t    y_size_t;
#define Y_SIZE_T    NM_INT32

#ifdef HAVE_STDBOOL_H
	#include <stdbool.h>
#else
	typedef char    bool;
	#define true    1;
	#define false   0;
#endif

/*
 * Data
 */

/*
 * Functions
 */

#endif
