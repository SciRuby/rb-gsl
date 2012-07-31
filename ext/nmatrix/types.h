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
// Definition of simple types used throughout NMatrix.

#ifndef NMATRIX_TYPES_H
#define NMATRIX_TYPES_H

/*
 * Standard Includes
 */

#include <stdint.h>

/*
 * Project Includes
 */

/*
 * Macros
 */

#define EPSILON 1E-10
#define FP_IS_ZERO(n) (-EPSILON < n && n < EPSILON)
#define FP_EQUAL(a, b) FP_IS_ZERO((a - b))

/*
 * Types
 */

typedef float		float32_t;
typedef double	float64_t;


#ifndef HAVE_SIZE_T /// If you modify this, make sure to modify the definition of y_size_t and Y_SIZE_T!
typedef u_int64_t    size_t;
# define SIZE_T   INT64
#else
# if SIZEOF_SIZE_T == 8
#  define SIZE_T  INT64
# else
#  if SIZEOF_SIZE_T == 4
#   define SIZE_T INT32
#  else
---->> Please define size_t and y_size_t manually because sizeof(size_t) is neither 8 nor 4. <<----
#  endif
# endif
#endif


typedef uint32_t y_size_t;

/*
 * Data
 */

/*
 * Functions
 */

#endif
