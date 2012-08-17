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
// == dense.h
//
// Dense n-dimensional matrix storage.

#ifndef DENSE_H
#define DENSE_H

/*
 * Standard Includes
 */

#include <stdlib.h>

/*
 * Project Includes
 */

#include "types.h"
//#include "util/math.h"

#include "data/data.h"

#include "common.h"

/*
 * Macros
 */

/*
 * Types
 */

struct DENSE_STORAGE : STORAGE {
	size_t*	stride;
	int			count;
	void*		src;
	void*		elements;
};

/*
 * Data
 */
 

/*
 * Functions
 */

///////////////
// Lifecycle //
///////////////

DENSE_STORAGE*	dense_storage_create(dtype_t dtype, size_t* shape, size_t rank, void* elements, size_t elements_length);
void						dense_storage_delete(STORAGE* s);
void						dense_storage_delete_ref(STORAGE* s);
void						dense_storage_mark(void*);

///////////////
// Accessors //
///////////////

void*	dense_storage_get(STORAGE* s, SLICE* slice);
void*	dense_storage_ref(STORAGE* s, SLICE* slice);
void	dense_storage_set(STORAGE* s, SLICE* slice, void* val);

///////////
// Tests //
///////////

bool dense_storage_eqeq(const STORAGE* left, const STORAGE* right);
bool dense_storage_is_symmetric(const DENSE_STORAGE* mat, int lda);
bool dense_storage_is_hermitian(const DENSE_STORAGE* mat, int lda);

//////////
// Math //
//////////

STORAGE* dense_storage_ew_add(const STORAGE* left, const STORAGE* right);
STORAGE* dense_storage_ew_subtract(const STORAGE* left, const STORAGE* right);
STORAGE* dense_storage_ew_multiply(const STORAGE* left, const STORAGE* right);
STORAGE* dense_storage_ew_divide(const STORAGE* left, const STORAGE* right);
//STORAGE* dense_storage_ew_mod(const STORAGE* left, const STORAGE* right);

STORAGE* dense_storage_matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector);

/////////////
// Utility //
/////////////

size_t dense_storage_pos(const DENSE_STORAGE* s, const size_t* coords);

/////////////////////////
// Copying and Casting //
/////////////////////////

DENSE_STORAGE*  dense_storage_copy(const DENSE_STORAGE* rhs);
STORAGE*        dense_storage_copy_transposed(const STORAGE* rhs_base);
STORAGE*        dense_storage_cast_copy(const STORAGE* rhs, dtype_t new_dtype);

#endif // DENSE_H
