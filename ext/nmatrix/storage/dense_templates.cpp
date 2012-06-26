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
// == dense_templates.c
//
// Templates for dense n-dimensional matrix storage.

/*
 * Standard Includes
 */

#include <ruby.h>

/*
 * Project Includes
 */

#include "dense.h"
#include "dense_templates.h"

/*
 * Macros
 */

#define TEMPLATE_SWITCH(fun, dtype, ...) switch(dtype) {	\
	case CHAR:																							\
		return fun<char>(__VA_ARGS__);												\
																													\
	case UCHAR:																							\
		return fun<unsigned char>(__VA_ARGS__);								\
																													\
	case SHORT:																							\
		return fun<short>(__VA_ARGS__);												\
																													\
	case USHORT:																						\
		return fun<unsigned short>(__VA_ARGS__);							\
																													\
	case INT:																								\
		return fun<int>(__VA_ARGS__);													\
																													\
	case UINT:																							\
		return fun<unsigned int>(__VA_ARGS__);								\
																													\
	case FLOAT:																							\
		return fun<float>(__VA_ARGS__);												\
																													\
	case DOUBLE:																						\
		return fun<double>(__VA_ARGS__);											\
	}

/*
 * Global Variables
 */

extern bool				(*ElemEqEq[NM_TYPES][2])(const void*, const void*, const int, const int);
extern const int	nm_sizeof[NM_TYPES];

/*
 * Forward Declarations
 */

/*
 * Functions
 */

/////////////////////////
// C Wrapper Functions //
/////////////////////////

extern "C" {
	/*
	 * Do these two dense matrices of the same dtype have exactly the same
	 * contents?
	 */
	bool dense_storage_eqeq(const DENSE_STORAGE* left, const DENSE_STORAGE* right) {
		TEMPLATE_SWITCH(dense_storage_eqeq_template, left->dtype, left, right)
	}

	/*
	 * Is this dense matrix symmetric about the diagonal?
	 *
	 * FIXME: Add templating.
	 */
	bool dense_storage_is_symmetric(const DENSE_STORAGE* mat, int lda, bool hermitian) {
		TEMPLATE_SWITCH(dense_storage_is_semetric_template, mat->dtype, mat, lda, hermitian)
	}
}

/////////////////////////
// Templated Functions //
/////////////////////////

template <typename DType>
bool dense_storage_eqeq_template(const DENSE_STORAGE* left, const DENSE_STORAGE* right) {
	int index;
	
	DType* left_els		= (DType*)left->elements;
	DType* right_els	= (DType*)right->elements;
	
	for (index = count_dense_storage_elements(left); index-- > 0;) {
		if (left_els[index] != right_els[index]) {
			return false;
		}
	}
	
	return true;
}

template <typename DType>
bool dense_storage_is_symmetric_template(const DENSE_STORAGE* mat, int lda, bool hermitian) {
	unsigned int i, j;
	const DType* els = (DType*) mat->elements;
	const DType* a, * b;
	
	for (i = mat->shape[0]; i-- > 0;) {
		for (j = i + 1; j < mat->shape[1]; ++j) {
			a = mat->elements + (i*lda+j);
	  	b = mat->elements + (j*lda+i);
	  	
	  	if (*a != *b) {
	      return false;
	    }
		}
	}
	
	return true;
}
