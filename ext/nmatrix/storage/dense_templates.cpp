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
// == dense_templates.cpp
//
// Templates for dense n-dimensional matrix storage.

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

#include "data/data.h"

#include "dense.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

/*
 * Forward Declarations
 */

template <typename DType>
bool dense_storage_eqeq_template(const DENSE_STORAGE* left, const DENSE_STORAGE* right);

template <typename DType>
bool dense_storage_is_hermitian_template(const DENSE_STORAGE* mat, int lda);

template <typename DType>
bool dense_storage_is_symmetric_template(const DENSE_STORAGE* mat, int lda);

/*
 * Functions
 */

/////////////////////////
// C Wrapper Functions //
/////////////////////////

extern "C" {
	
	///////////
	// Tests //
	///////////	
	
	/*
	 * Do these two dense matrices of the same dtype have exactly the same
	 * contents?
	 */
	bool dense_storage_eqeq(const DENSE_STORAGE* left, const DENSE_STORAGE* right) {
		DTYPE_TEMPLATE_TABLE(dense_storage_eqeq_template, bool, const DENSE_STORAGE*, const DENSE_STORAGE*);
		
		return ttable[left->dtype](left, right);
	}
	
	bool dense_storage_is_hermitian(const DENSE_STORAGE* mat, int lda) {
		if (mat->dtype == COMPLEX64) {
			return dense_storage_is_hermitian_template<Complex64>(mat, lda);
			
		} else if (mat->dtype == COMPLEX128) {
			return dense_storage_is_hermitian_template<Complex128>(mat, lda);
			
		} else {
			return dense_storage_is_symmetric(mat, lda);
		}
	}
	
	/*
	 * Is this dense matrix symmetric about the diagonal?
	 */
	bool dense_storage_is_symmetric(const DENSE_STORAGE* mat, int lda) {
		DTYPE_TEMPLATE_TABLE(dense_storage_is_symmetric_template, bool, const DENSE_STORAGE*, int);
		
		return ttable[mat->dtype];
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
	
	for (index = storage_count_max_elements(left->rank, left->shape); index-- > 0;) {
		if (left_els[index] != right_els[index]) {
			return false;
		}
	}
	
	return true;
}

template <typename DType>
bool dense_storage_is_hermitian_template(const DENSE_STORAGE* mat, int lda) {
	unsigned int i, j;
	register DType complex_conj;
	
	const DType* els = (DType*) mat->elements;
	
	for (i = mat->shape[0]; i-- > 0;) {
		for (j = i + 1; j < mat->shape[1]; ++j) {
			complex_conj		= els[j*lda + 1];
			complex_conj.i	= -complex_conj.i;
			
			if (els[i*lda+j] != complex_conj) {
	      return false;
	    }
		}
	}
	
	return true;
}

template <typename DType>
bool dense_storage_is_symmetric_template(const DENSE_STORAGE* mat, int lda) {
	unsigned int i, j;
	const DType* els = (DType*) mat->elements;
	
	for (i = mat->shape[0]; i-- > 0;) {
		for (j = i + 1; j < mat->shape[1]; ++j) {
			if (els[i*lda+j] != els[j*lda+i]) {
	      return false;
	    }
		}
	}
	
	return true;
}
