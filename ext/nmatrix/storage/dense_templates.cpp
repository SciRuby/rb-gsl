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

template <typename DType, typename NewDType>
DENSE_STORAGE* dense_storage_cast_copy_template(const DENSE_STORAGE* rhs, dtype_t new_dtype);

template <typename LDType, typename RDType>
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
	 * Do these two dense matrices have the same contents?
	 */
	bool dense_storage_eqeq(const DENSE_STORAGE* left, const DENSE_STORAGE* right) {
		LR_DTYPE_TEMPLATE_TABLE(dense_storage_eqeq_template, bool, const DENSE_STORAGE*, const DENSE_STORAGE*);
		
		return ttable[left->dtype][right->dtype](left, right);
	}
	
	/*
	 * Test to see if the matrix is Hermitian.  If the matrix does not have a
	 * dtype of Complex64 or Complex128 this is the same as testing for symmetry.
	 */
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
		
		return ttable[mat->dtype](mat, lda);
	}
	
	/////////////////////////
	// Casting and Copying //
	/////////////////////////
	
	/*
	 * Documentation goes here.
	 */
	DENSE_STORAGE* dense_storage_cast_copy(const DENSE_STORAGE* rhs, dtype_t new_dtype) {
		LR_DTYPE_TEMPLATE_TABLE(dense_storage_cast_copy_template, DENSE_STORAGE*, const DENSE_STORAGE*, dtype_t);
		
		return ttable[new_dtype][rhs->dtype](rhs, new_dtype);
	}
}

/////////////////////////
// Templated Functions //
/////////////////////////

template <typename DType, typename NewDType>
DENSE_STORAGE* dense_storage_cast_copy_template(const DENSE_STORAGE* rhs, dtype_t new_dtype) {
  size_t  count = storage_count_max_elements(rhs->rank, rhs->shape), p;
  size_t* shape = ALLOC_N(size_t, rhs->rank);
  
  DType*		rhs_els = (DType*)rhs->elements;
  NewDType* lhs_els;
  
  DENSE_STORAGE* lhs;
  
  if (!shape) {
  	return NULL;
  }
	
  // Copy shape array.
  for (p = rhs->rank; p-- > 0;) {
  	shape[p] = rhs->shape[p];
  }
	
  lhs			= dense_storage_create(new_dtype, shape, rhs->rank, NULL, 0);
  lhs_els	= (NewDType*)lhs->elements;

	// Ensure that allocation worked before copying.
  if (lhs && count) {
    if (lhs->dtype == rhs->dtype) {
      memcpy(lhs->elements, rhs->elements, DTYPE_SIZES[rhs->dtype] * count);
      
    } else {
    	while (count-- > 0) {
    		lhs_els[count] = rhs_els[count];
      }
    }
  }
	
  return lhs;
}

template <typename LDType, typename RDType>
bool dense_storage_eqeq_template(const DENSE_STORAGE* left, const DENSE_STORAGE* right) {
	int index;
	
	LDType* left_els	= (LDType*)left->elements;
	RDType* right_els	= (RDType*)right->elements;
	
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
