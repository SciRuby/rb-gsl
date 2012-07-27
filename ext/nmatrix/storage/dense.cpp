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
// == dense.c
//
// Dense n-dimensional matrix storage.

/*
 * Standard Includes
 */

#include <ruby.h>

/*
 * Project Includes
 */

#include "types.h"

#include "data/data.h"

#include "common.h"
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

///////////////
// Lifecycle //
///////////////

/*
 * Note that elements and elements_length are for initial value(s) passed in.
 * If they are the correct length, they will be used directly. If not, they
 * will be concatenated over and over again into a new elements array. If
 * elements is NULL, the new elements array will not be initialized.
 */
DENSE_STORAGE* dense_storage_create(dtype_t dtype, size_t* shape, size_t rank, void* elements, size_t elements_length) {
  DENSE_STORAGE* s;
  size_t count, i, copy_length = elements_length;

  s = ALLOC( DENSE_STORAGE );

  s->rank       = rank;
  s->shape      = shape;
  s->dtype      = dtype;
  s->offset     = (size_t*) calloc(sizeof(size_t), rank);
  s->stride     = dense_storage_stride(shape, rank);
  s->count      = 1;
  s->src        = s;
	
	count         = storage_count_max_elements(s->rank, s->shape);

  if (elements_length == count) {
  	s->elements = elements;
  	
  } else {
    s->elements = ALLOC_N(char, DTYPE_SIZES[dtype]*count);

    if (elements_length > 0) {
      // Repeat elements over and over again until the end of the matrix.
      for (i = 0; i < count; i += elements_length) {
        
        if (i + elements_length > count) {
        	copy_length = count - i;
        }
        
        memcpy((char*)(s->elements)+i*DTYPE_SIZES[dtype], (char*)(elements)+(i % elements_length)*DTYPE_SIZES[dtype], copy_length*DTYPE_SIZES[dtype]);
      }

      // Get rid of the init_val.
      free(elements);
    }
  }

  return s;
}

/*
 * Documentation goes here.
 */
void dense_storage_delete(DENSE_STORAGE* s) {
  // Sometimes Ruby passes in NULL storage for some reason (probably on copy construction failure).
  if (s) {
    if(s->count-- == 1) {
      free(s->shape);
      free(s->offset);
      free(s->stride);
      free(s->elements);
      free(s);
    }
  }
}

/*
 * Documentation goes here.
 */
void dense_storage_delete_ref(DENSE_STORAGE* s) {
  // Sometimes Ruby passes in NULL storage for some reason (probably on copy construction failure).
  if (s) {
    dense_storage_delete(s->src);
    free(s->shape);
    free(s->offset);
    free(s);
  }
}

/*
 * Mark values in a dense matrix for garbage collection. This may not be necessary -- further testing required.
 */
void dense_storage_mark(DENSE_STORAGE* storage) {
  size_t i;

  if (storage && storage->dtype == RUBYOBJ) {
  	for (i = storage_count_max_elements(storage->rank, storage->shape); i-- > 0;) {
      rb_gc_mark(*((VALUE*)((char*)(storage->elements) + i*DTYPE_SIZES[RUBYOBJ])));
    }
  }
}

///////////////
// Accessors //
///////////////

/*
 * Get a slice or one element, using copying.
 */
void* dense_storage_get(DENSE_STORAGE* s, SLICE* slice) {
  DENSE_STORAGE *ns;
  size_t count;

  if (slice->is_one_el)
    return (char*)(s->elements) + dense_storage_pos(s, slice->coords) * DTYPE_SIZES[s->dtype];
  else { // Make references
    ns = ALLOC( DENSE_STORAGE );

    NM_CHECK_ALLOC(ns);

    ns->rank       = s->rank;
    ns->shape      = slice->lengths;
    ns->dtype      = s->dtype;

    ns->offset     = calloc(sizeof(size_t),ns->rank);
    NM_CHECK_ALLOC(ns->offset);

    ns->strides    = dense_storage_stride(ns->shape, ns->rank);
    ns->count      = 1;
    ns->src        = ns;

    count          = storage_count_max_elements( s->rank, s->shape );

    ns->elements   = ALLOC_N(char, DTYPE_SIZES[ns->dtype] * count);
    NM_CHECK_ALLOC(ns->elements);

    dense_storage_slice_copy(s, ns, slice->lengths, dense_storage_pos(s, slice->coords), 0, 0);
    return ns;
  }
}


/*
 * Documentation goes here.
 */
void* dense_storage_ref(DENSE_STORAGE* s, SLICE* slice) {
  DENSE_STORAGE *ns;

  if (slice->is_one_el)
    return (char*)(s->elements) + dense_storage_pos(s, slice->coords) * DTYPE_SIZES[s->dtype];
    
  else {
    ns = ALLOC( DENSE_STORAGE );
    NM_CHECK_ALLOC(ns);

    ns->rank       = s->rank;
    ns->dtype      = s->dtype;

    ns->offset     = calloc(sizeof(*ns->offset), ns->rank);
    NM_CHECK_ALLOC(ns->offset);

    ns->shape      = calloc(sizeof(*ns->offset), ns->rank);
    NM_CHECK_ALLOC(ns->shape);

    for (index = 0; index < ns->rank; ++index) {
      ns->offset[i] = slice->coords[i] + s->offset[i];
      ns->shape[i]  = slice->lengths[i];
    }

    ns->stride     = s->stride;
    ns->elements   = s->elements;
    
    ((DENSE_STORAGE*)s->src)->count++;
    ns->src = s->src;

    return ns;
  }
}


/*
 * Does not free passed-in value! Different from list_storage_insert.
 */
void dense_storage_set(DENSE_STORAGE* s, SLICE* slice, void* val) {
  memcpy((char*)(s->elements) + dense_storage_pos(s, slice->coords) * DTYPE_SIZES[s->dtype], val, DTYPE_SIZES[s->dtype]);
}

///////////
// Tests //
///////////

/*
 * Do these two dense matrices have the same contents?
 *
 * TODO: Test the shape of the two matrices.
 * TODO: See if using memcmp is faster when the left- and right-hand matrices
 *				have the same dtype.
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

/////////////
// Utility //
/////////////

/*
 * Determine the linear array position (in elements of s) of some set of coordinates
 * (given by slice).
 */
size_t dense_storage_pos(DENSE_STORAGE* s, const size_t* coords) {
  size_t index, pos = 0;

  for (index = 0; index < s->rank; ++index)
    pos += (coords[i] + s->offset[i]) * s->stride[i];

  return pos;
}


/*
 * Calculate the stride length.
 */
size_t* dense_storage_stride(size_t* shape, size_t rank) {
  size_t i, j;
  size_t* stride = calloc(sizeof(*shape), rank);

  NM_CHECK_ALLOC(stride);

  for (i = 0; i < rank; ++i) {
    stride[i] = 1;
    for (j = i+1; j < rank; ++j) {
      stride[i] *= shape[j];
    }
  }

  return stride;
}


/*
 * Recursive slicing for N-dimensional matrix.
 */
void dense_storage_slice_copy(
    DENSE_STORAGE *src, DENSE_STORAGE *dest,
    size_t* lengths,
    size_t psrc, size_t pdest,
    size_t n)
{
  size_t index;

  if (src->rank - n > 1) {
    for (index = 0; index < lengths[n]; ++index) {
      dense_storage_slice_copy(src, dest, lengths,
                                    psrc + src->stride[n]*i, pdest + dest->stride[n]*i,
                                    n + 1);
    }
  } else {
    memcpy((char*)dest->elements + pdest*DTYPE_SIZES[dest->dtype],
        (char*)src->elements + psrc*DTYPE_SIZES[src->dtype],
        src->shape[n]*DTYPE_SIZES[dest->dtype]);
  }

}

/////////////////////////
// Copying and Casting //
/////////////////////////

/*
 * Documentation goes here.
 */
DENSE_STORAGE* dense_storage_cast_copy(const DENSE_STORAGE* rhs, dtype_t new_dtype) {
	LR_DTYPE_TEMPLATE_TABLE(dense_storage_cast_copy_template, DENSE_STORAGE*, const DENSE_STORAGE*, dtype_t);
	
	return ttable[new_dtype][rhs->dtype](rhs, new_dtype);
}

/*
 * Documentation goes here.
 */
DENSE_STORAGE* dense_storage_copy(const DENSE_STORAGE* rhs) {
  DENSE_STORAGE* lhs;

  size_t  count = storage_count_max_elements(rhs->rank, rhs->shape);
  size_t* shape = ALLOC_N(size_t, rhs->rank);
  NM_CHECK_ALLOC(shape);
  memcpy(shape, rhs->shape, sizeof(*shape) * rhs->rank);

  lhs = dense_storage_create(rhs->dtype, shape, rhs->rank, NULL, 0);

	// Ensure that allocation worked before copying.
  if (lhs && count) {
    if (rhs == rhs->ref) // not a reference
      memcpy(lhs->elements, rhs->elements, DTYPE_SIZES[rhs->dtype] * count);
    else // slice whole matrix
      dense_storage_slice_copy(lhs, rhs->src, rhs->shape, dense_storage_pos(rhs->src, rhs->offset), 0, 0);
  }

  return lhs;
}

/////////////////////////
// Templated Functions //
/////////////////////////

template <typename DType, typename NewDType>
DENSE_STORAGE* dense_storage_cast_copy_template(const DENSE_STORAGE* rhs, dtype_t new_dtype) {
  size_t  count = storage_count_max_elements(rhs->rank, rhs->shape);
  size_t* shape = ALLOC_N(size_t, rhs->rank);
  NM_CHECK_ALLOC(shape);
  memcpy(shape, rhs->shape, sizeof(*shape) * rhs->rank);
  
  DType*		rhs_els = (DType*)rhs->elements;
  NewDType* lhs_els;
  
  DENSE_STORAGE *lhs, *tmp;
	
  lhs			= dense_storage_create(new_dtype, shape, rhs->rank, NULL, 0);
  lhs_els	= (NewDType*)lhs->elements;

	// Ensure that allocation worked before copying.
  if (lhs && count) {
    if (NM_STORAGE(self)->src != NM_STORAGE(self)) {
      tmp = dense_storage_copy(rhs);
      while (count-- > 0)         lhs_els[count] = tmp->elements[count];
      dense_storage_delete(tmp);
    } else {
    	while (count-- > 0)     		lhs_els[count] = rhs_els[count];
    }
  }
	
  return lhs;
}

template <typename LDType, typename RDType>
bool dense_storage_eqeq_template(const DENSE_STORAGE* left, const DENSE_STORAGE* right) {
	int index;

  /* FIXME: Very strange behavior! The GC calls the method directly with non-initialized data. */
  if (left->rank != right->rank) return false;

	LDType* left_elements	  = (LDType*)left->elements;
	RDType* right_elements	= (RDType*)right->elements;
	
	for (index = storage_count_max_elements(left->rank, left->shape); index-- > 0;) {
		if (left_elements[index] != right_elements[index]) return false;
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

