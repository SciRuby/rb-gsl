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
// #include "types.h"
#include "util/math.h"

#include "data/data.h"

#include "common.h"
#include "dense.h"

/*
 * Macros
 */
#ifndef NM_CHECK_ALLOC
#define NM_CHECK_ALLOC(x) if (!x) rb_raise(rb_eNoMemError, "insufficient memory");
#endif

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

static size_t* dense_storage_stride(size_t* shape, size_t rank);
static void dense_storage_slice_copy(DENSE_STORAGE *dest, const DENSE_STORAGE *src, size_t* lengths, size_t psrc, size_t pdest, size_t n);

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
  DENSE_STORAGE* s = ALLOC( DENSE_STORAGE );

  s->rank       = rank;
  s->shape      = shape;
  s->dtype      = dtype;
  s->offset     = ALLOC_N(size_t, rank);

  memset(s->offset, 0, sizeof(size_t)*rank);

  s->stride     = dense_storage_stride(shape, rank);
  s->count      = 1;
  s->src        = s;
	
	size_t count  = storage_count_max_elements(s);

  if (elements_length == count) {
  	s->elements = elements;
  	
  } else {
    s->elements = ALLOC_N(char, DTYPE_SIZES[dtype]*count);
    size_t copy_length = elements_length;

    if (elements_length > 0) {
      // Repeat elements over and over again until the end of the matrix.
      for (size_t i = 0; i < count; i += elements_length) {

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
 * Destructor for dense storage
 */
void dense_storage_delete(STORAGE* s) {
  // Sometimes Ruby passes in NULL storage for some reason (probably on copy construction failure).
  if (s) {
    DENSE_STORAGE* storage = reinterpret_cast<DENSE_STORAGE*>(s);
    if(storage->count-- == 1) {
      free(storage->shape);
      free(storage->offset);
      free(storage->stride);
      free(storage->elements);
      free(storage);
    }
  }
}

/*
 * Destructor for dense storage references (slicing).
 */
void dense_storage_delete_ref(STORAGE* s) {
  // Sometimes Ruby passes in NULL storage for some reason (probably on copy construction failure).
  if (s) {
    DENSE_STORAGE* storage = reinterpret_cast<DENSE_STORAGE*>(s);
    dense_storage_delete( reinterpret_cast<STORAGE*>(storage->src) );
    free(storage->shape);
    free(storage->offset);
    free(storage);
  }
}

/*
 * Mark values in a dense matrix for garbage collection. This may not be necessary -- further testing required.
 */
void dense_storage_mark(void* storage_base) {
  DENSE_STORAGE* storage = reinterpret_cast<DENSE_STORAGE*>(storage_base);
  size_t index;
	
	VALUE* els = (VALUE*)storage->elements;
	
  if (storage && storage->dtype == RUBYOBJ) {
  	for (index = storage_count_max_elements(storage); index-- > 0;) {
      rb_gc_mark(els[index]);
    }
  }
}

///////////////
// Accessors //
///////////////

/*
 * Get a slice or one element, using copying.
 *
 * FIXME: Template the first condition.
 */
void* dense_storage_get(STORAGE* storage, SLICE* slice) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)storage;
  DENSE_STORAGE* ns;
  size_t count;

  if (slice->is_one_el)
    return (char*)(s->elements) + dense_storage_pos(s, slice->coords) * DTYPE_SIZES[s->dtype];
  else { // Make references
    ns = ALLOC( DENSE_STORAGE );

    NM_CHECK_ALLOC(ns);

    ns->rank       = s->rank;
    ns->shape      = slice->lengths;
    ns->dtype      = s->dtype;

    ns->offset     = ALLOC_N(size_t, ns->rank);
    NM_CHECK_ALLOC(ns->offset);

    ns->stride     = dense_storage_stride(ns->shape, ns->rank);
    ns->count      = 1;
    ns->src        = ns;

    count          = storage_count_max_elements(s);

    ns->elements   = ALLOC_N(char, DTYPE_SIZES[ns->dtype] * count);
    NM_CHECK_ALLOC(ns->elements);

    dense_storage_slice_copy(ns, s, slice->lengths, dense_storage_pos(s, slice->coords), 0, 0);
    return ns;
  }
}



/*
 * Get a slice or one element by reference (no copy).
 *
 * FIXME: Template the first condition.
 */
void* dense_storage_ref(STORAGE* storage, SLICE* slice) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)storage;

  if (slice->is_one_el)
    return (char*)(s->elements) + dense_storage_pos(s, slice->coords) * DTYPE_SIZES[s->dtype];
    
  else {
    DENSE_STORAGE* ns = ALLOC( DENSE_STORAGE );
    NM_CHECK_ALLOC(ns);

    ns->rank       = s->rank;
    ns->dtype      = s->dtype;

    ns->offset     = ALLOC_N(size_t, ns->rank);
    NM_CHECK_ALLOC(ns->offset);

    ns->shape      = ALLOC_N(size_t, ns->rank);
    NM_CHECK_ALLOC(ns->shape);

    for (size_t i = 0; i < ns->rank; ++i) {
      ns->offset[i] = slice->coords[i] + s->offset[i];
      ns->shape[i]  = slice->lengths[i];
    }

    ns->stride     = s->stride;
    ns->elements   = s->elements;
    
    reinterpret_cast<DENSE_STORAGE*>((DENSE_STORAGE*)s->src)->count++;
    ns->src = s->src;

    return ns;
  }
}


/*
 * Does not free passed-in value! Different from list_storage_insert.
 */
void dense_storage_set(STORAGE* storage, SLICE* slice, void* val) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)storage;
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
bool dense_storage_eqeq(const STORAGE* left, const STORAGE* right) {
	LR_DTYPE_TEMPLATE_TABLE(dense_storage_eqeq_template, bool, const DENSE_STORAGE*, const DENSE_STORAGE*);
	
	return ttable[left->dtype][right->dtype]((const DENSE_STORAGE*)left, (const DENSE_STORAGE*)right);
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
size_t dense_storage_pos(const DENSE_STORAGE* s, const size_t* coords) {
  size_t pos = 0;

  for (size_t i = 0; i < s->rank; ++i)
    pos += (coords[i] + s->offset[i]) * s->stride[i];

  return pos;
}


/*
 * Calculate the stride length.
 */
static size_t* dense_storage_stride(size_t* shape, size_t rank) {
  size_t i, j;
  size_t* stride = ALLOC_N(size_t, rank);

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
static void dense_storage_slice_copy(DENSE_STORAGE *dest, const DENSE_STORAGE *src, size_t* lengths, size_t psrc, size_t pdest, size_t n) {
  if (src->rank - n > 1) {
    for (size_t i = 0; i < lengths[n]; ++i) {
      dense_storage_slice_copy(dest, src, lengths,
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
 * Copy dense storage, changing dtype if necessary.
 */
STORAGE* dense_storage_cast_copy(const STORAGE* rhs, dtype_t new_dtype) {
	LR_DTYPE_TEMPLATE_TABLE(dense_storage_cast_copy_template, DENSE_STORAGE*, const DENSE_STORAGE*, dtype_t);
	
	return (STORAGE*)ttable[new_dtype][rhs->dtype]( (DENSE_STORAGE*)rhs, new_dtype );
}

/*
 * Copy dense storage without a change in dtype.
 */
DENSE_STORAGE* dense_storage_copy(const DENSE_STORAGE* rhs) {
  DENSE_STORAGE* lhs;

  size_t  count = storage_count_max_elements(rhs);
  size_t* shape = ALLOC_N(size_t, rhs->rank);
  NM_CHECK_ALLOC(shape);

  // copy shape
  for (size_t i = 0; i < rhs->rank; ++i)
    shape[i] = rhs->shape[i];

  lhs = dense_storage_create(rhs->dtype, shape, rhs->rank, NULL, 0);

	// Ensure that allocation worked before copying.
  if (lhs && count) {
    if (rhs == rhs->src) // not a reference
      memcpy(lhs->elements, rhs->elements, DTYPE_SIZES[rhs->dtype] * count);
    else // slice whole matrix
      dense_storage_slice_copy(lhs,
                               reinterpret_cast<const DENSE_STORAGE*>(rhs->src),
                               rhs->shape,
                               dense_storage_pos(reinterpret_cast<const DENSE_STORAGE*>(rhs->src), rhs->offset),
                               0,
                               0);
  }

  return lhs;
}

/////////////////////////
// Templated Functions //
/////////////////////////

template <typename LDType, typename RDType>
DENSE_STORAGE* dense_storage_cast_copy_template(const DENSE_STORAGE* rhs, dtype_t new_dtype) {
  size_t  count = storage_count_max_elements(rhs);
  size_t* shape = ALLOC_N(size_t, rhs->rank);
  NM_CHECK_ALLOC(shape);

  memcpy(shape, rhs->shape, sizeof(*shape) * rhs->rank);

  DENSE_STORAGE* lhs			= dense_storage_create(new_dtype, shape, rhs->rank, NULL, 0);
  RDType*	rhs_els = reinterpret_cast<RDType*>(rhs->elements);
  LDType* lhs_els	= reinterpret_cast<LDType*>(lhs->elements);

	// Ensure that allocation worked before copying.
  if (lhs && count) {
    if (rhs->src != rhs) {
      DENSE_STORAGE* tmp = dense_storage_copy(rhs);
      RDType* tmp_els    = reinterpret_cast<RDType*>(tmp->elements);
      while (count-- > 0)         lhs_els[count] = tmp_els[count];
      dense_storage_delete(reinterpret_cast<STORAGE*>(tmp));
    } else {
    	while (count-- > 0)     		lhs_els[count] = rhs_els[count];
    }
  }
	
  return lhs;
}

template <typename LDType, typename RDType>
bool dense_storage_eqeq_template(const DENSE_STORAGE* left, const DENSE_STORAGE* right) {
  size_t index;
  
  /* FIXME: Very strange behavior! The GC calls the method directly with non-initialized data. */
  if (left->rank != right->rank) return false;

	LDType* left_elements	  = (LDType*)left->elements;
	RDType* right_elements	= (RDType*)right->elements;
	
	for (index = storage_count_max_elements(left); index-- > 0;) {
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


template <typename DType>
static DENSE_STORAGE* dense_storage_matrix_multiply_template(STORAGE_PAIR casted_storage, size_t* resulting_shape, bool vector) {
  DENSE_STORAGE *left  = (DENSE_STORAGE*)(casted_storage.left),
                *right = (DENSE_STORAGE*)(casted_storage.right);

  // Create result storage.
  DENSE_STORAGE* result = dense_storage_create(left->dtype, resulting_shape, 2, NULL, 0);

  DType *pAlpha = new DType(1),
        *pBeta  = new DType(0);

  // Do the multiplication
  bool succ;

  if (vector) succ = gemv<DType>(CblasNoTrans, left->shape[0], left->shape[1], pAlpha,
                                 reinterpret_cast<DType*>(left->elements), left->shape[1],
                                 reinterpret_cast<DType*>(right->elements), 1, pBeta,
                                 reinterpret_cast<DType*>(result->elements), 1);
  else        succ = gemm<DType>(CblasNoTrans, CblasNoTrans, right->shape[1], left->shape[0], left->shape[1], pAlpha,
                                 reinterpret_cast<DType*>(right->elements), right->shape[1],
                                 reinterpret_cast<DType*>(left->elements), left->shape[1], pBeta,
                                 reinterpret_cast<DType*>(result->elements), result->shape[1]);

  delete pAlpha;
  delete pBeta;

  if (!succ)
    rb_raise(rb_eStandardError, "gemm/gemv failed for an unknown reason");

  return result;
}

STORAGE* dense_storage_matrix_multiply(STORAGE_PAIR casted_storage, size_t* resulting_shape, bool vector) {
  NAMED_DTYPE_TEMPLATE_TABLE(ttable, dense_storage_matrix_multiply_template, DENSE_STORAGE*, STORAGE_PAIR casted_storage, size_t* resulting_shape, bool vector);

  return (STORAGE*)ttable[reinterpret_cast<DENSE_STORAGE*>(casted_storage.left)->dtype](casted_storage, resulting_shape, vector);
}
