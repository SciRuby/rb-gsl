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
    if(s->count <= 1) {
      free(s->shape);
      free(s->offset);
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
    ((DENSE_STORAGE*)s->src)->count--;
    free(s->shape);
    free(s->offset);
    free(s);
  }
}

/*
 * Documentation goes here.
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
 * Documentation goes here.
 */
void* dense_storage_get(DENSE_STORAGE* s, SLICE* slice) {
  DENSE_STORAGE *ns;

  if (slice->is_one_el) {
    return (char*)(s->elements) + dense_storage_pos(s, slice) * DTYPE_SIZES[s->dtype];
    
  } else {
    ns = ALLOC( DENSE_STORAGE );

    ns->rank       = s->rank;
    ns->shape      = slice->lens;
    ns->dtype      = s->dtype;
    ns->offset     = slice->coords;
    ns->elements   = s->elements;
    
    s->count++;
    ns->src = (void*)s;

    return ns;
  }
}


/*
 * Does not free passed-in value! Different from list_storage_insert.
 */
void dense_storage_set(DENSE_STORAGE* s, SLICE* slice, void* val) {
  memcpy((char*)(s->elements) + dense_storage_pos(s, slice) * DTYPE_SIZES[s->dtype], val, DTYPE_SIZES[s->dtype]);
}

///////////
// Tests //
///////////

/////////////
// Utility //
/////////////

/*
 * Documentation goes here.
 */
size_t dense_storage_pos(DENSE_STORAGE* s, SLICE* slice) {
  size_t k, l;
  size_t inner, outer = 0;
  
  for (k = s->rank; k-- > 0;) {
  	inner = slice->coords[k] + s->offset[k];
    
    for (l = k+1; l < s->rank; ++l) {
      inner *= ((DENSE_STORAGE*)s->src)->shape[l];
    }
    
    outer += inner;
  }
  
  return outer;
}

/////////////////////////
// Copying and Casting //
/////////////////////////

/*
 * Documentation goes here.
 */
DENSE_STORAGE* dense_storage_copy(DENSE_STORAGE* rhs) {
  DENSE_STORAGE* lhs;
  
  size_t  count = storage_count_max_elements(rhs->rank, rhs->shape), p;
  size_t* shape = ALLOC_N(size_t, rhs->rank);
  
  if (!shape) {
  	return NULL;
  }

  // copy shape array
  for (p = rhs->rank; p-- > 0;) {
    shape[p] = rhs->shape[p];
  }

  lhs = dense_storage_create(rhs->dtype, shape, rhs->rank, NULL, 0);

	// Ensure that allocation worked before copying.
  if (lhs && count) {
    memcpy(lhs->elements, rhs->elements, DTYPE_SIZES[rhs->dtype] * count);
  }

  return lhs;
}

