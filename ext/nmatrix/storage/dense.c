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

#include "dense.h"
#include "list.h"
#include "yale.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

extern bool				(*ElemEqEq[NM_TYPES][2])(const void*, const void*, const int, const int);
extern const int	nm_sizeof[NM_TYPES];

/*
 * Forward Declarations
 */

static void dense_storage_cast_copy_list_contents(void* lhs, const LIST* rhs, void* default_val, dtype_t l_dtype, dtype_t r_dtype,
	size_t* pos, const size_t* shape, size_t rank, size_t max_elements, size_t recursions);

static void dense_storage_cast_copy_list_default(void* lhs, void* default_val, dtype_t l_dtype, dtype_t r_dtype,
	size_t* pos, const size_t* shape, size_t rank, size_t max_elements, size_t recursions);

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
  s->offset     = calloc(sizeof(size_t),rank);
  s->count      = 1;
  s->src        = s;
	
	count         = count_dense_storage_elements(s);

  if (elements_length == count) {
  	s->elements = elements;
  	
  } else {
    s->elements = ALLOC_N(char, nm_sizeof[dtype]*count);

    if (elements_length > 0) {
      // Repeat elements over and over again until the end of the matrix.
      for (i = 0; i < count; i += elements_length) {
        
        if (i + elements_length > count) {
        	copy_length = count - i;
        }
        
        memcpy((char*)(s->elements)+i*nm_sizeof[dtype], (char*)(elements)+(i % elements_length)*nm_sizeof[dtype], copy_length*nm_sizeof[dtype]);
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
void dense_storage_mark(void* m) {
  size_t i;
  DENSE_STORAGE* storage;

  if (m) {
    storage = (DENSE_STORAGE*)(((NMATRIX*)m)->storage);
    
    if (storage && storage->dtype == NM_ROBJ) {
      for (i = dense_storage_count_elements(storage); i-- > 0;) {
        rb_gc_mark(*((VALUE*)((char*)(storage->elements) + i*nm_sizeof[NM_ROBJ])));
      }
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
    return (char*)(s->elements) + dense_storage_pos(s, slice) * nm_sizeof[s->dtype];
    
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
  memcpy((char*)(s->elements) + dense_storage_pos(s, slice) * nm_sizeof[s->dtype], val, nm_sizeof[s->dtype]);
}

///////////
// Tests //
///////////

/////////////
// Utility //
/////////////

/*
 * Calculate the number of elements in the dense storage structure, based on
 * shape and rank.
 */
size_t dense_storage_count_elements(const DENSE_STORAGE* s) {
  size_t i;
  size_t count = 1;
  
  for (i = s->rank; i-- > 0;) {
  	count *= s->shape[i];
  }
  
  return count;
}

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
  
  size_t count = count_dense_storage_elements(rhs), p;
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
    memcpy(lhs->elements, rhs->elements, nm_sizeof[rhs->dtype] * count);
  }

  return lhs;
}

/*
 * Documentation goes here.
 */
DENSE_STORAGE* dense_storage_cast_copy(DENSE_STORAGE* rhs, dtype_t new_dtype) {
  DENSE_STORAGE* lhs;
  
  size_t count = count_dense_storage_elements(rhs), p;
  size_t* shape = ALLOC_N(size_t, rhs->rank);
  
  if (!shape) {
  	return NULL;
  }

  // Copy shape array.
  for (p = rhs->rank; p-- > 0;) {
  	shape[p] = rhs->shape[p];
  }

  lhs = dense_storage_create(new_dtype, shape, rhs->rank, NULL, 0);

	// Ensure that allocation worked before copying.
  if (lhs && count) {
    if (lhs->dtype == rhs->dtype) {
      memcpy(lhs->elements, rhs->elements, nm_sizeof[rhs->dtype] * count);
      
    } else {
      SetFuncs[lhs->dtype][rhs->dtype](count, lhs->elements, nm_sizeof[lhs->dtype], rhs->elements, nm_sizeof[rhs->dtype]);
    }
  }


  return lhs;
}

//////////////////////
// Helper Functions //
//////////////////////

/*
 * Copy list contents into dense recursively.
 */
static void dense_storage_cast_copy_list_contents(void* lhs, const LIST* rhs, void* default_val, dtype_t l_dtype, dtype_t r_dtype,size_t* pos, const size_t* shape, size_t rank, size_t max_elements, size_t recursions) {

  NODE *curr = rhs->first;
  int last_key = -1;
  size_t i = 0;

	for (i = shape[rank - 1 - recursions]; i-- > 0; ++(*pos) {

    if (!curr || (curr->key > (size_t)(last_key+1))) {
      //fprintf(stderr, "pos = %u, dim = %u, curr->key XX, last_key+1 = %d\t", *pos, shape[rank-1-recursions], last_key+1);
      
      if (recursions == 0) {
      	cast_copy_value_single((char*)lhs + (*pos)*nm_sizeof[l_dtype], default_val, l_dtype, r_dtype);
    		//fprintf(stderr, "zero\n");
      
      } else {
     		dense_storage_cast_copy_list_default(lhs, default_val, l_dtype, r_dtype, pos, shape, rank, max_elements, recursions-1);
    		//fprintf(stderr, "column of zeros\n");
			}

      ++last_key;
      
    } else {
      //fprintf(stderr, "pos = %u, dim = %u, curr->key = %u, last_key+1 = %d\t", *pos, shape[rank-1-recursions], curr->key, last_key+1);
      
      if (recursions == 0) {
      	cast_copy_value_single((char*)lhs + (*pos)*nm_sizeof[l_dtype], curr->val, l_dtype, r_dtype);
    		//fprintf(stderr, "zero\n");
      	
      } else {
      	dense_storage_cast_copy_list_default(lhs, curr->val, default_val, l_dtype, r_dtype, pos, shape, rank, max_elements, recursions-1);
    		//fprintf(stderr, "column of zeros\n");
      }

      last_key = curr->key;
      curr     = curr->next;
    }
  }
  
  --(*pos);
}

/*
 * Copy a set of default values into dense.
 */
static void dense_storage_cast_copy_list_default(void* lhs, void* default_val, dtype_t l_dtype, dtype_t r_dtype, size_t* pos, const size_t* shape, size_t rank, size_t max_elements, size_t recursions) {
  size_t i;

	for (i = shape[rank - 1 - recursions]; i-- > 0; ++(*pos)) {
    //fprintf(stderr, "default: pos = %u, dim = %u\t", *pos, shape[rank-1-recursions]);

    if (recursions == 0) {
    	cast_copy_value_single((char*)lhs + (*pos)*nm_sizeof[l_dtype], default_val, l_dtype, r_dtype);
    	//fprintf(stderr, "zero\n");
    	
    } else {
    	dense_storage_cast_copy_list_default(lhs, default_val, l_dtype, r_dtype, pos, shape, rank, max_elements, recursions-1);
    	//fprintf(stderr, "column of zeros\n");
  	}
  }
  
  --(*pos);
}

