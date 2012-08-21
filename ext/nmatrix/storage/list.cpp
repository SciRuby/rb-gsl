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
// == list.c
//
// List-of-lists n-dimensional matrix storage. Uses singly-linked
// lists.

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
#include "list.h"

#include "util/math.h"
#include "util/sl_list.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

namespace nm { namespace list_storage {

/*
 * Forward Declarations
 */

template <typename LDType, typename RDType>
static LIST_STORAGE* cast_copy(const LIST_STORAGE* rhs, dtype_t new_dtype);

template <typename LDType, typename RDType>
static bool eqeq(const LIST_STORAGE* left, const LIST_STORAGE* right);

template <typename LDType, typename RDType>
static void* ew_add(LIST* dest, const LIST* left, const void* l_default, const LIST* right, const void* r_default, const size_t* shape, size_t rank);

template <typename LDType, typename RDType>
static void* ew_subtract(LIST* dest, const LIST* left, const void* l_default, const LIST* right, const void* r_default, const size_t* shape, size_t rank);

template <typename LDType, typename RDType>
static void* ew_multiply(LIST* dest, const LIST* left, const void* l_default, const LIST* right, const void* r_default, const size_t* shape, size_t rank);

template <typename LDType, typename RDType>
static void* ew_divide(LIST* dest, const LIST* left, const void* l_default, const LIST* right, const void* r_default, const size_t* shape, size_t rank);

template <typename LDType, typename RDType>
static void ew_add_prime(LIST* dest, LDType d_default, const LIST* left, LDType l_default, const LIST* right, RDType r_default, const size_t* shape, size_t last_level, size_t level);

template <typename LDType, typename RDType>
static void ew_subtract_prime(LIST* dest, LDType d_default, const LIST* left, LDType l_default, const LIST* right, RDType r_default, const size_t* shape, size_t last_level, size_t level);

template <typename LDType, typename RDType>
static void ew_multiply_prime(LIST* dest, LDType d_default, const LIST* left, LDType l_default, const LIST* right, RDType r_default, const size_t* shape, size_t last_level, size_t level);

template <typename LDType, typename RDType>
static void ew_divide_prime(LIST* dest, LDType d_default, const LIST* left, LDType l_default, const LIST* right, RDType r_default, const size_t* shape, size_t last_level, size_t level);

} // end of namespace list_storage

extern "C" {
/*
 * Functions
 */


////////////////
// Lifecycle //
///////////////

/*
 * Creates a list-of-lists(-of-lists-of-lists-etc) storage framework for a
 * matrix.
 *
 * Note: The pointers you pass in for shape and init_val become property of our
 * new storage. You don't need to free them, and you shouldn't re-use them.
 */
LIST_STORAGE* list_storage_create(dtype_t dtype, size_t* shape, size_t rank, void* init_val) {
  LIST_STORAGE* s;

  s = ALLOC( LIST_STORAGE );

  s->rank  = rank;
  s->shape = shape;
  s->dtype = dtype;

  s->rows  = list::create();

  s->default_val = init_val;

  return s;
}

/*
 * Documentation goes here.
 */
void list_storage_delete(STORAGE* s) {
  if (s) {
    LIST_STORAGE* storage = (LIST_STORAGE*)s;

    list::del( storage->rows, storage->rank - 1 );

    free(storage->shape);
    free(storage->default_val);
    free(s);
  }
}

/*
 * Documentation goes here.
 */
void list_storage_mark(void* storage_base) {
  LIST_STORAGE* storage = (LIST_STORAGE*)storage_base;

  if (storage && storage->dtype == RUBYOBJ) {
    rb_gc_mark(*((VALUE*)(storage->default_val)));
    list::mark(storage->rows, storage->rank - 1);
  }
}

///////////////
// Accessors //
///////////////


/*
 * Documentation goes here.
 */
void* list_storage_get(STORAGE* storage, SLICE* slice) {
  //LIST_STORAGE* s = (LIST_STORAGE*)storage;
  rb_raise(rb_eNotImpError, "This type of slicing not supported yet");
}


/*
 * Get the contents of some set of coordinates. Note: Does not make a copy!
 * Don't free!
 */
void* list_storage_ref(STORAGE* storage, SLICE* slice) {
  LIST_STORAGE* s = (LIST_STORAGE*)storage;
  size_t r;
  NODE*  n;
  LIST*  l = s->rows;

  for (r = s->rank; r > 1; --r) {
    n = list::find(l, slice->coords[s->rank - r]);
    if (n)  l = reinterpret_cast<LIST*>(n->val);
    else return s->default_val;
  }

  n = list::find(l, slice->coords[s->rank - r]);
  if (n) return n->val;
  else   return s->default_val;
}

/*
 * Documentation goes here.
 *
 * TODO: Allow this function to accept an entire row and not just one value -- for slicing
 */
void* list_storage_insert(STORAGE* storage, SLICE* slice, void* val) {
  LIST_STORAGE* s = (LIST_STORAGE*)storage;
  // Pretend ranks = 2
  // Then coords is going to be size 2
  // So we need to find out if some key already exists
  size_t r;
  NODE*  n;
  LIST*  l = s->rows;

  // drill down into the structure
  for (r = s->rank; r > 1; --r) {
    n = list::insert(l, false, slice->coords[s->rank - r], list::create());
    l = reinterpret_cast<LIST*>(n->val);
  }

  n = list::insert(l, true, slice->coords[s->rank - r], val);
  return n->val;
}

/*
 * Documentation goes here.
 *
 * TODO: Speed up removal.
 */
void* list_storage_remove(STORAGE* storage, SLICE* slice) {
  LIST_STORAGE* s = (LIST_STORAGE*)storage;
  int r;
  NODE  *n = NULL;
  LIST*  l = s->rows;
  void*  rm = NULL;

  // keep track of where we are in the traversals
  NODE** stack = ALLOCA_N( NODE*, s->rank - 1 );

  for (r = (int)(s->rank); r > 1; --r) {
  	// does this row exist in the matrix?
    n = list::find(l, slice->coords[s->rank - r]);

    if (!n) {
    	// not found
      free(stack);
      return NULL;
      
    } else {
    	// found
      stack[s->rank - r]    = n;
      l                     = reinterpret_cast<LIST*>(n->val);
    }
  }

  rm = list::remove(l, slice->coords[s->rank - r]);

  // if we removed something, we may now need to remove parent lists
  if (rm) {
    for (r = (int)(s->rank) - 2; r >= 0; --r) {
    	// walk back down the stack
      
      if (((LIST*)(stack[r]->val))->first == NULL)
        free(list::remove(reinterpret_cast<LIST*>(stack[r]->val), slice->coords[r]));
      else break; // no need to continue unless we just deleted one.

    }
  }

  return rm;
}

///////////
// Tests //
///////////

/*
 * Comparison of contents for list storage.
 */
bool list_storage_eqeq(const STORAGE* left, const STORAGE* right) {
	NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::list_storage::eqeq, bool, const LIST_STORAGE* left, const LIST_STORAGE* right);

	return ttable[left->dtype][right->dtype]((const LIST_STORAGE*)left, (const LIST_STORAGE*)right);
}

//////////
// Math //
//////////

// TODO: Condense all of the element-wise operations.

/*
 * Element-wise addition for list storage.
 */
STORAGE* list_storage_ew_add(const STORAGE* left, const STORAGE* right) {
	LR_DTYPE_TEMPLATE_TABLE(nm::list_storage::ew_add, void*, LIST*, const LIST*, const void*, const LIST*, const void*, const size_t*, size_t);
	
	dtype_t new_dtype = Upcast[left->dtype][right->dtype];
	
	const LIST_STORAGE* l = reinterpret_cast<const LIST_STORAGE*>(left),
										* r = reinterpret_cast<const LIST_STORAGE*>(right);
	
	LIST_STORAGE* new_l = NULL;
	
	// Allocate a new shape array for the resulting matrix.
	size_t* new_shape = (size_t*)calloc(l->rank, sizeof(size_t));
	memcpy(new_shape, left->shape, sizeof(size_t) * l->rank);
	
	// Create the result matrix.
	LIST_STORAGE* result = list_storage_create(new_dtype, new_shape, left->rank, NULL); 
	
	/*
	 * Call the templated elementwise multiplication function and set the default
	 * value for the resulting matrix.
	 */
	if (new_dtype != left->dtype) {
		// Upcast the left-hand side if necessary.
		new_l = reinterpret_cast<LIST_STORAGE*>(list_storage_cast_copy(l, new_dtype));
		
		result->default_val =
			ttable[left->dtype][right->dtype](result->rows, new_l->rows, new_l->default_val, r->rows, r->default_val, result->shape, result->rank);
		
		// Delete the temporary left-hand side matrix.
		list_storage_delete(reinterpret_cast<STORAGE*>(new_l));
			
	} else {
		result->default_val =
			ttable[left->dtype][right->dtype](result->rows, l->rows, l->default_val, r->rows, r->default_val, result->shape, result->rank);
	}
	
	return result;
}

/*
 * Element-wise subtraction for list storage.
 */
STORAGE* list_storage_ew_subtract(const STORAGE* left, const STORAGE* right) {
	LR_DTYPE_TEMPLATE_TABLE(nm::list_storage::ew_subtract, void*, LIST*, const LIST*, const void*, const LIST*, const void*, const size_t*, size_t);
	
	dtype_t new_dtype = Upcast[left->dtype][right->dtype];
	
	const LIST_STORAGE* l = reinterpret_cast<const LIST_STORAGE*>(left),
										* r = reinterpret_cast<const LIST_STORAGE*>(right);
	
	LIST_STORAGE* new_l = NULL;
	
	// Allocate a new shape array for the resulting matrix.
	size_t* new_shape = (size_t*)calloc(l->rank, sizeof(size_t));
	memcpy(new_shape, left->shape, sizeof(size_t) * l->rank);
	
	// Create the result matrix.
	LIST_STORAGE* result = list_storage_create(new_dtype, new_shape, left->rank, NULL); 
	
	/*
	 * Call the templated elementwise multiplication function and set the default
	 * value for the resulting matrix.
	 */
	if (new_dtype != left->dtype) {
		// Upcast the left-hand side if necessary.
		new_l = reinterpret_cast<LIST_STORAGE*>(list_storage_cast_copy(l, new_dtype));
		
		result->default_val =
			ttable[left->dtype][right->dtype](result->rows, new_l->rows, new_l->default_val, r->rows, r->default_val, result->shape, result->rank);
		
		// Delete the temporary left-hand side matrix.
		list_storage_delete(reinterpret_cast<STORAGE*>(new_l));
			
	} else {
		result->default_val =
			ttable[left->dtype][right->dtype](result->rows, l->rows, l->default_val, r->rows, r->default_val, result->shape, result->rank);
	}
	
	return result;
}

/*
 * Element-wise multiplication for list storage.
 */
STORAGE* list_storage_ew_multiply(const STORAGE* left, const STORAGE* right) {
	LR_DTYPE_TEMPLATE_TABLE(nm::list_storage::ew_multiply, void*, LIST*, const LIST*, const void*, const LIST*, const void*, const size_t*, size_t);
	
	dtype_t new_dtype = Upcast[left->dtype][right->dtype];
	
	const LIST_STORAGE* l = reinterpret_cast<const LIST_STORAGE*>(left),
										* r = reinterpret_cast<const LIST_STORAGE*>(right);
	
	LIST_STORAGE* new_l = NULL;
	
	// Allocate a new shape array for the resulting matrix.
	size_t* new_shape = (size_t*)calloc(l->rank, sizeof(size_t));
	memcpy(new_shape, left->shape, sizeof(size_t) * l->rank);
	
	// Create the result matrix.
	LIST_STORAGE* result = list_storage_create(new_dtype, new_shape, left->rank, NULL); 
	
	/*
	 * Call the templated elementwise multiplication function and set the default
	 * value for the resulting matrix.
	 */
	if (new_dtype != left->dtype) {
		// Upcast the left-hand side if necessary.
		new_l = reinterpret_cast<LIST_STORAGE*>(list_storage_cast_copy(l, new_dtype));
		
		result->default_val =
			ttable[left->dtype][right->dtype](result->rows, new_l->rows, new_l->default_val, r->rows, r->default_val, result->shape, result->rank);
		
		// Delete the temporary left-hand side matrix.
		list_storage_delete(reinterpret_cast<STORAGE*>(new_l));
			
	} else {
		result->default_val =
			ttable[left->dtype][right->dtype](result->rows, l->rows, l->default_val, r->rows, r->default_val, result->shape, result->rank);
	}
	
	return result;
}

/*
 * Element-wise division for list storage.
 */
STORAGE* list_storage_ew_divide(const STORAGE* left, const STORAGE* right) {
	LR_DTYPE_TEMPLATE_TABLE(nm::list_storage::ew_divide, void*, LIST*, const LIST*, const void*, const LIST*, const void*, const size_t*, size_t);
	
	dtype_t new_dtype = Upcast[left->dtype][right->dtype];
	
	const LIST_STORAGE* l = reinterpret_cast<const LIST_STORAGE*>(left),
										* r = reinterpret_cast<const LIST_STORAGE*>(right);
	
	LIST_STORAGE* new_l = NULL;
	
	// Allocate a new shape array for the resulting matrix.
	size_t* new_shape = (size_t*)calloc(l->rank, sizeof(size_t));
	memcpy(new_shape, left->shape, sizeof(size_t) * l->rank);
	
	// Create the result matrix.
	LIST_STORAGE* result = list_storage_create(new_dtype, new_shape, left->rank, NULL); 
	
	/*
	 * Call the templated elementwise multiplication function and set the default
	 * value for the resulting matrix.
	 */
	if (new_dtype != left->dtype) {
		// Upcast the left-hand side if necessary.
		new_l = reinterpret_cast<LIST_STORAGE*>(list_storage_cast_copy(l, new_dtype));
		
		result->default_val =
			ttable[left->dtype][right->dtype](result->rows, new_l->rows, new_l->default_val, r->rows, r->default_val, result->shape, result->rank);
		
		// Delete the temporary left-hand side matrix.
		list_storage_delete(reinterpret_cast<STORAGE*>(new_l));
			
	} else {
		result->default_val =
			ttable[left->dtype][right->dtype](result->rows, l->rows, l->default_val, r->rows, r->default_val, result->shape, result->rank);
	}
	
	return result;
}


/*
 * List storage matrix multiplication.
 */
STORAGE* list_storage_matrix_multiply(const STORAGE_PAIR& casted_storage, size_t* resulting_shape, bool vector) {
  free(resulting_shape);
  rb_raise(rb_eNotImpError, "multiplication not implemented for list-of-list matrices");
  return NULL;
  //DTYPE_TEMPLATE_TABLE(dense_storage::matrix_multiply, NMATRIX*, STORAGE_PAIR, size_t*, bool);

  //return ttable[reinterpret_cast<DENSE_STORAGE*>(casted_storage.left)->dtype](casted_storage, resulting_shape, vector);
}

/////////////
// Utility //
/////////////

/*
 * Recursively count the non-zero elements in a list storage object.
 */
size_t list_storage_count_elements_r(const LIST* l, size_t recursions) {
  size_t count = 0;
  NODE* curr = l->first;
  
  if (recursions) {
    while (curr) {
      count += list_storage_count_elements_r(reinterpret_cast<const LIST*>(curr->val), recursions - 1);
      curr   = curr->next;
    }
    
  } else {
    while (curr) {
      ++count;
      curr = curr->next;
    }
  }
  
  return count;
}

/*
 * Count non-diagonal non-zero elements.
 */
size_t list_storage_count_nd_elements(const LIST_STORAGE* s) {
  NODE *i_curr, *j_curr;
  size_t count = 0;
  
  if (s->rank != 2) {
  	rb_raise(rb_eNotImpError, "non-diagonal element counting only defined for rank = 2");
  }

  for (i_curr = s->rows->first; i_curr; i_curr = i_curr->next) {
    for (j_curr = ((LIST*)(i_curr->val))->first; j_curr; j_curr = j_curr->next) {
      if (i_curr->key != j_curr->key) {
      	++count;
      }
    }
  }
  
  return count;
}

/////////////////////////
// Copying and Casting //
/////////////////////////

/*
 * List storage copy constructor C access.
 */
STORAGE* list_storage_cast_copy(const STORAGE* rhs, dtype_t new_dtype) {
  NAMED_LR_DTYPE_TEMPLATE_TABLE(ttable, nm::list_storage::cast_copy, LIST_STORAGE*, const LIST_STORAGE* rhs, dtype_t new_dtype);

  return (STORAGE*)ttable[new_dtype][rhs->dtype]((LIST_STORAGE*)rhs, new_dtype);
}


/*
 * List storage copy constructor for transposing.
 */
STORAGE* list_storage_copy_transposed(const STORAGE* rhs_base) {
  rb_raise(rb_eNotImpError, "list storage transpose not yet implemented");
  return NULL;
}


} // end of extern "C" block


/////////////////////////
// Templated Functions //
/////////////////////////

namespace list_storage {

/*
 * List storage copy constructor for changing dtypes.
 */
template <typename LDType, typename RDType>
static LIST_STORAGE* cast_copy(const LIST_STORAGE* rhs, dtype_t new_dtype) {

  // allocate and copy shape
  size_t* shape = ALLOC_N(size_t, rhs->rank);
  memcpy(shape, rhs->shape, rhs->rank * sizeof(size_t));

  // copy default value
  LDType* default_val = ALLOC_N(LDType, 1);
  *default_val = *reinterpret_cast<RDType*>(rhs->default_val);

  LIST_STORAGE* lhs = list_storage_create(new_dtype, shape, rhs->rank, default_val);
  lhs->rows         = list::create();
  list::cast_copy_contents<LDType, RDType>(lhs->rows, rhs->rows, rhs->rank - 1);

  return lhs;
}


/*
 * Do these two dense matrices of the same dtype have exactly the same
 * contents?
 */
template <typename LDType, typename RDType>
bool eqeq(const LIST_STORAGE* left, const LIST_STORAGE* right) {

  // in certain cases, we need to keep track of the number of elements checked.
  size_t num_checked  = 0,

	max_elements = storage_count_max_elements(left);

  if (!left->rows->first) {
    // Easy: both lists empty -- just compare default values
    if (!right->rows->first) {
    	return *reinterpret_cast<LDType*>(left->default_val) == *reinterpret_cast<RDType*>(right->default_val);
    	
    } else if (!list::eqeq_value<RDType,LDType>(right->rows, reinterpret_cast<LDType*>(left->default_val), left->rank-1, num_checked)) {
    	// Left empty, right not empty. Do all values in right == left->default_val?
    	return false;
    	
    } else if (num_checked < max_elements) {
    	// If the matrix isn't full, we also need to compare default values.
    	return *reinterpret_cast<LDType*>(left->default_val) == *reinterpret_cast<RDType*>(right->default_val);
    }

  } else if (!right->rows->first) {
    // fprintf(stderr, "!right->rows true\n");
    // Right empty, left not empty. Do all values in left == right->default_val?
    if (!list::eqeq_value<LDType,RDType>(left->rows, reinterpret_cast<RDType*>(right->default_val), left->rank-1, num_checked)) {
    	return false;
    	
    } else if (num_checked < max_elements) {
   		// If the matrix isn't full, we also need to compare default values.
    	return *reinterpret_cast<LDType*>(left->default_val) == *reinterpret_cast<RDType*>(right->default_val);
    }

  } else {
    // fprintf(stderr, "both matrices have entries\n");
    // Hardest case. Compare lists node by node. Let's make it simpler by requiring that both have the same default value
    if (!list::eqeq<LDType,RDType>(left->rows, right->rows, reinterpret_cast<LDType*>(left->default_val), reinterpret_cast<RDType*>(right->default_val), left->rank-1, num_checked)) {
    	return false;
    	
    } else if (num_checked < max_elements) {
      return *reinterpret_cast<LDType*>(left->default_val) == *reinterpret_cast<RDType*>(right->default_val);
    }
  }

  return true;
}

/*
 * List storage element-wise addition.
 */
template <typename LDType, typename RDType>
static void* ew_add(LIST* dest, const LIST* left, const void* l_default, const LIST* right, const void* r_default, const size_t* shape, size_t rank) {
	
	/*
	 * Allocate space for, and calculate, the default value for the destination
	 * matrix.
	 */
	LDType* d_default_mem = ALLOC(LDType);
	*d_default_mem = *reinterpret_cast<const LDType*>(l_default) + *reinterpret_cast<const RDType*>(r_default);
	
	// Now that setup is done call the actual elementwise multiplication function.
	ew_add_prime<LDType, RDType>(dest, *reinterpret_cast<const LDType*>(d_default_mem),
		left, *reinterpret_cast<const LDType*>(l_default), right, *reinterpret_cast<const RDType*>(r_default), shape, rank - 1, 0);
	
	// Return a pointer to the destination matrix's default value.
	return d_default_mem;
}

/*
 * List storage element-wise subtraction.
 */
template <typename LDType, typename RDType>
static void* ew_subtract(LIST* dest, const LIST* left, const void* l_default, const LIST* right, const void* r_default, const size_t* shape, size_t rank) {
	
	/*
	 * Allocate space for, and calculate, the default value for the destination
	 * matrix.
	 */
	LDType* d_default_mem = ALLOC(LDType);
	*d_default_mem = *reinterpret_cast<const LDType*>(l_default) - *reinterpret_cast<const RDType*>(r_default);
	
	// Now that setup is done call the actual elementwise multiplication function.
	ew_subtract_prime<LDType, RDType>(dest, *reinterpret_cast<const LDType*>(d_default_mem),
		left, *reinterpret_cast<const LDType*>(l_default), right, *reinterpret_cast<const RDType*>(r_default), shape, rank - 1, 0);
	
	// Return a pointer to the destination matrix's default value.
	return d_default_mem;
}

/*
 * List storage element-wise multiplication.
 */
template <typename LDType, typename RDType>
static void* ew_multiply(LIST* dest, const LIST* left, const void* l_default, const LIST* right, const void* r_default, const size_t* shape, size_t rank) {
	
	/*
	 * Allocate space for, and calculate, the default value for the destination
	 * matrix.
	 */
	LDType* d_default_mem = ALLOC(LDType);
	*d_default_mem = *reinterpret_cast<const LDType*>(l_default) * *reinterpret_cast<const RDType*>(r_default);
	
	// Now that setup is done call the actual elementwise multiplication function.
	ew_multiply_prime<LDType, RDType>(dest, *reinterpret_cast<const LDType*>(d_default_mem),
		left, *reinterpret_cast<const LDType*>(l_default), right, *reinterpret_cast<const RDType*>(r_default), shape, rank - 1, 0);
	
	// Return a pointer to the destination matrix's default value.
	return d_default_mem;
}

/*
 * List storage element-wise division.
 */
template <typename LDType, typename RDType>
static void* ew_divide(LIST* dest, const LIST* left, const void* l_default, const LIST* right, const void* r_default, const size_t* shape, size_t rank) {
	
	/*
	 * Allocate space for, and calculate, the default value for the destination
	 * matrix.
	 */
	LDType* d_default_mem = ALLOC(LDType);
	*d_default_mem = *reinterpret_cast<const LDType*>(l_default) / *reinterpret_cast<const RDType*>(r_default);
	
	// Now that setup is done call the actual elementwise multiplication function.
	ew_divide_prime<LDType, RDType>(dest, *reinterpret_cast<const LDType*>(d_default_mem),
		left, *reinterpret_cast<const LDType*>(l_default), right, *reinterpret_cast<const RDType*>(r_default), shape, rank - 1, 0);
	
	// Return a pointer to the destination matrix's default value.
	return d_default_mem;
}

/*
 * List storage element-wise addition, recursive helper.
 */
template <typename LDType, typename RDType>
static void ew_add_prime(LIST* dest, LDType d_default, const LIST* left, LDType l_default, const LIST* right, RDType r_default, const size_t* shape, size_t last_level, size_t level) {
	
	static LIST EMPTY_LIST = {NULL};
	
	size_t index;
	
	LDType tmp_result;
	
	LIST* new_level = NULL;
	
	NODE* l_node		= left->first,
			* r_node		= right->first,
			* dest_node	= NULL;
	
	for (index = 0; index < shape[level]; ++index) {
		if (l_node == NULL and r_node == NULL) {
			/*
			 * Both source lists are now empty.  Because the default value of the
			 * destination is already set appropriately we can now return.
			 */
			
			return;
			
		} else {
			// At least one list still has entries.
			
			if (l_node == NULL and (l_default == 0 and d_default == 0)) {
				/* 
				 * The left hand list has run out of elements.  We don't need to add new
				 * values to the destination if l_default and d_default are both 0.
				 */
				
				return;
			
			} else if (r_node == NULL and (r_default == 0 and d_default == 0)) {
				/*
				 * The right hand list has run out of elements.  We don't need to add new
				 * values to the destination if r_default and d_default are both 0.
				 */
				
				return;
			}
			
			// We need to continue processing the lists.
			
			if (l_node == NULL and r_node->key == index) {
				/*
				 * One source list is empty, but the index has caught up to the key of
				 * the other list.
				 */
				
				if (level == last_level) {
					tmp_result = l_default + *reinterpret_cast<RDType*>(r_node->val);
					
					if (tmp_result != d_default) {
						dest_node = list::insert_helper<LDType>(dest, dest_node, index, tmp_result);
					}
					
				} else {
					new_level = list::create();
					dest_node = list::insert_helper<LIST*>(dest, dest_node, index, new_level);
				
					ew_add_prime<LDType, RDType>(new_level, d_default,
						&EMPTY_LIST, l_default,
						reinterpret_cast<LIST*>(r_node->val), r_default,
						shape, last_level, level + 1);
				}
				
				r_node = r_node->next;
				
			} else if (r_node == NULL and l_node->key == index) {
				/*
				 * One source list is empty, but the index has caught up to the key of
				 * the other list.
				 */
				
				if (level == last_level) {
					tmp_result = *reinterpret_cast<LDType*>(l_node->val) + r_default;
					
					if (tmp_result != d_default) {
						dest_node = list::insert_helper<LDType>(dest, dest_node, index, tmp_result);
					}
					
				} else {
					new_level = list::create();
					dest_node = list::insert_helper<LIST*>(dest, dest_node, index, new_level);
				
					ew_add_prime<LDType, RDType>(new_level, d_default,
						reinterpret_cast<LIST*>(r_node->val), l_default,
						&EMPTY_LIST, r_default,
						shape, last_level, level + 1);
				}
				
				l_node = l_node->next;
				
			} else if (l_node != NULL and r_node != NULL and index == NM_MIN(l_node->key, r_node->key)) {
				/*
				 * Neither list is empty and our index has caught up to one of the
				 * source lists.
				 */
				
				if (l_node->key == r_node->key) {
					
					if (level == last_level) {
						tmp_result = *reinterpret_cast<LDType*>(l_node->val) + *reinterpret_cast<RDType*>(r_node->val);
						
						if (tmp_result != d_default) {
							dest_node = list::insert_helper<LDType>(dest, dest_node, index, tmp_result);
						}
						
					} else {
						new_level = list::create();
						dest_node = list::insert_helper<LIST*>(dest, dest_node, index, new_level);
					
						ew_add_prime<LDType, RDType>(new_level, d_default,
							reinterpret_cast<LIST*>(l_node->val), l_default,
							reinterpret_cast<LIST*>(r_node->val), r_default,
							shape, last_level, level + 1);
					}
				
					l_node = l_node->next;
					r_node = r_node->next;
			
				} else if (l_node->key < r_node->key) {
					// Advance the left node knowing that the default value is OK.
			
					l_node = l_node->next;
					 
				} else /* if (l_node->key > r_node->key) */ {
					// Advance the right node knowing that the default value is OK.
			
					r_node = r_node->next;
				}
				
			} else {
				/*
				 * Our index needs to catch up but the default value is OK.  This
				 * conditional is here only for documentation and should be optimized
				 * out.
				 */
			}
		}
	}
}

/*
 * List storage element-wise subtraction, recursive helper.
 */
template <typename LDType, typename RDType>
static void ew_subtract_prime(LIST* dest, LDType d_default, const LIST* left, LDType l_default, const LIST* right, RDType r_default, const size_t* shape, size_t last_level, size_t level) {
	
	static LIST EMPTY_LIST = {NULL};
	
	size_t index;
	
	LDType tmp_result;
	
	LIST* new_level = NULL;
	
	NODE* l_node		= left->first,
			* r_node		= right->first,
			* dest_node	= NULL;
	
	for (index = 0; index < shape[level]; ++index) {
		if (l_node == NULL and r_node == NULL) {
			/*
			 * Both source lists are now empty.  Because the default value of the
			 * destination is already set appropriately we can now return.
			 */
			
			return;
			
		} else {
			// At least one list still has entries.
			
			if (l_node == NULL and (l_default == 0 and d_default == 0)) {
				/* 
				 * The left hand list has run out of elements.  We don't need to add new
				 * values to the destination if l_default and d_default are both 0.
				 */
				
				return;
			
			} else if (r_node == NULL and (r_default == 0 and d_default == 0)) {
				/*
				 * The right hand list has run out of elements.  We don't need to add new
				 * values to the destination if r_default and d_default are both 0.
				 */
				
				return;
			}
			
			// We need to continue processing the lists.
			
			if (l_node == NULL and r_node->key == index) {
				/*
				 * One source list is empty, but the index has caught up to the key of
				 * the other list.
				 */
				
				if (level == last_level) {
					tmp_result = l_default - *reinterpret_cast<RDType*>(r_node->val);
					
					if (tmp_result != d_default) {
						dest_node = list::insert_helper<LDType>(dest, dest_node, index, tmp_result);
					}
					
				} else {
					new_level = list::create();
					dest_node = list::insert_helper<LIST*>(dest, dest_node, index, new_level);
				
					ew_subtract_prime<LDType, RDType>(new_level, d_default,
						&EMPTY_LIST, l_default,
						reinterpret_cast<LIST*>(r_node->val), r_default,
						shape, last_level, level + 1);
				}
				
				r_node = r_node->next;
				
			} else if (r_node == NULL and l_node->key == index) {
				/*
				 * One source list is empty, but the index has caught up to the key of
				 * the other list.
				 */
				
				if (level == last_level) {
					tmp_result = *reinterpret_cast<LDType*>(l_node->val) - r_default;
					
					if (tmp_result != d_default) {
						dest_node = list::insert_helper<LDType>(dest, dest_node, index, tmp_result);
					}
					
				} else {
					new_level = list::create();
					dest_node = list::insert_helper<LIST*>(dest, dest_node, index, new_level);
				
					ew_subtract_prime<LDType, RDType>(new_level, d_default,
						reinterpret_cast<LIST*>(r_node->val), l_default,
						&EMPTY_LIST, r_default,
						shape, last_level, level + 1);
				}
				
				l_node = l_node->next;
				
			} else if (l_node != NULL and r_node != NULL and index == NM_MIN(l_node->key, r_node->key)) {
				/*
				 * Neither list is empty and our index has caught up to one of the
				 * source lists.
				 */
				
				if (l_node->key == r_node->key) {
					
					if (level == last_level) {
						tmp_result = *reinterpret_cast<LDType*>(l_node->val) - *reinterpret_cast<RDType*>(r_node->val);
						
						if (tmp_result != d_default) {
							dest_node = list::insert_helper<LDType>(dest, dest_node, index, tmp_result);
						}
						
					} else {
						new_level = list::create();
						dest_node = list::insert_helper<LIST*>(dest, dest_node, index, new_level);
					
						ew_subtract_prime<LDType, RDType>(new_level, d_default,
							reinterpret_cast<LIST*>(l_node->val), l_default,
							reinterpret_cast<LIST*>(r_node->val), r_default,
							shape, last_level, level + 1);
					}
				
					l_node = l_node->next;
					r_node = r_node->next;
			
				} else if (l_node->key < r_node->key) {
					// Advance the left node knowing that the default value is OK.
			
					l_node = l_node->next;
					 
				} else /* if (l_node->key > r_node->key) */ {
					// Advance the right node knowing that the default value is OK.
			
					r_node = r_node->next;
				}
				
			} else {
				/*
				 * Our index needs to catch up but the default value is OK.  This
				 * conditional is here only for documentation and should be optimized
				 * out.
				 */
			}
		}
	}
}

/*
 * List storage element-wise multiplication, recursive helper.
 */
template <typename LDType, typename RDType>
static void ew_multiply_prime(LIST* dest, LDType d_default, const LIST* left, LDType l_default, const LIST* right, RDType r_default, const size_t* shape, size_t last_level, size_t level) {
	
	static LIST EMPTY_LIST = {NULL};
	
	size_t index;
	
	LDType tmp_result;
	
	LIST* new_level = NULL;
	
	NODE* l_node		= left->first,
			* r_node		= right->first,
			* dest_node	= NULL;
	
	for (index = 0; index < shape[level]; ++index) {
		if (l_node == NULL and r_node == NULL) {
			/*
			 * Both source lists are now empty.  Because the default value of the
			 * destination is already set appropriately we can now return.
			 */
			
			return;
			
		} else {
			// At least one list still has entries.
			
			if (l_node == NULL and (l_default == 0 and d_default == 0)) {
				/* 
				 * The left hand list has run out of elements.  We don't need to add new
				 * values to the destination if l_default and d_default are both 0.
				 */
				
				return;
			
			} else if (r_node == NULL and (r_default == 0 and d_default == 0)) {
				/*
				 * The right hand list has run out of elements.  We don't need to add new
				 * values to the destination if r_default and d_default are both 0.
				 */
				
				return;
			}
			
			// We need to continue processing the lists.
			
			if (l_node == NULL and r_node->key == index) {
				/*
				 * One source list is empty, but the index has caught up to the key of
				 * the other list.
				 */
				
				if (level == last_level) {
					tmp_result = l_default * *reinterpret_cast<RDType*>(r_node->val);
					
					if (tmp_result != d_default) {
						dest_node = list::insert_helper<LDType>(dest, dest_node, index, tmp_result);
					}
					
				} else {
					new_level = list::create();
					dest_node = list::insert_helper<LIST*>(dest, dest_node, index, new_level);
				
					ew_multiply_prime<LDType, RDType>(new_level, d_default,
						&EMPTY_LIST, l_default,
						reinterpret_cast<LIST*>(r_node->val), r_default,
						shape, last_level, level + 1);
				}
				
				r_node = r_node->next;
				
			} else if (r_node == NULL and l_node->key == index) {
				/*
				 * One source list is empty, but the index has caught up to the key of
				 * the other list.
				 */
				
				if (level == last_level) {
					tmp_result = *reinterpret_cast<LDType*>(l_node->val) * r_default;
					
					if (tmp_result != d_default) {
						dest_node = list::insert_helper<LDType>(dest, dest_node, index, tmp_result);
					}
					
				} else {
					new_level = list::create();
					dest_node = list::insert_helper<LIST*>(dest, dest_node, index, new_level);
				
					ew_multiply_prime<LDType, RDType>(new_level, d_default,
						reinterpret_cast<LIST*>(r_node->val), l_default,
						&EMPTY_LIST, r_default,
						shape, last_level, level + 1);
				}
				
				l_node = l_node->next;
				
			} else if (l_node != NULL and r_node != NULL and index == NM_MIN(l_node->key, r_node->key)) {
				/*
				 * Neither list is empty and our index has caught up to one of the
				 * source lists.
				 */
				
				if (l_node->key == r_node->key) {
					
					if (level == last_level) {
						tmp_result = *reinterpret_cast<LDType*>(l_node->val) * *reinterpret_cast<RDType*>(r_node->val);
						
						if (tmp_result != d_default) {
							dest_node = list::insert_helper<LDType>(dest, dest_node, index, tmp_result);
						}
						
					} else {
						new_level = list::create();
						dest_node = list::insert_helper<LIST*>(dest, dest_node, index, new_level);
					
						ew_multiply_prime<LDType, RDType>(new_level, d_default,
							reinterpret_cast<LIST*>(l_node->val), l_default,
							reinterpret_cast<LIST*>(r_node->val), r_default,
							shape, last_level, level + 1);
					}
				
					l_node = l_node->next;
					r_node = r_node->next;
			
				} else if (l_node->key < r_node->key) {
					// Advance the left node knowing that the default value is OK.
			
					l_node = l_node->next;
					 
				} else /* if (l_node->key > r_node->key) */ {
					// Advance the right node knowing that the default value is OK.
			
					r_node = r_node->next;
				}
				
			} else {
				/*
				 * Our index needs to catch up but the default value is OK.  This
				 * conditional is here only for documentation and should be optimized
				 * out.
				 */
			}
		}
	}
}

/*
 * List storage element-wise multiplication, recursive helper.
 */
template <typename LDType, typename RDType>
static void ew_divide_prime(LIST* dest, LDType d_default, const LIST* left, LDType l_default, const LIST* right, RDType r_default, const size_t* shape, size_t last_level, size_t level) {
	
	static LIST EMPTY_LIST = {NULL};
	
	size_t index;
	
	LDType tmp_result;
	
	LIST* new_level = NULL;
	
	NODE* l_node		= left->first,
			* r_node		= right->first,
			* dest_node	= NULL;
	
	for (index = 0; index < shape[level]; ++index) {
		if (l_node == NULL and r_node == NULL) {
			/*
			 * Both source lists are now empty.  Because the default value of the
			 * destination is already set appropriately we can now return.
			 */
			
			return;
			
		} else {
			// At least one list still has entries.
			
			if (l_node == NULL and (l_default == 0 and d_default == 0)) {
				/* 
				 * The left hand list has run out of elements.  We don't need to add new
				 * values to the destination if l_default and d_default are both 0.
				 */
				
				return;
			
			} else if (r_node == NULL and (r_default == 0 and d_default == 0)) {
				/*
				 * The right hand list has run out of elements.  We don't need to add new
				 * values to the destination if r_default and d_default are both 0.
				 */
				
				return;
			}
			
			// We need to continue processing the lists.
			
			if (l_node == NULL and r_node->key == index) {
				/*
				 * One source list is empty, but the index has caught up to the key of
				 * the other list.
				 */
				
				if (level == last_level) {
					tmp_result = l_default / *reinterpret_cast<RDType*>(r_node->val);
					
					if (tmp_result != d_default) {
						dest_node = list::insert_helper<LDType>(dest, dest_node, index, tmp_result);
					}
					
				} else {
					new_level = list::create();
					dest_node = list::insert_helper<LIST*>(dest, dest_node, index, new_level);
				
					ew_divide_prime<LDType, RDType>(new_level, d_default,
						&EMPTY_LIST, l_default,
						reinterpret_cast<LIST*>(r_node->val), r_default,
						shape, last_level, level + 1);
				}
				
				r_node = r_node->next;
				
			} else if (r_node == NULL and l_node->key == index) {
				/*
				 * One source list is empty, but the index has caught up to the key of
				 * the other list.
				 */
				
				if (level == last_level) {
					tmp_result = *reinterpret_cast<LDType*>(l_node->val) / r_default;
					
					if (tmp_result != d_default) {
						dest_node = list::insert_helper<LDType>(dest, dest_node, index, tmp_result);
					}
					
				} else {
					new_level = list::create();
					dest_node = list::insert_helper<LIST*>(dest, dest_node, index, new_level);
				
					ew_divide_prime<LDType, RDType>(new_level, d_default,
						reinterpret_cast<LIST*>(r_node->val), l_default,
						&EMPTY_LIST, r_default,
						shape, last_level, level + 1);
				}
				
				l_node = l_node->next;
				
			} else if (l_node != NULL and r_node != NULL and index == NM_MIN(l_node->key, r_node->key)) {
				/*
				 * Neither list is empty and our index has caught up to one of the
				 * source lists.
				 */
				
				if (l_node->key == r_node->key) {
					
					if (level == last_level) {
						tmp_result = *reinterpret_cast<LDType*>(l_node->val) / *reinterpret_cast<RDType*>(r_node->val);
						
						if (tmp_result != d_default) {
							dest_node = list::insert_helper<LDType>(dest, dest_node, index, tmp_result);
						}
						
					} else {
						new_level = list::create();
						dest_node = list::insert_helper<LIST*>(dest, dest_node, index, new_level);
					
						ew_divide_prime<LDType, RDType>(new_level, d_default,
							reinterpret_cast<LIST*>(l_node->val), l_default,
							reinterpret_cast<LIST*>(r_node->val), r_default,
							shape, last_level, level + 1);
					}
				
					l_node = l_node->next;
					r_node = r_node->next;
			
				} else if (l_node->key < r_node->key) {
					// Advance the left node knowing that the default value is OK.
			
					l_node = l_node->next;
					 
				} else /* if (l_node->key > r_node->key) */ {
					// Advance the right node knowing that the default value is OK.
			
					r_node = r_node->next;
				}
				
			} else {
				/*
				 * Our index needs to catch up but the default value is OK.  This
				 * conditional is here only for documentation and should be optimized
				 * out.
				 */
			}
		}
	}
}

}} // end of namespace nm::list_storage
