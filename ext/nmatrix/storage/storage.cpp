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
// == storage.cpp
//
// Code that is used by or involves more then one storage type.

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

#include "data/data.h"

#include "storage.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

const char* const STYPE_NAMES[NUM_STYPES] = {
	"dense",
	"list",
	"yale"
};



/*
 * Forward Declarations
 */

// FIXME: Template
//static void dense_storage_cast_copy_list_contents(void* lhs, const LIST* rhs, void* default_val, dtype_t l_dtype, dtype_t r_dtype,
//	size_t* pos, const size_t* shape, size_t rank, size_t max_elements, size_t recursions);

// FIXME: Template
//static void dense_storage_cast_copy_list_default(void* lhs, void* default_val, dtype_t l_dtype, dtype_t r_dtype,
//	size_t* pos, const size_t* shape, size_t rank, size_t max_elements, size_t recursions);

/*
 * Functions
 */

/////////////////////////
// Templated Functions //
/////////////////////////


#if false

// To be fixed once each storage class is compiling on its own.

/*
 * Convert (by creating a copy) from list storage to dense storage.
 *
 * FIXME: Add templating.
 */
DENSE_STORAGE* dense_storage_from_list(const LIST_STORAGE* rhs, dtype_t l_dtype) {
  DENSE_STORAGE* lhs;
  
  // Position in lhs->elements.
  size_t pos = 0;

  // allocate and copy shape
  size_t* shape = ALLOC_N(size_t, rhs->rank);
  memcpy(shape, rhs->shape, rhs->rank * sizeof(size_t));

  lhs = dense_storage_create(l_dtype, shape, rhs->rank, NULL, 0);

  // recursively copy the contents
  dense_storage_cast_copy_list_contents(lhs->elements, rhs->rows, rhs->default_val, l_dtype, rhs->dtype, &pos, shape, lhs->rank, count_storage_max_elements((STORAGE*)rhs), rhs->rank-1);

  return lhs;
}

/*
 * Create/allocate dense storage, copying into it the contents of a Yale matrix.
 */
template <typename LDType, typename RDType, typename RIType>
DENSE_STORAGE* dense_storage_from_yale_template(const YALE_STORAGE* rhs, dtype_t l_dtype) {

  // Position in rhs->elements.
  RIType* rhs_ija = reinterpret_cast<RIType*>(rhs->ija);
  RDType* rhs_a   = reinterpret_cast<RDType*>(rhs->a);
  
  // Allocate and set shape.
  size_t* shape = ALLOC_N(size_t, rhs->rank);
  memcpy(shape, rhs->shape, rhs->rank * sizeof(size_t));

  DENSE_STORAGE* lhs = dense_storage_create(l_dtype, shape, rhs->rank, NULL, 0);
  LDType* lhs_elements = reinterpret_cast<LDType*>(lhs->elements);

  // Position in yale array
  RIType ija      = rhs_ija[i],
         ija_next = rhs_ija[i+1];

  // Position in dense to write to.
  size_t pos = 0;

  LDType LCAST_ZERO = static_cast<LDType>(rhs_a[rhs->shape[0]]);

  // Walk through rows. For each entry we set in dense, increment pos.
  for (size_t i = rhs->shape[0]; i-- > 0;) {

    if (ija == ija_next) { // Check boundaries of row: is row empty?

			// Write zeros in each column.
			for (size_t j = rhs->shape[1]; j-- > 0;) { // Move to next dense position.

        // Fill in zeros (except for diagonal)
        if (i == j) lhs_elements[pos] = static_cast<LDType>(rhs_a[i]);
				else        lhs_elements[pos] = LCAST_ZERO;

				++pos;
      }

    } else {
      // Row contains entries: write those in each column, interspersed with zeros.
      RIType jj = rhs_ija[ija];

			for (size_t j = rhs->shape[1]; j-- > 0;) {
        if (i == j) {
          lhs_elements[pos] = static_cast<LDType>(rhs_a[i]);

        } else if (j == jj) {
          lhs_elements[pos] = static_cast<LDType>(rhs_a[ija]); // Copy from rhs.

          // Get next.
          ++ija;

          // Increment to next column ID (or go off the end).
          if (ija < ija_next) jj = rhs_ija[ija];
          else               	jj = rhs->shape[1];
          
        } else { // j < jj

          // Insert zero.
          lhs_elements[pos] = LCAST_ZERO;
        }
        
        // Move to next dense position.
        ++pos;
      }
    }
  }

  return lhs;
}

/*
 * Documentation goes here.
 *
 * FIXME: Add templating.
 */
LIST_STORAGE* list_storage_from_dense(const DENSE_STORAGE* rhs, int8_t l_dtype) {
  LIST_STORAGE* lhs;
  size_t pos = 0;
  void* l_default_val = ALLOC_N(char, DTYPE_SIZES[l_dtype]);
  void* r_default_val = ALLOCA_N(char, DTYPE_SIZES[rhs->dtype]); // clean up when finished with this function

  // allocate and copy shape and coords
  size_t *shape = ALLOC_N(size_t, rhs->rank), *coords = ALLOC_N(size_t, rhs->rank);
  memcpy(shape, rhs->shape, rhs->rank * sizeof(size_t));
  memset(coords, 0, rhs->rank * sizeof(size_t));

  // set list default_val to 0
  if (l_dtype == NM_ROBJ) {
  	*(VALUE*)l_default_val = INT2FIX(0);
  	
  } else {
  	memset(l_default_val, 0, DTYPE_SIZES[l_dtype]);
  }

  // need test default value for comparing to elements in dense matrix
  if (rhs->dtype == l_dtype) {
  	r_default_val = l_default_val;
  	
  } else if (rhs->dtype == NM_ROBJ) {
  	*(VALUE*)r_default_val = INT2FIX(0);
  	
  } else {
  	memset(r_default_val, 0, DTYPE_SIZES[rhs->dtype]);
  }

  lhs = create_list_storage(l_dtype, shape, rhs->rank, l_default_val);

  lhs->rows = create_list();
  cast_copy_list_contents_dense(lhs->rows, rhs->elements, r_default_val, l_dtype, rhs->dtype, &pos, coords, rhs->shape, rhs->rank, rhs->rank - 1);

  return lhs;
}

/*
 * Documentation goes here.
 *
 * FIXME: Add templating.
 */
LIST_STORAGE* list_storage_from_yale(const YALE_STORAGE* rhs, int8_t l_dtype) {
  LIST_STORAGE* lhs;
  NODE *last_added, *last_row_added = NULL;
  LIST* curr_row;
  y_size_t ija, ija_next, i, jj;
  bool add_diag;
  void* default_val = ALLOC_N(char, DTYPE_SIZES[l_dtype]);
  void* R_ZERO = (char*)(rhs->a) + rhs->shape[0]*DTYPE_SIZES[rhs->dtype];
  void* insert_val;

  // allocate and copy shape
  size_t *shape = ALLOC_N(size_t, rhs->rank);
  shape[0] = rhs->shape[0]; shape[1] = rhs->shape[1];

  // copy default value from the zero location in the Yale matrix
  SetFuncs[l_dtype][rhs->dtype](1, default_val, 0, R_ZERO, 0);

  lhs = create_list_storage(l_dtype, shape, rhs->rank, default_val);

  if (rhs->rank != 2) {
    rb_raise(nm_eStorageTypeError, "Can only convert matrices of rank 2 from yale.");
  }

  // Walk through rows and columns as if RHS were a dense matrix
  for (i = rhs->shape[0]; i-- > 0;) {

    // Get boundaries of beginning and end of row
    YaleGetIJA(ija, rhs, i);
    YaleGetIJA(ija_next, rhs, i+1);

    // Are we going to need to add a diagonal for this row?
    if (ElemEqEq[rhs->dtype][0]((char*)(rhs->a) + i*DTYPE_SIZES[rhs->dtype], R_ZERO, 1, DTYPE_SIZES[rhs->dtype])) {
    	// zero
    	add_diag = false;
    	
    } else {
    	// nonzero diagonal
    	add_diag = true;
    }
		
    if (ija < ija_next || add_diag) {

      curr_row = create_list();
      last_added = NULL;

      while (ija < ija_next) {
        YaleGetIJA(jj, rhs, ija); // what column number is this?

        // Is there a nonzero diagonal item between the previously added item and the current one?
        if (jj > i && add_diag) {
          // Allocate and copy insertion value
          insert_val = ALLOC_N(char, DTYPE_SIZES[l_dtype]);
          SetFuncs[l_dtype][rhs->dtype](1, insert_val, 0, (char*)(rhs->a) + i*DTYPE_SIZES[rhs->dtype], 0);
					
          // insert the item in the list at the appropriate location
          if (last_added) {
          	last_added = list_insert_after(last_added, i, insert_val);
          	
          } else {
          	last_added = list_insert(curr_row, false, i, insert_val);
          }
					
					// don't add again!
          add_diag = false;
        }

        // now allocate and add the current item
        insert_val = ALLOC_N(char, DTYPE_SIZES[l_dtype]);
        SetFuncs[l_dtype][rhs->dtype](1, insert_val, 0, (char*)(rhs->a) + ija*DTYPE_SIZES[rhs->dtype], 0);

        if (last_added) {
        	last_added = list_insert_after(last_added, jj, insert_val);
        	
        } else {
        	last_added = list_insert(curr_row, false, jj, insert_val);
        }

        ++ija; // move to next entry in Yale matrix
      }

      if (add_diag) {
      	// still haven't added the diagonal.
      	
        insert_val = ALLOC_N(char, DTYPE_SIZES[l_dtype]);
        SetFuncs[l_dtype][rhs->dtype](1, insert_val, 0, (char*)(rhs->a) + i*DTYPE_SIZES[rhs->dtype], 0);

        // insert the item in the list at the appropriate location
        if (last_added) {
        	last_added = list_insert_after(last_added, i, insert_val);
        	
        } else {
        	last_added = list_insert(curr_row, false, i, insert_val);
        }
      }

      // Now add the list at the appropriate location
      if (last_row_added) {
      	last_row_added = list_insert_after(last_row_added, i, curr_row);
      	
      } else {
      	last_row_added = list_insert(lhs->rows, false, i, curr_row);
      }
    }
	
		// end of walk through rows
  }

  return lhs;
}

/*
 * Creation of yale storage from dense storage.
 */
template <typename LDType, typename LIType, typename RDType>
YALE_STORAGE* yale_storage_from_dense_template(const DENSE_STORAGE* rhs, dtype_t l_dtype) {
  LIType pos = 0, ndnz = 0;

  RDType R_ZERO; // need zero for easier comparisons
  if (rhs->dtype == RUBYOBJ)  R_ZERO = INT2FIX(0);
  else                        R_ZERO = 0;

  if (rhs->rank != 2) rb_raise(nm_eStorageTypeError, "can only convert matrices of rank 2 to yale");

  RDType* rhs_elements = reinterpret_cast<RDType*>(rhs->elements);

  // First, count the non-diagonal nonzeros
	for (size_t i = rhs->shape[0]; i-- > 0;) {
		for (size_t j = rhs->shape[1]; j-- > 0;) {
		  if (i != j && rhs_elements[pos] != R_ZERO)	++ndnz;

      // move forward 1 position in dense matrix elements array
      ++pos;
    }
  }

  // Copy shape for yale construction
  size_t* shape = ALLOC_N(size_t, 2);
  shape[0] = rhs->shape[0];
  shape[1] = rhs->shape[1];

  // Create with minimum possible capacity -- just enough to hold all of the entries
  YALE_STORAGE* lhs = yale_storage_create(l_dtype, shape, 2, shape[0] + ndnz + 1);
  LDType* lhs_a     = reinterpret_cast<LDType*>(lhs->a);
  LIType* lhs_ija   = reinterpret_cast<LIType*>(lhs->ija);

  // Set the zero position in the yale matrix
  lhs_a[ shape[0] ] = static_cast<LDType>(R_ZERO);

  // Start just after the zero position.
  LIType ija = lhs->shape[0]+1;
  pos        = 0;

  // Copy contents
  for (LIType i = 0; i < rhs->shape[0]; ++i) {
    // indicate the beginning of a row in the IJA array
    lhs_ija[i] = ija;

    for (LIType j = 0; j < rhs->shape[1]; ++j) {

      if (i == j) { // copy to diagonal
        lhs_a[i] = static_cast<LDType>(rhs_elements[pos]);
      } else if (rhs->elements[pos] != R_ZERO) { // copy nonzero to LU
        lhs_ija[ija] = j; // write column index

        lhs_a[ija] = static_cast<LDType>(rhs->elements[pos]);

        ++ija;
      }
      ++pos;
    }
  }
  lhs_ija[i] = ija; // indicate the end of the last row

  lhs->ndnz = ndnz;

  return lhs;
}

/*
 * Documentation goes here.
 *
 * FIXME: Add templating.
 */
YALE_STORAGE* yale_storage_from_list(const LIST_STORAGE* rhs, int8_t l_dtype) {
  YALE_STORAGE* lhs;
  size_t* shape;
  NODE *i_curr, *j_curr;
  y_size_t ija;
  size_t ndnz = count_list_storage_nd_elements(rhs);

  if (rhs->rank != 2) {
    rb_raise(nm_eStorageTypeError, "can only convert matrices of rank 2 to yale");
  }

  if ((rhs->dtype == NM_ROBJ && *(VALUE*)(rhs->default_val) == INT2FIX(0)) || strncmp(rhs->default_val, "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", DTYPE_SIZES[rhs->dtype])) {
    rb_raise(nm_eStorageTypeError, "list matrix must have default value of 0 to convert to yale");
  }

  // Copy shape for yale construction
  shape = ALLOC_N(size_t, 2);
  shape[0] = rhs->shape[0];
  shape[1] = rhs->shape[1];

  lhs = yale_storage_create(l_dtype, shape, 2, shape[0] + ndnz + 1);
  clear_diagonal_and_zero(lhs); // clear the diagonal and the zero location.

  ija = lhs->shape[0]+1;

  for (i_curr = rhs->rows->first; i_curr; i_curr = i_curr->next) {

    // indicate the beginning of a row in the IJA array
    YaleSetIJA(i_curr->key, lhs, ija);

    for (j_curr = ((LIST*)(i_curr->val))->first; j_curr; j_curr = j_curr->next) {
      if (i_curr->key == j_curr->key) {
        // set diagonal
        SetFuncs[l_dtype][rhs->dtype](1, (char*)(lhs->a) + (i_curr->key)*DTYPE_SIZES[l_dtype], 0, j_curr->val, 0);
        
      } else {
        // set column value
        YaleSetIJA(ija, lhs, j_curr->key); // write column index

        // set cell value
        SetFuncs[l_dtype][rhs->dtype](1, (char*)(lhs->a) + ija*DTYPE_SIZES[l_dtype], 0, j_curr->val, 0);

        ++ija;
      }
    }

    if (!i_curr->next) {
    	// indicate the end of the last row
    	YaleSetIJA(i_curr->key, lhs, ija);
    }
  }

  lhs->ndnz = ndnz;
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
      	cast_copy_value_single((char*)lhs + (*pos)*DTYPE_SIZES[l_dtype], default_val, l_dtype, r_dtype);
    		//fprintf(stderr, "zero\n");
      
      } else {
     		dense_storage_cast_copy_list_default(lhs, default_val, l_dtype, r_dtype, pos, shape, rank, max_elements, recursions-1);
    		//fprintf(stderr, "column of zeros\n");
			}

      ++last_key;
      
    } else {
      //fprintf(stderr, "pos = %u, dim = %u, curr->key = %u, last_key+1 = %d\t", *pos, shape[rank-1-recursions], curr->key, last_key+1);
      
      if (recursions == 0) {
      	cast_copy_value_single((char*)lhs + (*pos)*DTYPE_SIZES[l_dtype], curr->val, l_dtype, r_dtype);
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
    	cast_copy_value_single((char*)lhs + (*pos)*DTYPE_SIZES[l_dtype], default_val, l_dtype, r_dtype);
    	//fprintf(stderr, "zero\n");
    	
    } else {
    	dense_storage_cast_copy_list_default(lhs, default_val, l_dtype, r_dtype, pos, shape, rank, max_elements, recursions-1);
    	//fprintf(stderr, "column of zeros\n");
  	}
  }
  
  --(*pos);
}

#endif

