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
// == sl_list.h
//
// Singly-linked list implementation

#ifndef SL_LIST_H
#define SL_LIST_H

/*
 * Standard Includes
 */

#include <cstdlib>

/*
 * Project Includes
 */

#include "types.h"

#include "data/data.h"

/*
 * Macros
 */

/*
 * Types
 */

/* Singly-linked ordered list
 * - holds keys and values
 * - no duplicate keys
 * - keys are ordered
 * - values may be lists themselves
 */
typedef struct l_node {
  size_t key;
  void*  val;
  struct l_node* next;
} NODE;

typedef struct {
  NODE* first;
} LIST;

/*
 * Data
 */
 

/*
 * Functions
 */
 
////////////////
// Lifecycle //
///////////////

LIST*	list_create(void);
void	list_delete(LIST* list, size_t recursions);
void	list_mark(LIST* list, size_t recursions);

///////////////
// Accessors //
///////////////

NODE* list_insert(LIST* list, bool replace, size_t key, void* val);
NODE* list_insert_after(NODE* node, size_t key, void* val);
void* list_remove(LIST* list, size_t key);

template <typename Type>
inline NODE* list_insert_val_helper(LIST* list, NODE* node, size_t key, Type val) {
	Type* val_mem = ALLOC(Type);
	*val_mem = val;
	
	if (node == NULL) {
		return list_insert(list, false, key, val_mem);
		
	} else {
		return list_insert_after(node, key, val_mem);
	}
}

inline NODE* list_insert_ptr_helper(LIST* list, NODE* node, size_t key, void* ptr) {
	if (node == NULL) {
		return list_insert(list, false, key, ptr);
		
	} else {
		return list_insert_after(node, key, ptr);
	}
}

///////////
// Tests //
///////////

/*
 * Do all values in a list == some value?
 *
 * Note that the template parameters should line up with the first two function parameters. This differs from most
 * other eqeq functions, which use left and right dtypes.
 */
template <typename ListDType, typename ValueDType>
bool list_eqeq_value_template(const LIST* l, const ValueDType* v, size_t recursions, size_t& checked) {
  NODE *next, *curr = l->first;

  while (curr) {
    next = curr->next;

    if (recursions == 0) {
      ++checked;

      if (*reinterpret_cast<ListDType*>(curr->val) != *v) return false;

    } else if (!list_eqeq_value_template<ListDType,ValueDType>((LIST*)curr->val, v, recursions - 1, checked)) {
      return false;
    }

    curr = next;
  }

  return true;
}


/*
 * Are all values in the two lists equal? If one is missing a value, but the
 * other isn't, does the value in the list match the default value?
 */
template <typename LDType, typename RDType>
bool list_eqeq_list_template(const LIST* left, const LIST* right, const LDType* left_val, const RDType* right_val, size_t recursions, size_t& checked) {
  NODE *lnext = NULL, *lcurr = left->first, *rnext = NULL, *rcurr = right->first;

  if (lcurr) lnext = lcurr->next;
  if (rcurr) rnext = rcurr->next;

  while (lcurr && rcurr) {

    if (lcurr->key == rcurr->key) {
    	// MATCHING KEYS

      if (recursions == 0) {
        ++checked;

        if (*reinterpret_cast<LDType*>(lcurr->val) != *reinterpret_cast<RDType*>(rcurr->val)) return false;

      } else if (!list_eqeq_list_template<LDType,RDType>(reinterpret_cast<LIST*>(lcurr->val), (LIST*)rcurr->val, left_val, right_val, recursions - 1, checked)) {
        return false;
      }

      // increment both iterators
      rcurr = rnext;
      if (rcurr) rnext = rcurr->next;
      lcurr = lnext;
      if (lcurr) lnext = lcurr->next;

    } else if (lcurr->key < rcurr->key) {
    	// NON-MATCHING KEYS

      if (recursions == 0) {
        // compare left entry to right default value
        ++checked;

        if (*reinterpret_cast<LDType*>(lcurr->val) != *right_val) return false;

      } else if (!list_eqeq_value_template<LDType,RDType>(reinterpret_cast<LIST*>(lcurr->val), right_val, recursions - 1, checked)) {
        return false;
      }

      // increment left iterator
      lcurr = lnext;
      if (lcurr) lnext = lcurr->next;

    } else {
			// if (rcurr->key < lcurr->key)

      if (recursions == 0) {
        // compare right entry to left default value
        ++checked;

        if (*reinterpret_cast<RDType*>(rcurr->val) != *left_val) return false;

      } else if (!list_eqeq_value_template<RDType,LDType>(reinterpret_cast<LIST*>(rcurr->val), left_val, recursions - 1, checked)) {
        return false;
      }

      // increment right iterator
      rcurr = rnext;
      if (rcurr) rnext = rcurr->next;
    }

  }

  /*
   * One final check, in case we get to the end of one list but not the other
   * one.
   */
  if (lcurr) {
  	// nothing left in right-hand list
  	if (*reinterpret_cast<LDType*>(lcurr->val) != *right_val) return false;

  } else if (rcurr) {
  	// nothing left in left-hand list
  	if (*reinterpret_cast<RDType*>(rcurr->val) != *left_val) return false;

  }

  /*
   * Nothing different between the two lists -- but make sure after this return
   * that you compare the default values themselves, if we haven't visited
   * every value in the two matrices.
   */
  return true;
}

/////////////
// Utility //
/////////////

NODE* list_find(LIST* list, size_t key);
NODE* list_find_preceding_from(NODE* prev, size_t key);
NODE* list_find_nearest(LIST* list, size_t key);
NODE* list_find_nearest_from(NODE* prev, size_t key);

/////////////////////////
// Copying and Casting //
/////////////////////////

template <typename LDType, typename RDType>
void list_cast_copy_contents_template(LIST* lhs, const LIST* rhs, size_t recursions);
void list_cast_copy_contents(LIST* lhs, const LIST* rhs, dtype_t lhs_dtype, dtype_t rhs_dtype, size_t recursions);

#endif // SL_LIST_H
