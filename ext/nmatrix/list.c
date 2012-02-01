#ifndef LIST_C
# define LIST_C

#include <ruby.h>

#include "nmatrix.h"


/* Get the contents of some set of coordinates. Note: Does not make a copy! Don't free! */
void* list_storage_get(LIST_STORAGE* s, size_t* coords) {
  //LIST_STORAGE* s = (LIST_STORAGE*)(t);
  size_t r;
  NODE*  n;
  LIST*  l = s->rows;

  for (r = s->rank; r > 1; --r) {
    n = list_find(l, coords[s->rank - r]);
    if (n)  l = n->val;
    else return s->default_val;
  }

  n = list_find(l, coords[s->rank - r]);
  if (n) return n->val;
  else   return s->default_val;
}

/// TODO: Speed up removal.
void* list_storage_remove(LIST_STORAGE* s, size_t* coords) {
  int r;
  NODE  *n = NULL;
  LIST*  l = s->rows;
  void*  rm = NULL;

  // keep track of where we are in the traversals
  NODE** stack = malloc((s->rank - 1) * sizeof(NODE*));

  for (r = (int)(s->rank); r > 1; --r) {
    n = list_find(l, coords[s->rank - r]); // does this row exist in the matrix?

    if (!n) { // not found
      free(stack);
      return NULL;
    } else { // found
      stack[s->rank - r]    = n;
      l                     = n->val;
    }
  }

  rm = list_remove(l, coords[s->rank - r]);

  // if we removed something, we may now need to remove parent lists
  if (rm) {
    for (r = (int)(s->rank) - 2; r >= 0; --r) { // walk back down the stack
      if (((LIST*)(stack[r]->val))->first == NULL)
        free(list_remove(stack[r]->val, coords[r]));
      else
        break; // no need to continue unless we just deleted one.
    }
  }

  free(stack);
  return rm;
}


// TODO: Allow this function to accept an entire row and not just one value -- for slicing
void* list_storage_insert(LIST_STORAGE* s, size_t* coords, void* val) {
  // Pretend ranks = 2
  // Then coords is going to be size 2
  // So we need to find out if some key already exists
  size_t r;
  NODE*  n;
  LIST*  l = s->rows;

  // drill down into the structure
  for (r = s->rank; r > 1; --r) {
    n = list_insert(l, false, coords[s->rank - r], create_list());
    l = n->val;
  }

  n = list_insert(l, true, coords[s->rank - r], val);
  return n->val;
}

// Creates a list-of-lists(-of-lists-of-lists-etc) storage framework for a matrix.
//
// Note: The pointers you pass in for shape and init_val become property of our new
// storage. You don't need to free them, and you shouldn't re-use them.
LIST_STORAGE* create_list_storage(int8_t dtype, size_t* shape, size_t rank, void* init_val) {
  LIST_STORAGE* s;

  if (!(s = malloc(sizeof(LIST_STORAGE)))) return NULL;

  //fprintf(stderr, "Creating list storage at %p\n", s);

  s->rank  = rank;
  s->shape = shape;
  s->dtype = dtype;

  if (!(s->rows  = create_list())) {
    free(s);
    return NULL;
  }

  s->default_val = init_val;

  return s;
}


LIST_STORAGE* copy_list_storage(LIST_STORAGE* rhs) {
  LIST_STORAGE* lhs;
  size_t* shape;
  void* default_val = malloc(nm_sizeof[rhs->dtype]);

  //fprintf(stderr, "copy_list_storage\n");

  // allocate and copy shape
  shape = malloc( sizeof(size_t) * rhs->rank );
  memcpy(shape, rhs->shape, rhs->rank * sizeof(size_t));
  memcpy(default_val, rhs->default_val, nm_sizeof[rhs->dtype]);

  lhs = create_list_storage(rhs->dtype, shape, rhs->rank, default_val);

  if (lhs) {
    lhs->rows = create_list();
    copy_list_contents(lhs->rows, rhs->rows, rhs->dtype, rhs->rank - 1);
  } else free(shape);

  return lhs;
}


void delete_list_storage(LIST_STORAGE* s) {
  if (s) {
    //fprintf(stderr, "* Deleting list storage rows at %p\n", s->rows);
    delete_list( s->rows, s->rank - 1 );

    //fprintf(stderr, "  Deleting list storage shape at %p\n", s->shape);
    free(s->shape);
    //fprintf(stderr, "  Deleting list storage default_val at %p\n", s->default_val);
    free(s->default_val);
    //fprintf(stderr, "  Deleting list storage at %p\n", s);
    free(s);
  }
}


/* Creates an empty linked list */
LIST* create_list() {
  LIST* list;
  if (!(list = malloc(sizeof(LIST)))) return NULL;

  //fprintf(stderr, "    create_list LIST: %p\n", list);

  list->first = NULL;
  return list;
}


void copy_list_contents(LIST* lhs, LIST* rhs, int8_t dtype, size_t recursions) {
  NODE *lcurr = NULL, *rcurr = rhs->first;

  if (rhs->first) {
    // copy head node
    rcurr = rhs->first;
    lcurr = lhs->first = malloc(sizeof(NODE));

    while (rcurr != NULL) {
      lcurr->key = rcurr->key;

      if (recursions == 0) { // contents is some kind of value
        lcurr->val = malloc(nm_sizeof[dtype]);
        //fprintf(stderr, "    create_val: %p\n", lcurr->val);

        memcpy(lcurr->val, rcurr->val, nm_sizeof[dtype]);

      } else { // contents is a list
        lcurr->val = malloc(sizeof(LIST));
        //fprintf(stderr, "    create_list: %p\n", lcurr->val);

        copy_list_contents(lcurr->val, rcurr->val, nm_sizeof[dtype], recursions-1);
      }
      if (rcurr->next) lcurr->next = malloc(sizeof(NODE));
      else             lcurr->next = NULL;

      lcurr = lcurr->next;
      rcurr = rcurr->next;
    }
  } else {
    lhs->first = NULL;
  }
}


/* Deletes the linked list and all of its contents. If you want to delete a list inside of a list,
 * set recursions to 1. For lists inside of lists inside of the list, set it to 2; and so on.
 * Setting it to 0 is for no recursions.
 */
void delete_list(LIST* list, size_t recursions) {
  NODE* next;
  NODE* curr = list->first;

  while (curr != NULL) {
    next = curr->next;

    if (recursions == 0) {
      //fprintf(stderr, "    free_val: %p\n", curr->val);
      free(curr->val);
    } else {
      //fprintf(stderr, "    free_list: %p\n", list);
      delete_list(curr->val, recursions-1);
    }

    free(curr);
    curr = next;
  }
  //fprintf(stderr, "    free_list: %p\n", list);
  free(list);
}


/* Find some element in the list and return the node ptr for that key. */
NODE* list_find(LIST* list, size_t key) {
  NODE* f;
  if (!list->first) return NULL; // empty list -- does not exist

  // see if we can find it.
  f = list_find_nearest_from(list->first, key);
  if (!f || f->key == key) return f;
  return NULL;
}

/* Finds a node or the one immediately preceding it if it doesn't exist */
NODE* list_find_nearest_from(NODE* prev, size_t key) {
  NODE* f;

  if (prev && prev->key == key) return prev;

  f = list_find_preceding_from(prev, key);

  if (!f->next) return f;
  else if (key == f->next->key) return f->next;
  else return prev;
}

/* Finds the node that should go before whatever key we request, whether or not that key is present */
NODE* list_find_preceding_from(NODE* prev, size_t key) {
  NODE* curr = prev->next;

  if (!curr || key <= curr->key) return prev;
  return list_find_preceding_from(curr, key);
}


/* Finds the node or, if not present, the node that it should follow.
 * NULL indicates no preceding node. */
NODE* list_find_nearest(LIST* list, size_t key) {
  return list_find_nearest_from(list->first, key);
}


/* Returns the value pointer (not the node) for some key. Note that it doesn't free the memory
 * for the value stored in the node -- that pointer gets returned! Only the node is destroyed.
 */
void* list_remove(LIST* list, size_t key) {
  NODE *f, *rm;
  void* val;

  if (!list->first || list->first->key > key) return NULL; // empty list or def. not present

  if (list->first->key == key) {
    val = list->first->val;
    rm  = list->first;
    list->first = rm->next;
    free(rm);
    return val;
  }

  f = list_find_preceding_from(list->first, key);
  if (!f || !f->next) return NULL; // not found, end of list

  if (f->next->key == key) {
    // remove the node
    rm      = f->next;
    f->next = rm->next;

    // get the value and free the memory for the node
    val = rm->val;
    free(rm);
    return val;
  }

  return NULL; // not found, middle of list
}


NODE* list_insert_after(NODE* node, size_t key, void* val) {
  NODE* ins;

  if (!(ins = malloc(sizeof(NODE)))) return NULL;

  // insert 'ins' between 'node' and 'node->next'
  ins->next  = node->next;
  node->next = ins;

  // initialize our new node
  ins->key  = key;
  ins->val  = val;

  return ins;
}



/* Given a list and a key/value-ptr pair, create a node (and return that node).
 * If NULL is returned, it means insertion failed.
 * If the key already exists in the list, replace tells it to delete the old value
 * and put in your new one. !replace means delete the new value.
 */
NODE* list_insert(LIST* list, bool replace, size_t key, void* val) {
  NODE *ins;

  if (list->first == NULL) {                        // List is empty
    if (!(ins = malloc(sizeof(NODE)))) return NULL;
    ins->next             = NULL;
    ins->val              = val;
    ins->key              = key;
    list->first           = ins;
    return ins;

  } else if (key < list->first->key) {              // Goes at the beginning of the list
    if (!(ins = malloc(sizeof(NODE)))) return NULL;
    ins->next             = list->first;
    ins->val              = val;
    ins->key              = key;
    list->first           = ins;
    return ins;
  }

  // Goes somewhere else in the list.
  ins = list_find_nearest_from(list->first, key);

  if (ins->key == key) {
    // key already exists
    if (replace) {
      free(ins->val);
      ins->val = val;
    } else free(val);
    return ins;

  } else return list_insert_after(ins, key, val);

}

#endif