#ifndef LIST_C
# define LIST_C

#include <ruby.h>

#include "nmatrix.h"


/* Get the contents of some set of coordinates. Note: Does not make a copy! Don't free! */
void* list_storage_get(LIST_STORAGE* s, size_t* coords) {
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

/* Creates a LIL matrix of n dimensions */
LIST_STORAGE* create_list_storage(size_t elem_size, size_t* shape, size_t rank, void* init_val) {
  LIST_STORAGE* s;
  if (!(s = malloc(sizeof(LIST_STORAGE)))) return NULL;
  s->rows  = create_list();
  s->shape = shape;
  s->rank  = rank;
  if (!(s->default_val = malloc(sizeof(elem_size)))) {
    delete_list( s->rows, s->rank - 1 );
    free(s);
    return NULL;
  }
  memcpy(s->default_val, init_val, elem_size);

  return s;
}


void delete_list_storage(LIST_STORAGE* s) {
  delete_list( s->rows, s->rank - 1 );
  free(s->default_val);
  free(s);
}


/* Creates an empty linked list */
LIST* create_list() {
  LIST* list;
  if (!(list = malloc(sizeof(LIST)))) return NULL;
  list->first = NULL;
  return list;
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

    if (recursions == 0)
      free(curr->val); // TODO: Check that this doesn't create a memory leak, since val is a void*.
    else
      delete_list(curr->val, recursions-1);

    free(curr);
    curr = next;
  }
  free(list);
}


void list_print_list(LIST* list) {
  NODE* curr = list->first;

  printf("[ ");
  while (curr != NULL) {
    printf("%d->", (int)(curr->key));
    list_print_int(curr->val);
    printf(" ");
    curr = curr->next;
  }

  printf(" ]");
}


void list_print_int(LIST* list) {
  NODE* curr = list->first;

  printf("[ ");
  while (curr != NULL) {
    printf("%d->%d ", (int)(curr->key), *((int *)(curr->val)));
    curr = curr->next;
  }
  printf(" ]");
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

/* int main() {
    int *v1, *v2, *v3, *v4, *v5;
    LIST_STORAGE* s;
    size_t c0[] = {1,3};
    size_t c1[] = {1,2};
    size_t c2[] = {1,4};
    size_t c3[] = {0,3};
    LIST* l = create_list();

    printf("Test\n");

    v1 = malloc(sizeof(int));
    *v1 = 5;

    v2 = malloc(sizeof(int));
    *v2 = 4;

    v3 = malloc(sizeof(int));
    *v3 = -2;

    v4 = malloc(sizeof(int));
    *v4 = -50;

    v5 = malloc(sizeof(int));
    *v5 = 10;

    if (list_insert(l, true, 2, v1))
        printf("First insertion successful\n");
    list_print_int(l);

    list_insert(l, true, 2, v2);
    list_print_int(l);
    list_insert(l, true, 0, v3);
    list_print_int(l);
    list_insert(l, true, 1, v4);
    list_print_int(l);
    list_insert(l, true, 3, v5);
    list_print_int(l);

    NODE* n = list_find_preceding_from(l->first, 2);
    if (n) printf("Found key %d, value=%d\n", (int)(n->key), *((int*)(n->val)));


    v1 = list_remove(l, 2);
    if (v1) printf("Successfully removed value: %d\n", *v1);

    v2 = list_remove(l, 0);
    if (v2) printf("Successfully removed value: %d\n", *v2);

    v3 = list_remove(l, 3);
    if (v3) printf("Successfully removed value: %d\n", *v3);

    v4 = list_remove(l, 1);
    if (v4) printf("Successfully removed value: %d\n", *v4);

    delete_list(l, 0);

    s = create_list_storage(2);
    list_storage_insert(s, c0, v4);
    list_print_list(s->rows); printf("\n");
    list_storage_insert(s, c1, v3);
    list_print_list(s->rows); printf("\n");
    list_storage_insert(s, c2, v2);
    list_print_list(s->rows); printf("\n");
    list_storage_insert(s, c3, v1);
    list_print_list(s->rows); printf("\n");

    return 0;
}
*/

#endif