#ifndef LIST_H
# define LIST_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>
#endif

/* Singly-linked ordered list
 * - holds keys and values
 * - no duplicate keys
 * - keys are ordered
 * - values may be lists themselves
 */

typedef struct l_node { /* Linked list node */
  size_t key;
  void*  val;
  struct l_node * next; // next
} NODE;

typedef struct l_list {
  NODE* first;
} LIST;


typedef struct list_s {
  LIST* rows;
  size_t rank;
  void* default_val;
} LIST_STORAGE;

LIST_STORAGE*   create_list_storage(size_t elem_size, size_t rank, void* init_val);
void            delete_list_storage(LIST_STORAGE* s);

void*           list_storage_get(LIST_STORAGE* s, size_t* coords);
void*           list_storage_insert(LIST_STORAGE* s, size_t* coords, void* val);

LIST*           create_list();
void            delete_list(LIST* list, size_t recursions);
NODE*           list_find(LIST* list, size_t key);
NODE*           list_find_preceding_from(NODE* prev, size_t key);
NODE*           list_find_nearest_from(NODE* prev, size_t key);
void*           list_remove(LIST* list, size_t key);
NODE*           list_insert(LIST* list, bool replace, size_t key, void* val);
NODE*           list_insert_after(NODE* node, size_t key, void* val);
void            list_print_int(LIST* list);


#endif