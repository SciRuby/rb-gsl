// yale.c - "new yale" storage format for 2D matrices
//
// Specifications:
// * dtype and index dtype must necessarily differ
//      * index dtype is defined by whatever unsigned type can store max(rows,cols)
//      * that means vector ija stores only index dtype, but a stores dtype
// * vectors must be able to grow as necessary
//      * maximum size is rows*cols+1

#ifndef YALE_C
# define YALE_C

#include <ruby.h> // mostly for exceptions

#include "nmatrix.h"


extern const char *nm_dtypestring[];



void print_vectors(YALE_STORAGE* s) {
  size_t i;
  fprintf(stderr, "------------------------------\n");
  fprintf(stderr, "dtype:%s\tshape:%dx%d\tndnz:%d\tcapacity:%d\tindex_dtype:%s\n", nm_dtypestring[s->dtype], s->shape[0], s->shape[1], s->ndnz, s->capacity, nm_dtypestring[s->index_dtype]);


  if (s->capacity > 60) rb_raise(rb_eArgError, "overflow in print_vectors; cannot handle that large of a vector");
  // print indices

  fprintf(stderr, "i:\t");
  for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5u ", i);

  fprintf(stderr, "\nija:\t");
  if (YALE_MAX_SIZE(s) < UINT8_MAX)
    for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5u ", *(u_int8_t*)YALE_IJA(s,nm_sizeof[s->index_dtype],i));
  else if (YALE_MAX_SIZE(s) < UINT16_MAX)
    for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5u ", *(u_int16_t*)YALE_IJA(s,nm_sizeof[s->index_dtype],i));
  else if (YALE_MAX_SIZE(s) < UINT32_MAX)
    for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5u ", *(u_int32_t*)YALE_IJA(s,nm_sizeof[s->index_dtype],i));
  else
    for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5u ", *(u_int64_t*)YALE_IJA(s,nm_sizeof[s->index_dtype],i));
  fprintf(stderr, "\n");

  // print values
  fprintf(stderr, "a:\t");
  for (i = 0; i < s->capacity; ++i)
    fprintf(stderr, "%-*.3g ", 5, *(double*)((char*)(s->a) + nm_sizeof[NM_FLOAT64]*i));
  fprintf(stderr, "\n");

  fprintf(stderr, "------------------------------\n");
}


// Determine the index dtype (which will be used for the ija vector). This is determined by matrix shape, not IJA/A vector capacity.
int8_t yale_index_dtype(YALE_STORAGE* s) {
  if (YALE_MAX_SIZE(s) < UINT8_MAX) return NM_INT8;
  else if (YALE_MAX_SIZE(s) < UINT16_MAX) return NM_INT16;
  else if (YALE_MAX_SIZE(s) < UINT32_MAX) return NM_INT32;
  else return NM_INT64;
}



char yale_vector_replace(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n) {
  fprintf(stderr, "yale_vector_replace at pos=%u, n=%u, sizeof_dtype=%u, val=%f\n", pos, n, nm_sizeof[s->dtype], *(double*)val);

  // Now insert the new values
  SetFuncs[s->index_dtype][Y_SIZE_T](n,
                                     pos*nm_sizeof[s->index_dtype] + (char*)(s->ija),
                                     nm_sizeof[s->index_dtype],
                                     j,
                                     sizeof(y_size_t));

  SetFuncs[s->dtype][s->dtype](n,
                               pos*nm_sizeof[s->dtype] + (char*)(s->a),
                               nm_sizeof[s->dtype],
                               val,
                               nm_sizeof[s->dtype]);

  return 'r';
}


char yale_vector_insert_resize(YALE_STORAGE* s, y_size_t current_size, y_size_t pos, y_size_t* j, void* val, y_size_t n) {
  void *new_ija, *new_a;
  // Determine the new capacity for the IJA and A vectors.
  size_t new_capacity = s->capacity * YALE_GROWTH_CONSTANT;
  if (new_capacity > YALE_MAX_SIZE(s)) {
    new_capacity = YALE_MAX_SIZE(s);
    if (current_size + n > YALE_MAX_SIZE(s)) rb_raise(rb_eNoMemError, "insertion size exceeded maximum yale matrix size");
  }
  if (new_capacity < current_size + n) new_capacity = current_size + n;

  fprintf(stderr, "yale_vector_insert_resize: resizing to %u\n", new_capacity);
  fprintf(stderr, "\tpos=%u,j=%u,n=%u\n", pos, j[0], n);


  // Allocate the new vectors.
  new_ija     = malloc( nm_sizeof[s->index_dtype] * new_capacity );
  new_a       = malloc( nm_sizeof[s->dtype]       * new_capacity );

  // Check that allocation succeeded.
  if (!new_ija || !new_a) {
    fprintf(stderr, "Allocation error!\n");
    free(new_a); free(new_ija);
    rb_raise(rb_eNoMemError, "yale sparse vectors are full and there is insufficient memory for growing them");
    return (char)false;
  }

  // Copy all values prior to the insertion site to the new IJA and new A
  SetFuncs[s->index_dtype][s->index_dtype](pos, new_ija, nm_sizeof[s->index_dtype], s->ija, nm_sizeof[s->index_dtype]);
  SetFuncs[s->dtype      ][s->dtype      ](pos, new_a,   nm_sizeof[s->dtype],       s->a,   nm_sizeof[s->dtype]      );

  // Copy all values subsequent to the insertion site to the new IJA and new A, leaving room (size n) for insertion.
  fprintf(stderr, "Inserting at %d, copying %d values to offset %d from offset %d\n", pos, current_size-pos+n-1, pos+n, pos);
  SetFuncs[s->index_dtype][s->index_dtype](current_size-pos+n-1, (char*)new_ija + nm_sizeof[s->index_dtype]*(pos+n), nm_sizeof[s->index_dtype], (char*)(s->ija) + nm_sizeof[s->index_dtype]*pos, nm_sizeof[s->index_dtype]);
  SetFuncs[s->dtype      ][s->dtype      ](current_size-pos+n-1, (char*)new_a   + nm_sizeof[s->dtype]*(pos+n),       nm_sizeof[s->dtype],       (char*)(s->a)   + nm_sizeof[s->dtype      ]*pos, nm_sizeof[s->dtype      ]);

  s->capacity = new_capacity;

  free(s->ija);
  free(s->a);

  s->ija = new_ija;
  s->a   = new_a;

  return 'i';
}


// Insert a value or contiguous values in the ija and a vectors (after ja and diag). Does not free anything; you are responsible!
//
// TODO: Improve this so it can handle non-contiguous element insertions efficiently.
// (For now, we can just sort the elements in the row in question.)
char yale_vector_insert(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n) {
  y_size_t sz, i;

  if (pos < s->shape[0])
    rb_raise(rb_eArgError, "vector insert pos is before beginning of ja; this should not happen");

  YaleGetSize(sz, s);

  if (sz + n > s->capacity) {
    fprintf(stderr, "yvi: %d, %d, %d\n", sz, n, s->capacity);
    yale_vector_insert_resize(s, sz, pos, j, val, n);
  } else {

    fprintf(stderr, "yale_vector_insert: no resize needed\n");
    fprintf(stderr, "\tpos=%u,j=%u,n=%u\n", pos, j[0], n);

    // No resize required:
    // easy (but somewhat slow), just copy elements to the tail, starting at the end, one element at a time.
    // TODO: This can be made slightly more efficient, but only after the tests are written.
    for (i = 0; i < sz - pos; ++i) {
      SetFuncs[s->index_dtype][s->index_dtype](1, (char*)(s->ija) + (sz+n-1-i)*nm_sizeof[s->index_dtype], 0, (char*)(s->ija) + (sz-1-i)*nm_sizeof[s->index_dtype], 0);
      SetFuncs[s->dtype      ][s->dtype      ](1, (char*)(s->a)   + (sz+n-1-i)*nm_sizeof[s->dtype      ], 0, (char*)(s->a)   + (sz-1-i)*nm_sizeof[s->dtype      ], 0);
    }
  }

  // Now insert the new values.
  yale_vector_replace(s, pos, j, val, n);

  return 'i';
}


void delete_yale_storage(YALE_STORAGE* s) {
  if (s) {
    free(s->shape);
    free(s->ija);
    free(s->a);
    free(s);
  }
}


// copy constructor
YALE_STORAGE* copy_yale_storage(YALE_STORAGE* rhs) {
  YALE_STORAGE* lhs = malloc(sizeof(YALE_STORAGE));
  lhs->dtype        = rhs->dtype;
  lhs->rank         = rhs->rank;
  lhs->shape        = malloc( sizeof(size_t) * lhs->rank );
                      memcpy(lhs->shape, rhs->shape, lhs->rank * sizeof(size_t));
  lhs->ndnz         = rhs->ndnz;
  lhs->capacity     = rhs->capacity;
  lhs->index_dtype  = rhs->index_dtype;

  if (!(lhs->ija = malloc(nm_sizeof[lhs->index_dtype] * lhs->capacity))) {
    free(lhs->shape);
    free(lhs);
    return NULL;
  } else {
    if (!(lhs->a = malloc(nm_sizeof[lhs->dtype] * lhs->capacity))) {
      free(lhs->shape);
      free(lhs->ija);
      free(lhs);
      return NULL;
    }
  }

  return lhs;
}



YALE_STORAGE* create_yale_storage(int8_t dtype, size_t* shape, size_t rank, size_t init_capacity) {
  YALE_STORAGE* s;

  if (rank > 2) rb_raise(rb_eNotImpError, "Can only support 2D matrices");

  if (!(s = malloc(sizeof(YALE_STORAGE)))) return NULL;
  s->ndnz        = 0;
  s->dtype       = dtype;
  s->shape       = shape;
  s->rank        = rank;
  s->index_dtype = yale_index_dtype(s);

  if (init_capacity < YALE_MINIMUM(s)) init_capacity = YALE_MINIMUM(s);
  s->capacity    = init_capacity;

  if (!(s->ija = malloc(nm_sizeof[s->index_dtype] * init_capacity))) {
    free(s->shape);
    free(s);
    return NULL;
  } else {
    if (!(s->a = malloc(nm_sizeof[s->dtype] * init_capacity))) {
      free(s->shape);
      free(s->ija);
      free(s);
      return NULL;
    }
  }

  return s;
}


// Empty the matrix
void init_yale_storage(YALE_STORAGE* s) {
  y_size_t i, ia_init = YALE_IA_SIZE(s)+1;

  // clear out IJA vector
  for (i = 0; i < ia_init; ++i)
    SetFuncs[s->index_dtype][Y_SIZE_T](1, (char*)(s->ija) + i, 0, &ia_init, 0); // set initial values for IJA

  // Clear out the diagonal + one extra entry
  memset(s->a, 0, nm_sizeof[s->dtype] * ia_init);
}


char yale_storage_set_diagonal(YALE_STORAGE* s, y_size_t i, void* v) {
  memcpy(YALE_DIAG(s, nm_sizeof[s->dtype], i), v, nm_sizeof[s->dtype]); return 'r';
}


// Binary search for returning stored values
int yale_storage_binary_search(YALE_STORAGE* s, y_size_t left, y_size_t right, y_size_t key) {
  y_size_t mid = (left + right)/2, mid_j;

  if (left > right) return -1;

  YaleGetIJA(mid_j, s, mid);
  if (mid_j == key) return mid;
  else if (mid_j > key) return yale_storage_binary_search(s, left, mid-1, key);
  else return yale_storage_binary_search(s, left+1, mid, key);
}


// Binary search for returning insertion points
y_size_t yale_storage_insert_search(YALE_STORAGE* s, y_size_t left, y_size_t right, y_size_t key, bool* found) {
  y_size_t mid = (left + right)/2, mid_j;

  if (left > right) {
    *found = false;
    return right+1;
  }

  YaleGetIJA(mid_j, s, mid);
  if (mid_j == key) {
    *found = true;
    return mid;
  }
  else if (mid_j > key) return yale_storage_insert_search(s, left, mid-1, key, found);
  else return yale_storage_insert_search(s, left+1, mid, key, found);
}


// If we add n items to row i, we need to increment ija[i+1] and onward
// TODO: Add function pointer array for AddFuncs, the same way we have SetFuncs
void yale_storage_increment_ia_after(YALE_STORAGE* s, y_size_t ija_size, y_size_t i, y_size_t n) {
  y_size_t val;

  ++i;
  for (; i <= ija_size; ++i) {
    YaleGetIJA(val, s, i);
    val += n;
    YaleSetIJA(i, s, val);
  }
}


char yale_storage_set(YALE_STORAGE* s, size_t* coords, void* v) {
  y_size_t i_next = coords[0] + 1;
  y_size_t ija, ija_next, ija_size;
  y_size_t pos;
  bool found = false;
  char ins_type;

  if (coords[0] == coords[1]) return yale_storage_set_diagonal(s, coords[0], v);

  // Get IJA positions of the beginning and end of the row
  YaleGetIJA(ija,      s, coords[0]);
  YaleGetIJA(ija_next, s, i_next);

  if (ija == ija_next) { // empty row
    fprintf(stderr, "--- ija == ija_next (empty row)\n");
    ins_type = yale_vector_insert(s, ija, &(coords[1]), v, 1);
    yale_storage_increment_ia_after(s, YALE_IA_SIZE(s), coords[0], 1);
    s->ndnz++;
    return ins_type;
  }

  // non-empty row. search for coords[1] in the IJA array, between ija and ija_next
  // (including ija, not including ija_next)
  YaleGetSize(ija_size, s);
  --ija_next;

  // Do a binary search for the column
  pos = yale_storage_insert_search(s, ija, ija_next, coords[1], &found);

  if (found) {
    fprintf(stderr, "--- found replace\n");
    return yale_vector_replace(s, pos, &(coords[1]), v, 1);
  } else {
    fprintf(stderr, "--- not found insert: %d\n", pos);
    ins_type = yale_vector_insert(s, pos, &(coords[1]), v, 1);
    yale_storage_increment_ia_after(s, YALE_IA_SIZE(s), coords[0], 1);
    s->ndnz++;
    return ins_type;
  }
}


void* yale_storage_ref(YALE_STORAGE* s, size_t* coords) {
  y_size_t l, r, i_plus_one = coords[0] + 1, test_j, sz;
  int pos;

  if (coords[0] == coords[1]) return YALE_DIAG(s,nm_sizeof[s->dtype],coords[0]);

  YaleGetIJA(l, s, coords[0]);
  YaleGetIJA(r, s, i_plus_one);

  if (l == r) return YALE_A(s, nm_sizeof[s->dtype], s->shape[0]); // return zero pointer

  YaleGetSize(sz, s);
  pos = yale_storage_binary_search(s, l, r, coords[1]); // binary search for the column's location
  if (pos != -1) {
    // TODO: Technically, the above function should also tell us whether the exact value was found or not.
    YaleGetIJA(test_j, s, pos);
    if (test_j == coords[1])
      return YALE_A(s, nm_sizeof[s->dtype], pos); // Found exact value.
  }

  // return a pointer that happens to be zero
  return YALE_A(s, nm_sizeof[s->dtype], s->shape[0]);
}

#endif
