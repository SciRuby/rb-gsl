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

#include "nmatrix.h"


// for when we need to return array indices.
// TODO: Make it so this is less than or equal to size_t
typedef uint32_t    y_size_t;
#define Y_SIZE_T    NM_INT32

#define YALE_GROWTH_CONSTANT    1.5


//#define YALE_JA_START(sptr)             (((YALE_STORAGE*)(sptr))->shape[0]+1)
#define YALE_IJA(sptr,elem_size,i)          ( (((YALE_STORAGE*)(sptr))->ija) + i * elem_size )
//#define YALE_JA(sptr,dtype,j)           ((((dtype)*)((YALE_STORAGE*)(sptr))->ija)[(YALE_JA_START(sptr))+j])
#define YALE_ROW_LENGTH(sptr,elem_size,i)   (*(size_t*)YALE_IA((sptr),(elem_size),(i)+1) - *(size_t*)YALE_IJA((sptr),(elem_size),(i)))
#define YALE_A(sptr,elem_size,i)            (((YALE_STORAGE*)(sptr))->a + elem_size * i)
#define YALE_DIAG(sptr, elem_size, i)       ( YALE_A((sptr),(elem_size),(i)) )
#define YALE_LU(sptr,dtype,i,j)             (((dtype)*)(((YALE_STORAGE*)(sptr))->a)[ YALE_JA_START(sptr) +  ])
#define YALE_MINIMUM(sptr)                  ((YALE_STORAGE*)(sptr)->shape[0]*2 + 1) // arbitrarily defined
#define YALE_SIZE_PTR(sptr,elem_size)       ((YALE_STORAGE*)(sptr)->ija + (YALE_STORAGE*)(sptr)->shape[0]*elem_size )

// None of these next three return anything. They set a reference directly.
#define YaleGetIJA(victim,s,i)              (SetFuncs[Y_SIZE_T][(s)->index_dtype](1, &(victim), 0, YALE_IJA((s), nm_sizeof[s->index_dtype], (i)), 0))
#define YaleSetIJA(i,s,from)                (SetFuncs[s->index_dtype][Y_SIZE_T](1, YALE_IJA((s), nm_sizeof[s->index_dtype], (i)), 0, &(from), 0))
#define YaleGetSize(sz,s)                   (SetFuncs[Y_SIZE_T][s->index_dtype](1, &sz, 0, YALE_SIZE_PTR(s, nm_sizeof[s->index_dtype]), 0))

#define YALE_MAX_SIZE(sptr)                 ((YALE_STORAGE*)(sptr)->shape[0] * (YALE_STORAGE*)(sptr)->shape[1] + 1)
//#define YALE_FIRST_NZ_ROW_ENTRY(sptr,elem_size,i)


// Determine the index dtype (which will be used for the ija vector)
int8_t yale_index_dtype(YALE_STORAGE* s) {
  if (YALE_MAX_SIZE(s) < UINT8_MAX) return NM_INT8;
  else if (YALE_MAX_SIZE(s) < UINT16_MAX) return NM_INT16;
  else if (YALE_MAX_SIZE(s) < UINT32_MAX) return NM_INT32;
  else if (YALE_MAX_SIZE(s) < UINT64_MAX) return NM_INT64;
  else rb_raise(rb_eNotImpError, "new yale matrix maximum size exceeded");
}


// Create a zero-filled vector
bool yale_vector_grow(YALE_STORAGE* s, size_t new_capacity) {
  void *new_ija, *new_a;
  y_size_t current_size;

  new_ija   =       malloc( nm_sizeof[s->index_dtype]  * new_capacity );
  new_a     =       malloc( nm_sizeof[s->dtype]        * new_capacity );

  if (!new_ija || !new_a) {
    free(new_ija); free(new_a);
    return false;
  } else {
    // get current size
    YaleGetSize(sz, s);

    // copy contents of old vectors into new ones
    SetFuncs[s->index_dtype][s->index_dtype](sz, new_ija, nm_sizeof[s->index_dtype], s->ija, nm_sizeof[s->index_dtype]);
    SetFuncs[s->dtype      ][s->dtype      ](sz, new_a,   nm_sizeof[s->dtype],       s->a,   nm_sizeof[s->dtype]      );

    // free old vectors
    free(s->a);
    free(s->ija);

    // replace vectors
    s->a        = new_a;
    s->ija      = new_ija;
    s->capacity = new_capacity;

    return true;
  }
}


char yale_vector_replace(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n) {
  // Now insert the new values
  SetFuncs[s->index_dtype][Y_SIZE_T](n, s->ija + pos*nm_sizeof[s->index_dtype], nm_sizeof[s->index_dtype], j,   nm_sizeof[Y_SIZE_T]);
  SetFuncs[s->dtype      ][s->dtype](n, s->a   + pos*nm_sizeof[s->dtype      ], nm_sizeof[s->dtype      ], val, nm_sizeof[s->dtype]);
  return 'r';
}


// Insert a value or contiguous values in the ija and a vectors (after ja and diag). Does not free anything; you are responsible!
//
// TODO: Improve this so it can handle non-contiguous element insertions efficiently.
// (For now, we can just sort the elements in the row in question.)
char yale_vector_insert(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n) {
  void* new_ija, *new_a;
  y_size_t rel_pos = pos - s->shape[0], sz, new_capacity, i;

  if (pos < s->shape[0]) rb_raise(rb_eArgError, "vector insert pos is before beginning of ja; this should not happen");

  YaleGetSize(sz, s);

  if (sz + n > s->capacity) {

    // figure out how big to grow the vector
    new_capacity = s->capacity * YALE_GROWTH_CONSTANT;
    if (new_capacity > YALE_MAX_SIZE(s)) new_capacity = YALE_MAX_SIZE(s);

    // allocate the new vectors
    new_ija     = malloc( nm_sizeof[s->index_dtype] * new_capacity );
    new_a       = malloc( nm_sizeof[s->dtype])      * new_capacity );

    if (!new_ija || !new_a) {
      free(new_a); free(new_ija);
      rb_raise(rb_eNoMemError, "yale sparse vectors are full and there is insufficient memory for growing them");
      return (char)false;
    }

    // copy before insertion site
    SetFuncs[s->index_dtype][s->index_dtype](pos, new_ija, nm_sizeof[s->index_dtype], s->ija, nm_sizeof[s->index_dtype]);
    SetFuncs[s->dtype      ][s->dtype      ](pos, new_a,   nm_sizeof[s->dtype],       s->a,   nm_sizeof[s->dtype]      );

    // insert
    SetFuncs[s->index_dtype][Y_SIZE_T](n, new_ija + nm_sizeof[s->index_dtype]*pos, nm_sizeof[s->index_dtype], j,   nm_sizeof[s->index_dtype]);
    SetFuncs[s->dtype      ][s->dtype](n, new_a   + nm_sizeof[s->dtype]*pos,       nm_sizeof[s->dtype      ], val, nm_sizeof[s->dtype      ]);

    // copy after insertion site
    SetFuncs[s->index_dtype][s->index_dtype](pos, new_ija + nm_sizeof[s->index_dtype]*(pos+n), nm_sizeof[s->index_dtype], s->ija + nm_sizeof[s->index_dtype]*pos, nm_sizeof[s->index_dtype]);
    SetFuncs[s->dtype      ][s->dtype      ](pos, new_a   + nm_sizeof[s->dtype]*(pos+n),       nm_sizeof[s->dtype],       s->a   + nm_sizeof[s->dtype      ]*pos, nm_sizeof[s->dtype      ]);

    return 'i';
  }

  // No resize required:
  // easy (but somewhat slow), just copy elements to the tail, starting at the end, one element at a time.
  // TODO: This can be made slightly more efficient, but only after the tests are written.
  for (i = 0; i < sz - pos; ++i) {
    SetFuncs[s->index_dtype][s->index_dtype](1, s->ija + (sz+n-i)*nm_sizeof[s->index_dtype], 0, s->ija + (sz-1-i)*nm_sizeof[s->index_dtype], 0);
    SetFuncs[s->dtype      ][s->dtype      ](1, s->a   + (sz+n-i)*nm_sizeof[s->dtype      ], 0, s->a   + (sz-1-i)*nm_sizeof[s->dtype      ], 0);
  }

  // Now insert the new values
  yale_vector_replace(s, pos, j, val, n);

  return 'i';
}


YALE_STORAGE* create_yale_storage(int8_t dtype, size_t* shape, size_t rank, size_t init_capacity) {
  YALE_STORAGE* s;

  if (rank > 2) rb_raise(rb_eNotImpError, "Can only support 2D matrices");

  if (!(s = malloc(sizeof(YALE_STORAGE)))) return NULL;
  s->ndnz        = 0;
  s->capacity    = init_capacity;
  s->dtype       = dtype;
  s->shape       = shape;
  s->rank        = rank;
  s->index_dtype = yale_index_dtype(s);

  if (!(s->ija = malloc(nm_sizeof[s->index_dtype] * YALE_MINIMUM(s)))) {
    free(s->shape);
    free(s);
    return NULL;
  } else {
    if (!(s->a = malloc(nm_sizeof[s->dtype] * YALE_MINIMUM(s)))) {
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
  size_t i;

  // Clear out the diagonal + one extra entry, and also the list of rows.
  for (i = 0; i <= s->shape[0]; ++i) {
    memset( YALE_DIAG(s, nm_sizeof[s->dtype], i), 0, nm_sizeof[s->dtype] );
    memset( YALE_IA(s, nm_sizeof[s->dtype], i),   0, nm_sizeof[s->dtype] );
  }
}


void yale_storage_insert_diagonal(YALE_STORAGE* s, y_size_t i, void* v) {
  memcpy(YALE_DIAG(s, nm_sizeof[s->dtype], i), v, nm_sizeof[s->dtype]);
}

// bool yale_vector_insert(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, size_t n)


y_size_t yale_storage_find(YALE_STORAGE* s, y_size_t* coords) {
  y_size_t left, right, right_ia = coords[0]+1;
  YaleGetIJA(left, s, coords[0]);
  YaleGetIJA(right, s, right_ia);

  if (left == right) {
    fprintf(stderr, "left == right: %d\n", left);
  } else {

  }
}


char yale_storage_place_between(YALE_STORAGE* s, y_size_t l, y_size_t r, y_size_t j, void* v) {
  y_size_t m = (l + r) / 2, mj;
  YaleGetIJA(mj, s, m);

  if (j == mj) {
    return yale_vector_replace(s, m, j, v, 1);
  } else if (j < mj) {
    if (l == m) return yale_vector_insert(s, l, j, v, 1); // and j != mj. so, insert at l
    else return yale_storage_place_between(s, l, m, j, v);
  } else { // j > mj
    if (r == m) return yale_vector_insert(s, m, j, v, 1);
    else return yale_storage_place_between(s, m, r, j, v);
  }
}


// Insertion of one element.
void yale_storage_insert(YALE_STORAGE* s, y_size_t* coords, void* v) {
  y_size_t l, r, m, lj, rj, tmpj = coords[0]+1, tmpk;
  char ins_type;

  if (coords[0] == coords[1]) yale_storage_insert_diagonal(s, coords[0], v);
  else {
    // Find the left and right boundaries in the j's of IJA
    YaleGetIJA(lj, s, coords[0]);
    YaleGetIJA(rj, s, tmpj);
    if (lj != rj) { // no entries in that row.
      // binary search insertion
      ins_type = yale_storage_place_between(s, lj, rj, coords[1], v);

    } else {
      ins_type = yale_vector_insert(s, lj, coords[1], v, 1);

    }
    // if there was an insertion, increase the row length in IA.
    if (ins_type == 'i') {
      ++rj;
      YaleSetIJA(tmpj, s, rj); // mark that there is now an entry
    }
  }
}


void* yale_storage_ref_between(YALE_STORAGE* s, y_size_t l, y_size_t r, y_size_t j) {
  y_size_t m = (l + r) / 2, mj;
  YaleGetIJA(mj, s, m);

  if (l == m || r == m) return NULL;

  if (j == mj) {
    return YALE_A(s, nm_sizeof[s->dtype], m);
  } else if (j < mj) {
    return yale_storage_ref_between(s, l, m, j);
  } else { // j > mj
    return yale_storage_ref_between(s, m, r, j);
  }
}


void* yale_storage_ref(YALE_STORAGE* s, size_t* coords, size_t elem_size) {
  y_size_t l, r, tmpj = coords[0] + 1;

  // assume 2D! At some point need to look at how to make an ND yale matrix, maybe?
  if (coords[0] == coords[1]) return YALE_DIAG(s,elem_size,coords[0]);

  YaleGetIJA(l, s, coords[0]);
  YaleGetIJA(r, s, tmpj);
  return yale_storage_ref_between(s, l, r, coords[1]);
}

#endif
