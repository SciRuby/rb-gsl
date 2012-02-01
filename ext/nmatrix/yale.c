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

#define YALE_GROWTH_CONSTANT    1.5

//#define YALE_JA_START(sptr)             (((YALE_STORAGE*)(sptr))->shape[0]+1)
#define YALE_IJA(sptr,elem_size,i)          (void*)( (char*)(((YALE_STORAGE*)(sptr))->ija) + i * elem_size )
//#define YALE_JA(sptr,dtype,j)           ((((dtype)*)((YALE_STORAGE*)(sptr))->ija)[(YALE_JA_START(sptr))+j])
#define YALE_ROW_LENGTH(sptr,elem_size,i)   (*(size_t*)YALE_IA((sptr),(elem_size),(i)+1) - *(size_t*)YALE_IJA((sptr),(elem_size),(i)))
#define YALE_A(sptr,elem_size,i)            (void*)((char*)(((YALE_STORAGE*)(sptr))->a) + elem_size * i)
#define YALE_DIAG(sptr, elem_size, i)       ( YALE_A((sptr),(elem_size),(i)) )
//#define YALE_LU(sptr,dtype,i,j)             (((dtype)*)(((YALE_STORAGE*)(sptr))->a)[ YALE_JA_START(sptr) +  ])
#define YALE_MINIMUM(sptr)                  (((YALE_STORAGE*)(sptr))->shape[0]*2 + 1) // arbitrarily defined
#define YALE_SIZE_PTR(sptr,elem_size)       (void*)((char*)((YALE_STORAGE*)(sptr))->ija + ((YALE_STORAGE*)(sptr))->shape[0]*elem_size )
#define YALE_MAX_SIZE(sptr)                 (((YALE_STORAGE*)(sptr))->shape[0] * ((YALE_STORAGE*)(sptr))->shape[1] + 1)
#define YALE_IA_SIZE(sptr)                  ((YALE_STORAGE*)(sptr))->shape[0]

// None of these next three return anything. They set a reference directly.
#define YaleGetIJA(victim,s,i)              (SetFuncs[Y_SIZE_T][(s)->index_dtype](1, &(victim), 0, YALE_IJA((s), nm_sizeof[s->index_dtype], (i)), 0))
#define YaleSetIJA(i,s,from)                (SetFuncs[s->index_dtype][Y_SIZE_T](1, YALE_IJA((s), nm_sizeof[s->index_dtype], (i)), 0, &(from), 0))
#define YaleGetSize(sz,s)                   (SetFuncs[Y_SIZE_T][(s)->index_dtype](1, &sz, 0, (YALE_SIZE_PTR((s), nm_sizeof[(s)->index_dtype])), 0))
//#define YALE_FIRST_NZ_ROW_ENTRY(sptr,elem_size,i)


void print_vectors(YALE_STORAGE* s) {
  size_t i;
  fprintf(stderr, "------------------------------\n");
  fprintf(stderr, "dtype:%s\tshape:%dx%d\tndnz:%d\tcapacity:%d\tindex_dtype:%s\n", nm_dtypestring[s->dtype], s->shape[0], s->shape[1], s->ndnz, s->capacity, nm_dtypestring[s->index_dtype]);

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
    fprintf(stderr, "%-1.3f ", *(double*)((char*)(s->a) + nm_sizeof[NM_FLOAT64]*i));
  fprintf(stderr, "\n");

  fprintf(stderr, "------------------------------\n");
}


// Determine the index dtype (which will be used for the ija vector)
int8_t yale_index_dtype(YALE_STORAGE* s) {
  if (YALE_MAX_SIZE(s) < UINT8_MAX) return NM_INT8;
  else if (YALE_MAX_SIZE(s) < UINT16_MAX) return NM_INT16;
  else if (YALE_MAX_SIZE(s) < UINT32_MAX) return NM_INT32;
  else return NM_INT64;
}


/*bool yale_vector_grow(YALE_STORAGE* s, size_t new_capacity) {
  void *new_ija, *new_a;
  y_size_t current_size;

  new_ija   =       malloc( nm_sizeof[s->index_dtype]  * new_capacity );
  new_a     =       malloc( nm_sizeof[s->dtype]        * new_capacity );

  if (!new_ija || !new_a) {
    free(new_ija); free(new_a);
    return false;
  } else {
    // get current size
    YaleGetSize(current_size, s);

    // copy contents of old vectors into new ones
    SetFuncs[s->index_dtype][s->index_dtype](current_size, new_ija, nm_sizeof[s->index_dtype], s->ija, nm_sizeof[s->index_dtype]);
    SetFuncs[s->dtype      ][s->dtype      ](current_size, new_a,   nm_sizeof[s->dtype],       s->a,   nm_sizeof[s->dtype]      );

    // free old vectors
    free(s->a);
    free(s->ija);

    // replace vectors
    s->a        = new_a;
    s->ija      = new_ija;
    s->capacity = new_capacity;

    return true;
  }
}*/


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


// Insert a value or contiguous values in the ija and a vectors (after ja and diag). Does not free anything; you are responsible!
//
// TODO: Improve this so it can handle non-contiguous element insertions efficiently.
// (For now, we can just sort the elements in the row in question.)
char yale_vector_insert(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n) {
  void* new_ija, *new_a;
  y_size_t sz, new_capacity, i;

  if (pos < s->shape[0])
    rb_raise(rb_eArgError, "vector insert pos is before beginning of ja; this should not happen");

  YaleGetSize(sz, s);

  if (sz + n > s->capacity) {

    // figure out how big to grow the vector
    new_capacity = s->capacity * YALE_GROWTH_CONSTANT;
    if (new_capacity > YALE_MAX_SIZE(s)) new_capacity = YALE_MAX_SIZE(s);

    fprintf(stderr, "yale_vector_insert: resizing to %u\n", new_capacity);
    fprintf(stderr, "\tpos=%u,j=%u,n=%u\n", pos, j[0], n);

    // allocate the new vectors
    new_ija     = malloc( nm_sizeof[s->index_dtype] * new_capacity );
    new_a       = malloc( nm_sizeof[s->dtype]       * new_capacity );

    if (!new_ija || !new_a) {
      free(new_a); free(new_ija);
      rb_raise(rb_eNoMemError, "yale sparse vectors are full and there is insufficient memory for growing them");
      return (char)false;
    }

    // copy before insertion site
    SetFuncs[s->index_dtype][s->index_dtype](pos, new_ija, nm_sizeof[s->index_dtype], s->ija, nm_sizeof[s->index_dtype]);
    SetFuncs[s->dtype      ][s->dtype      ](pos, new_a,   nm_sizeof[s->dtype],       s->a,   nm_sizeof[s->dtype]      );

    // insert
    SetFuncs[s->index_dtype][Y_SIZE_T](n, (char*)new_ija + nm_sizeof[s->index_dtype]*pos, nm_sizeof[s->index_dtype], j,   sizeof(y_size_t));
    SetFuncs[s->dtype      ][s->dtype](n, (char*)new_a   + nm_sizeof[s->dtype]*pos,       nm_sizeof[s->dtype      ], val, nm_sizeof[s->dtype]);

    // copy after insertion site
    SetFuncs[s->index_dtype][s->index_dtype](s->capacity - pos, (char*)new_ija + nm_sizeof[s->index_dtype]*(pos+n), nm_sizeof[s->index_dtype], (char*)(s->ija) + nm_sizeof[s->index_dtype]*pos, nm_sizeof[s->index_dtype]);
    SetFuncs[s->dtype      ][s->dtype      ](s->capacity - pos, (char*)new_a   + nm_sizeof[s->dtype]*(pos+n),       nm_sizeof[s->dtype],       (char*)(s->a)   + nm_sizeof[s->dtype      ]*pos, nm_sizeof[s->dtype      ]);

    return 'i';
  }

  fprintf(stderr, "yale_vector_insert: no resize needed\n");
  fprintf(stderr, "\tpos=%u,j=%u,n=%u\n", pos, j[0], n);

  // No resize required:
  // easy (but somewhat slow), just copy elements to the tail, starting at the end, one element at a time.
  // TODO: This can be made slightly more efficient, but only after the tests are written.
  for (i = 0; i < sz - pos; ++i) {
    SetFuncs[s->index_dtype][s->index_dtype](1, (char*)(s->ija) + (sz+n-i)*nm_sizeof[s->index_dtype], 0, (char*)(s->ija) + (sz-1-i)*nm_sizeof[s->index_dtype], 0);
    SetFuncs[s->dtype      ][s->dtype      ](1, (char*)(s->a)   + (sz+n-i)*nm_sizeof[s->dtype      ], 0, (char*)(s->a)   + (sz-1-i)*nm_sizeof[s->dtype      ], 0);
  }

  // Now insert the new values
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

  print_vectors(s);

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

  print_vectors(s);
}


void yale_storage_set_diagonal(YALE_STORAGE* s, y_size_t i, void* v) {
  memcpy(YALE_DIAG(s, nm_sizeof[s->dtype], i), v, nm_sizeof[s->dtype]);
}



// Places value v in column j between indices l and r in ija/a.
// Returns 'i' for insert, 'r' for replace, or NULL for failure.
char yale_storage_place_between(YALE_STORAGE* s, y_size_t l, y_size_t r, y_size_t j, void* v) {
  y_size_t m = (l + r) / 2, mj;
  fprintf(stderr, "l=%d,r=%d\n",l,r);
  if (l == 39) rb_raise(rb_eArgError, "oops");
  YaleGetIJA(mj, s, m);

  if (j == mj) {
    return yale_vector_replace(s, m, &j, v, 1);
  } else if (j < mj) {
    if (l == m) return yale_vector_insert(s, l, &j, v, 1); // and j != mj. so, insert at l
    else return yale_storage_place_between(s, l, m, j, v);
  } else { // j > mj
    if (r == m) return yale_vector_insert(s, m, &j, v, 1);
    else return yale_storage_place_between(s, m, r, j, v);
  }
}


// Insertion of one element.
void yale_storage_set(YALE_STORAGE* s, y_size_t* coords, void* v) {
  y_size_t l, r, ja_l, ja_r, i_plus_one = coords[0]+1;
  char ins_type;

  if (coords[0] == coords[1]) yale_storage_set_diagonal(s, coords[0], v);
  else {
    // Find the left and right boundaries in the j's of IJA
    // l and r are the location in IJA and A of the beginning and end of the row for i=coords[0]
    YaleGetIJA(l, s, coords[0]);
    YaleGetIJA(r, s, i_plus_one);

    fprintf(stderr, "looking for i=%u (j=%u)\tfound l=%u, r=%u", coords[0], coords[1], l, r);

    if (l != r) { // binary search insertion
      YaleGetIJA(ja_l, s, l);
      YaleGetIJA(ja_r, s, r);
      fprintf(stderr, ", ja_l=%u, ja_r=%u\n", ja_l, ja_r);

      ins_type = yale_storage_place_between(s, ja_l, ja_r, coords[1], v);
    } else { // create row
      fprintf(stderr, "\n");
      ins_type = yale_vector_insert(s, l, &(coords[1]), v, 1);
    }


    // if there was an insertion, increase the row length in IA.
    if (ins_type == 'i') {
      ++r;
      s->ndnz++;

      // reuse l and r -- don't mean the same thing anymore
      YaleSetIJA(i_plus_one, s, r); // mark that there is now an entry
      for (l = i_plus_one+1; l <= YALE_IA_SIZE(s); ++l) {
        YaleGetIJA(r, s, l);
        fprintf(stderr, "Incrementing IJA position %u: currently %u\n", l, r);
        // TODO: Write increment function pointers or figure out a faster way to do this than "get, ++, set".
        ++r;
        YaleSetIJA(l, s, r);
      }
    }
  }

  print_vectors(s);
}


void* yale_storage_ref_between(YALE_STORAGE* s, y_size_t l, y_size_t r, y_size_t j) {
  y_size_t m = (l + r) / 2, mj;
  YaleGetIJA(mj, s, m);

  if (l == m || r == m) return NULL;

  if (j == mj) {
    return s->a + m*nm_sizeof[s->dtype];
  } else if (j < mj) {
    return yale_storage_ref_between(s, l, m, j);
  } else { // j > mj
    return yale_storage_ref_between(s, m, r, j);
  }
}


void* yale_storage_ref(YALE_STORAGE* s, size_t* coords) {
  y_size_t l, r, ja_l, ja_r, i_plus_one = coords[0] + 1;
  void* ref;

  print_vectors(s);

  // assume 2D! At some point need to look at how to make an ND yale matrix, maybe?
  if (coords[0] == coords[1]) return YALE_DIAG(s,nm_sizeof[s->dtype],coords[0]);

  YaleGetIJA(l, s, coords[0]);
  YaleGetIJA(r, s, i_plus_one);

  if (l == r) return YALE_IJA(s, nm_sizeof[s->dtype], s->shape[0]);

  YaleGetIJA(ja_l, s, l);
  YaleGetIJA(ja_r, s, r);

  ref = yale_storage_ref_between(s, ja_l, ja_r, coords[1]);
  if (ref) return ref;

  // return a pointer that happens to be zero
  return YALE_IJA(s, nm_sizeof[s->dtype], s->shape[0]);
}

#endif
