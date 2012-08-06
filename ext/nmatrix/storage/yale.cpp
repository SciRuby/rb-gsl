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
// == yale.c
//
// "new yale" storage format for 2D matrices (like yale, but with
// the diagonal pulled out for O(1) access).
//
// Specifications:
// * dtype and index dtype must necessarily differ
//      * index dtype is defined by whatever unsigned type can store
//        max(rows,cols)
//      * that means vector ija stores only index dtype, but a stores
//        dtype
// * vectors must be able to grow as necessary
//      * maximum size is rows*cols+1

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
#include "yale.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

extern bool (*ElemEqEq[NUM_DTYPES][2])(const void*, const void*, const int, const int);

/*
 * Forward Declarations
 */

template <typename DType, typename IType>
static bool						ndrow_is_empty_template(const YALE_STORAGE* s, y_size_t ija, const y_size_t ija_next);

template <typename LDType, typename RDType, typename IType>
static bool						ndrow_eqeq_ndrow_template(const YALE_STORAGE* l, const YALE_STORAGE* r, y_size_t l_ija, const y_size_t l_ija_next, y_size_t r_ija, const y_size_t r_ija_next);

template <typename LDType, typename RDType, typename IType>
static bool yale_storage_eqeq_template(const YALE_STORAGE* left, const YALE_STORAGE* right);

static char						yale_storage_vector_replace_j(YALE_STORAGE* s, y_size_t pos, y_size_t* j);
static YALE_STORAGE*	yale_storage_copy_alloc_struct(const YALE_STORAGE* rhs, const dtype_t new_dtype, const y_size_t new_capacity, const y_size_t new_size);
static void						clear_diagonal_and_zero(YALE_STORAGE* s);
static void						yale_storage_increment_ia_after(YALE_STORAGE* s, y_size_t ija_size, y_size_t i, y_size_t n);
static y_size_t				yale_storage_insert_search(YALE_STORAGE* s, y_size_t left, y_size_t right, y_size_t key, bool* found);
static YALE_STORAGE*	yale_storage_alloc(dtype_t dtype, size_t* shape, size_t rank);

/*
 * Functions
 */

///////////////
// Lifecycle //
///////////////

/*
 * Documentation goes here.
 */
YALE_STORAGE* yale_storage_create(dtype_t dtype, size_t* shape, size_t rank, size_t init_capacity) {
  YALE_STORAGE* s;
  size_t max_capacity;

	// FIXME: This error should be handled in the nmatrix.c file.
  if (rank != 2) {
   	rb_raise(rb_eNotImpError, "Can only support 2D matrices");
  }

  s = yale_storage_alloc(dtype, shape, rank);
  max_capacity = storage_count_max_elements(s->rank, s->shape) - s->shape[0] + 1;

  // Set matrix capacity (and ensure its validity)
  if (init_capacity < YALE_MINIMUM(s)) {
  	s->capacity = YALE_MINIMUM(s);

  } else if (init_capacity > max_capacity) {
  	// Don't allow storage to be created larger than necessary
  	s->capacity = max_capacity;

  } else {
  	s->capacity = init_capacity;
  }

  s->ija = ALLOC_N( char, ITYPE_SIZES[s->itype] * s->capacity );
  s->a   = ALLOC_N( char, DTYPE_SIZES[s->dtype] * s->capacity );

  return s;
}

/*
 * Documentation goes here.
 *
 * FIXME: Add templating.
 */
YALE_STORAGE* yale_storage_create_from_old_yale(dtype_t dtype, size_t* shape, char* ia, char* ja, char* a, dtype_t from_dtype, itype_t from_itype) {
  YALE_STORAGE* s;
  y_size_t i = 0, p, p_next, j, ndnz = 0, pp;

  // Read through ia and ja and figure out the ndnz count.

  // Walk down rows
  for (i = shape[0]; i-- > 0;) {
/*    SetFuncs[Y_SIZE_T][from_itype](1, &p, 0, ia + ITYPE_SIZES[from_itype]*i, 0);          // p = ia[i]*/
/*    SetFuncs[Y_SIZE_T][from_itype](1, &p_next, 0, ia + ITYPE_SIZES[from_itype]*(i+1), 0); // p_next = ia[i+1]*/

    // Now walk through columns
    for (; p < p_next; ++p) {
/*      SetFuncs[Y_SIZE_T][from_itype](1, &j, 0, ja + DTYPE_SIZES[from_itype]*p, 0);        // j = ja[p]*/

      if (i != j) {
      	// entry is non-diagonal and probably nonzero
      	++ndnz;
      }
    }
  }

  // Having walked through the matrix, we now go about allocating the space for it.
  s = yale_storage_alloc(dtype, shape, 2);
  s->capacity = shape[0] + ndnz + 1;
  s->ndnz     = ndnz;

  // Setup IJA and A arrays
  s->ija = ALLOC_N( char, ITYPE_SIZES[s->itype] * s->capacity );
  s->a   = ALLOC_N( char, DTYPE_SIZES[s->dtype] * s->capacity );

  // Figure out where to start writing JA in IJA:
  pp = s->shape[0]+1;

  // Find beginning of first row
/*  SetFuncs[Y_SIZE_T][from_itype](1, &p, 0, ia, 0);          // p = ia[i]*/

  // Now fill the arrays
  for (i = s->shape[0]; i-- > 0;) {

    // Find end of row (of input)
/*    SetFuncs[Y_SIZE_T][from_itype](1, &p_next, 0, ia + DTYPE_SIZES[from_itype]*(i+1), 0); // p_next = ia[i+1]*/

    // Set the beginning of the row (of output)
/*    SetFuncs[s->itype][Y_SIZE_T](1, (char*)(s->ija) + DTYPE_SIZES[s->itype]*i, 0, &pp, 0);*/

    // Now walk through columns
    for (; p < p_next; ++p, ++pp) {
/*      SetFuncs[Y_SIZE_T][from_itype](1, &j, 0, ja + DTYPE_SIZES[from_itype]*p, 0);        // j = ja[p]*/

      if (i == j) {
/*        SetFuncs[s->dtype][from_dtype](1, (char*)(s->a) + DTYPE_SIZES[s->dtype]*i, 0, a + DTYPE_SIZES[from_dtype]*p, 0);*/
        --pp;

      } else {
/*        SetFuncs[s->itype][from_itype](1, (char*)(s->ija) + DTYPE_SIZES[s->itype]*pp, 0, ja + DTYPE_SIZES[from_itype]*p, 0);*/
/*        SetFuncs[s->dtype][from_dtype](1,             (char*)(s->a)   + DTYPE_SIZES[s->dtype]*pp,       0, a  + DTYPE_SIZES[from_dtype]*p,       0);*/
      }
    }
  }

  // Set the end of the last row
/*  SetFuncs[s->itype][Y_SIZE_T](1, (char*)(s->ija) + DTYPE_SIZES[s->itype]*i, 0, &pp, 0);*/

  // Set the zero position for our output matrix
  if (dtype == RUBYOBJ) {
  	*(VALUE*)((char*)(s->a) + DTYPE_SIZES[s->dtype]*i) = INT2FIX(0);

  } else {
  	memset((char*)(s->a) + DTYPE_SIZES[s->dtype]*i, 0, DTYPE_SIZES[s->dtype]);
  }

  return s;
}

/*
 * Documentation goes here.
 */
YALE_STORAGE* yale_storage_create_merged(const YALE_STORAGE* template, const YALE_STORAGE* other) {
  y_size_t ija, ija_next, o_ija, o_ija_next;
  y_size_t ja, size, pos;
  YALE_STORAGE* s;
  bool found;
  char ins_type;
  size_t i;

  YaleGetSize(size, template);

  // FIXME: Something needs to happen here.  I'm not sure what SMMP_MAX is.
  //s = yale_storage_copy_alloc_struct(template, template->dtype, SMMP_MAX(template->capacity, other->capacity), size);

  // set the element between D and LU (the boundary in A), which should be 0.
  memcpy((char*)(s->a) + DTYPE_SIZES[s->dtype] * s->shape[0], (char*)(template->a) + DTYPE_SIZES[template->dtype] * template->shape[0], DTYPE_SIZES[s->dtype]);

  if (other && other != template) {
  	// some operations are unary and don't need this; others are x+x and don't need this

  	for (i = s->shape[0]; i-- > 0;) {
      YaleGetIJA(ija,        s, i);
      YaleGetIJA(ija_next,   s, i+1);

      YaleGetIJA(o_ija,      other, i);
      YaleGetIJA(o_ija_next, other, i+1);

      while (o_ija < o_ija_next) {
        YaleGetIJA(o_ja,     other, o_ija);
        YaleGetIJA(ja,       s,     ija);

        if (ija == ija_next) {
        	// destination row is empty

          ins_type = yale_vector_insert(s, ija, &ja, NULL, 1, true);
          yale_storage_increment_ia_after(s, YALE_IA_SIZE(s), i, 1);
          ++(s->ndnz);
          ++ija;

          if (ins_type == 'i') {
          	++ija_next;
          }

        } else {
        	// merge positions into destination row
          pos = yale_storage_insert_search(s, ija, ija_next-1, ja, &found);

          if (!found) {
            yale_vector_insert(s, pos, &ja, NULL, 1, true);
            yale_storage_increment_ia_after(s, YALE_IA_SIZE(s), i, 1);
            ++(s->ndnz);
            if (ins_type == 'i') {
            	++ija_next;
            }
          }

          // can now set a left boundary for the next search
          ija = pos + 1;
        }

        ++o_ija;
      }
    }
  }

  return s;
}

/*
 * Destructor for yale storage
 */
void yale_storage_delete(YALE_STORAGE* s) {
  if (s) {
    YALE_STORAGE* storage = reinterpret_cast<YALE_STORAGE*>(s);
    free(storage->shape);
    free(storage->ija);
    free(storage->a);
    free(storage);
  }
}

/*
 * Empty the matrix.
 *
 * FIXME: Add templating.
 */
void yale_storage_init(YALE_STORAGE* s) {
  y_size_t IA_INIT = YALE_IA_SIZE(s)+1, i;

  // clear out IJA vector
  for (i = YALE_IA_SIZE(s) + 1; i-- > 0;) {
    //SetFuncs[s->itype][Y_SIZE_T](1, (char*)(s->ija) + i*DTYPE_SIZES[s->itype], 0, &IA_INIT, 0); // set initial values for IJA
  }

  clear_diagonal_and_zero(s);
}

/*
 * Documentation goes here.
 */
void yale_storage_mark(void* storage_base) {
  YALE_STORAGE* storage = reinterpret_cast<YALE_STORAGE*>(storage_base);
  size_t i;

  if (storage && storage->dtype == RUBYOBJ) {
  	for (i = storage->capacity; i-- > 0;) {
      rb_gc_mark(*((VALUE*)((char*)(storage->a) + i*DTYPE_SIZES[RUBYOBJ])));
    }
  }
}

/*
 * Allocates and initializes the basic struct (but not the IJA or A vectors).
 */
static YALE_STORAGE* yale_storage_alloc(dtype_t dtype, size_t* shape, size_t rank) {
  YALE_STORAGE* s;

  s = ALLOC( YALE_STORAGE );

  s->ndnz        = 0;
  s->dtype       = dtype;
  s->shape       = shape;
  s->rank        = rank;
  s->itype       = yale_storage_itype(s);

  return s;
}

///////////////
// Accessors //
///////////////


/*
 * Documentation goes here.
 */
void yale_storage_get(STORAGE* storage, SLICE* slice) {
  YALE_STORAGE* s = (YALE_STORAGE*)storage;
  rb_raise(rb_eNotImpError, "This type of yale slicing not supported yet");
}


/*
 * Documentation goes here.
 */
void* yale_storage_ref(STORAGE* storage, SLICE* slice) {
  YALE_STORAGE* s = (YALE_STORAGE*)storage;
  size_t* coords = slice->coords;

  y_size_t l, r, i_plus_one = coords[0] + 1, test_j;
  int pos;

  if (coords[0] == coords[1]) {
  	return YALE_DIAG(s,DTYPE_SIZES[s->dtype],coords[0]);
  }

  YaleGetIJA(l, s, coords[0]);
  YaleGetIJA(r, s, i_plus_one);

  if (l == r) {
  	// return zero pointer
  	return YALE_A(s, DTYPE_SIZES[s->dtype], s->shape[0]);
  }

	// binary search for the column's location
  pos = yale_storage_binary_search(s, l, r-1, coords[1]);

  if (pos != -1) {
    YaleGetIJA(test_j, s, pos);
    if (test_j == coords[1]) {
    	// Found exact value.
      return YALE_A(s, DTYPE_SIZES[s->dtype], pos);
    }
  }

  // return a pointer that happens to be zero
  return YALE_A(s, DTYPE_SIZES[s->dtype], s->shape[0]);
}

/*
 * Documentation goes here.
 */
char yale_storage_set(STORAGE* storage, SLICE* slice, void* v) {
  YALE_STORAGE* s = (YALE_STORAGE*)storage;
  size_t* coords = slice->coords;

  y_size_t i_next = coords[0] + 1;
  y_size_t ija, ija_next, ija_size;
  y_size_t pos;
  bool found = false;
  char ins_type;

  if (coords[0] == coords[1]) {
  	return yale_storage_set_diagonal(s, coords[0], v);
  }

  // Get IJA positions of the beginning and end of the row
  YaleGetIJA(ija,      s, coords[0]);
  YaleGetIJA(ija_next, s, i_next);

  if (ija == ija_next) {
  	// empty row
    ins_type = yale_vector_insert(s, ija, &(coords[1]), v, 1, false);
    yale_storage_increment_ia_after(s, YALE_IA_SIZE(s), coords[0], 1);
    s->ndnz++;

    return ins_type;
  }

  // non-empty row. search for coords[1] in the IJA array, between ija and ija_next
  // (including ija, not including ija_next)
  YaleGetSize(ija_size, s);
  //--ija_next;

  // Do a binary search for the column
  pos = yale_storage_insert_search(s, ija, ija_next-1, coords[1], &found);

  if (found) {
  	return yale_storage_vector_replace(s, pos, &(coords[1]), v, 1);
  }

  ins_type = yale_vector_insert(s, pos, &(coords[1]), v, 1, false);
  yale_storage_increment_ia_after(s, YALE_IA_SIZE(s), coords[0], 1);
  s->ndnz++;

  return ins_type;
}

///////////
// Tests //
///////////

bool yale_storage_eqeq(const STORAGE* left, const STORAGE* right) {
  LRI_DTYPE_TEMPLATE_TABLE(yale_storage_eqeq_template, bool, const YALE_STORAGE*, const YALE_STORAGE*);

  return ttable[left->dtype][right->dtype][left->itype]((const YALE_STORAGE*)left, (const YALE_STORAGE*)right);
}

/*
 * Yale eql? -- for whole-matrix comparison returning a single value.
 */
template <typename LDType, typename RDType, typename IType>
static bool yale_storage_eqeq_template(const YALE_STORAGE* left, const YALE_STORAGE* right) {
  LDType* la = reinterpret_cast<LDType*>(left->a);
  RDType* ra = reinterpret_cast<RDType*>(right->a);

  // Compare the diagonals first.
  for (size_t index = 0; index < left->shape[0]; ++index) {
    if (la[index] != ra[index]) return false;
  }

  IType* lij = reinterpret_cast<IType*>(left->ija);
  IType* rij = reinterpret_cast<IType*>(right->ija);

  for (size_t i = 0; i < left->shape[0]; ++i) {

  // Get start and end positions of row
    IType l_ija = lij[i],
          l_ija_next = lij[i+1],
          r_ija = rij[i],
          r_ija_next = rij[i+1];

    // Check to see if one row is empty and the other isn't.
    if (ndrow_is_empty_template<LDType,IType>(left, l_ija, l_ija_next)) {
      if (!ndrow_is_empty_template<RDType,IType>(right, r_ija, r_ija_next)) {
      	return false;
      }

    } else if (ndrow_is_empty_template<RDType,IType>(right, r_ija, r_ija_next)) {
    	// one is empty but the other isn't
      return false;

    } else if (!ndrow_eqeq_ndrow_template<LDType,RDType,IType>(left, right, l_ija, l_ija_next, r_ija, r_ija_next)) {
    	// Neither row is empty. Must compare the rows directly.
      return false;
    }

  }

  return true;
}

/*
 * Are two non-diagonal rows the same? We already know.
 *
 * FIXME: Add templating.
 */
template <typename LDType, typename RDType, typename IType>
static bool ndrow_eqeq_ndrow_template(const YALE_STORAGE* l, const YALE_STORAGE* r, y_size_t l_ija, const y_size_t l_ija_next, y_size_t r_ija, const y_size_t r_ija_next) {
  bool l_no_more = false, r_no_more = false;

  IType *lij = reinterpret_cast<IType*>(l->ija),
        *rij = reinterpret_cast<IType*>(r->ija);

  LDType* la = reinterpret_cast<LDType*>(l->a);
  RDType* ra = reinterpret_cast<RDType*>(r->a);

  IType l_ja = lij[l_ija],
        r_ja = rij[r_ija],
        ja;

  while (!(l_no_more && r_no_more)) {
    if (l_ja == r_ja) {

      if (ra[r_ija] != la[l_ija]) return false; // Direct comparison

      ++l_ija;
      ++r_ija;

      if (l_ija < l_ija_next) {
      	l_ja = lij[l_ija];

      } else {
      	l_no_more = true;
      }

      if (r_ija < r_ija_next) {
      	r_ja = rij[r_ija];

      } else {
      	r_no_more = true;
      }

      ja = SMMP_MIN(l_ja, r_ja);

    } else if (l_no_more || ja < l_ja) {

      if (ra[r_ija] != 0) return false;

      ++r_ija;
      if (r_ija < r_ija_next) {
      	// get next column
      	r_ja = rij[r_ija];
        ja = SMMP_MIN(l_ja, r_ja);

      } else {
      	l_no_more = true;
      }

    } else if (r_no_more || ja < r_ja) {

      if (la[l_ija] != 0) return false;

      ++l_ija;
      if (l_ija < l_ija_next) {
      	// get next column
        l_ja = lij[l_ija];
        ja = SMMP_MIN(l_ja, r_ja);
      } else {
      	l_no_more = true;
      }

    } else {
      fprintf(stderr, "Unhandled in eqeq: l_ja=%d, r_ja=%d\n", l_ja, r_ja);
    }
  }

	// every item matched
  return true;
}

/*
 * Is the non-diagonal portion of the row empty?
 */
template <typename DType, typename IType>
static bool ndrow_is_empty_template(const YALE_STORAGE* s, IType ija, const IType ija_next) {
  if (ija == ija_next) return true;

  DType* a = reinterpret_cast<DType*>(s->a);

	// do all the entries = zero?
  for (; ija < ija_next; ++ija) {
    if (a[ija] != 0) return false;
  }

  return true;
}

/////////////
// Utility //
/////////////

/*
 * Determine the index dtype (which will be used for the ija vector). This is
 * determined by matrix shape, not IJA/A vector capacity. Note that it's MAX-2
 * because UINTX_MAX and UINTX_MAX-1 are both reserved for sparse matrix
 * multiplication.
 *
 * FIXME: Needs a default case.
 */
itype_t yale_storage_itype(YALE_STORAGE* s) {
  if (YALE_MAX_SIZE(s) < UINT8_MAX-2) {
  	return UINT8;

  } else if (YALE_MAX_SIZE(s) < UINT16_MAX - 2) {
  	return UINT16;

  } else if (YALE_MAX_SIZE(s) < UINT32_MAX - 2) {
  	return UINT32;

  } else if (YALE_MAX_SIZE(s) >= UINT64_MAX - 2) {
    fprintf(stderr, "WARNING: Matrix can contain no more than %llu non-diagonal non-zero entries, or results may be unpredictable\n", UINT64_MAX - SMMP_MIN(s->shape[0], s->shape[1]) - 2);
    // TODO: Turn this into an exception somewhere else. It's pretty unlikely, but who knows.
  	return UINT64;
  }
}

/*
 * Documentation goes here.
 */
void yale_storage_print_vectors(YALE_STORAGE* s) {
  size_t i;
  fprintf(stderr, "------------------------------\n");
  fprintf(stderr, "dtype:%s\tshape:%dx%d\tndnz:%d\tcapacity:%d\titype:%s\n",
  	DTYPE_NAMES[s->dtype], s->shape[0], s->shape[1], s->ndnz, s->capacity, nm_itypestring[s->itype]);


	// This needs to be handled somewhere else.
  //if (s->capacity > 60) rb_raise(rb_eArgError, "overflow in print_vectors; cannot handle that large of a vector");
  // print indices

  fprintf(stderr, "i:\t");
  for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5lu ", (unsigned long)i);

  fprintf(stderr, "\nija:\t");
  if (YALE_MAX_SIZE(s) < UINT8_MAX)
    for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5u ", *(uint8_t*)YALE_IJA(s,ITYPE_SIZES[s->itype],i));
  else if (YALE_MAX_SIZE(s) < UINT16_MAX)
    for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5u ", *(u_int16_t*)YALE_IJA(s,ITYPE_SIZES[s->itype],i));
  else if (YALE_MAX_SIZE(s) < UINT32_MAX)
    for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5u ", *(u_int32_t*)YALE_IJA(s,ITYPE_SIZES[s->itype],i));
  else
    for (i = 0; i < s->capacity; ++i) fprintf(stderr, "%-5llu ", *(u_int64_t*)YALE_IJA(s,ITYPE_SIZES[s->itype],i));
  fprintf(stderr, "\n");

  // print values
  fprintf(stderr, "a:\t");
  for (i = 0; i < s->capacity; ++i)
    fprintf(stderr, "%-*.3g ", 5, *(double*)((char*)(s->a) + ITYPE_SIZES[FLOAT64]*i));
  fprintf(stderr, "\n");

  fprintf(stderr, "------------------------------\n");
}

/*
 * Binary search for returning stored values.
 */
int yale_storage_binary_search(YALE_STORAGE* s, y_size_t left, y_size_t right, y_size_t key) {
  y_size_t mid = (left + right)/2, mid_j;

  if (left > right) {
  	return -1;
  }

  YaleGetIJA(mid_j, s, mid);
  if (mid_j == key) {
  	return mid;

  } else if (mid_j > key) {
  	return yale_storage_binary_search(s, left, mid - 1, key);

  } else {
  	return yale_storage_binary_search(s, mid + 1, right, key);
  }
}

/*
 * Documentation goes here.
 */
char yale_storage_set_diagonal(YALE_STORAGE* s, y_size_t i, void* v) {
  memcpy(YALE_DIAG(s, DTYPE_SIZES[s->dtype], i), v, DTYPE_SIZES[s->dtype]);

  return 'r';
}

/*
 * Documentation goes here.
 */
char yale_storage_vector_replace(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n) {

  // Now insert the new values
  SetFuncs[s->itype][Y_SIZE_T](n,
                                     pos*ITYPE_SIZES[s->itype] + (char*)(s->ija),
                                     ITYPE_SIZES[s->itype],
                                     j,
                                     sizeof(y_size_t));
  SetFuncs[s->dtype][s->dtype](n,
                               pos*DTYPE_SIZES[s->dtype] + (char*)(s->a),
                               DTYPE_SIZES[s->dtype],
                               val,
                               DTYPE_SIZES[s->dtype]);

  return 'r';
}


/*
 * Documentation goes here.
 */
char yale_storage_vector_insert_resize(YALE_STORAGE* s, y_size_t current_size, y_size_t pos, y_size_t* j, void* val, y_size_t n, bool struct_only) {
  void *new_ija, *new_a;
  // Determine the new capacity for the IJA and A vectors.
  size_t new_capacity = s->capacity * YALE_GROWTH_CONSTANT;

  if (new_capacity > YALE_MAX_SIZE(s)) {
    new_capacity = YALE_MAX_SIZE(s);

    if (current_size + n > YALE_MAX_SIZE(s)) {
    	rb_raise(rb_eNoMemError, "insertion size exceeded maximum yale matrix size");
    }
  }

  if (new_capacity < current_size + n) {
  	new_capacity = current_size + n;
  }

  // Allocate the new vectors.
  new_ija     = ALLOC_N( char, ITYPE_SIZES[s->itype] * new_capacity );
  new_a       = ALLOC_N( char, DTYPE_SIZES[s->dtype]       * new_capacity );

  // Check that allocation succeeded.
  if (!new_ija || !new_a) {
    free(new_a); free(new_ija);
    rb_raise(rb_eNoMemError, "yale sparse vectors are full and there is insufficient memory for growing them");
    return (char)false;
  }

  // Copy all values prior to the insertion site to the new IJA and new A
  SetFuncs[s->itype][s->itype](pos, new_ija, ITYPE_SIZES[s->itype], s->ija, ITYPE_SIZES[s->itype]);
  if (!struct_only) {
    SetFuncs[s->dtype][s->dtype](pos, new_a, DTYPE_SIZES[s->dtype], s->a, DTYPE_SIZES[s->dtype]);
  }

  // Copy all values subsequent to the insertion site to the new IJA and new A, leaving room (size n) for insertion.
  SetFuncs[s->itype][s->itype](current_size-pos+n-1, (char*)new_ija + ITYPE_SIZES[s->itype]*(pos+n), ITYPE_SIZES[s->itype],
  	(char*)(s->ija) + ITYPE_SIZES[s->itype]*pos, ITYPE_SIZES[s->itype]);

  if (!struct_only) {
    SetFuncs[s->dtype][s->dtype](current_size-pos+n-1, (char*)new_a + DTYPE_SIZES[s->dtype]*(pos+n), DTYPE_SIZES[s->dtype], (char*)(s->a) + DTYPE_SIZES[s->dtype]*pos, DTYPE_SIZES[s->dtype]);
  }

  s->capacity = new_capacity;

  free(s->ija);
  free(s->a);

  s->ija = new_ija;
  s->a   = new_a;

  return 'i';
}

/*
 * Insert a value or contiguous values in the ija and a vectors (after ja and
 * diag). Does not free anything; you are responsible!
 *
 * TODO: Improve this so it can handle non-contiguous element insertions
 *	efficiently. For now, we can just sort the elements in the row in
 *	question.)
 */
char yale_storage_vector_insert(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n, bool struct_only) {
  y_size_t sz, i;

  if (pos < s->shape[0]) {
    rb_raise(rb_eArgError, "vector insert pos is before beginning of ja; this should not happen");
  }

  YaleGetSize(sz, s);

  if (sz + n > s->capacity) {
  	yale_vector_insert_resize(s, sz, pos, j, val, n, struct_only);

  } else {

    /*
     * No resize required:
     * easy (but somewhat slow), just copy elements to the tail, starting at
     * the end, one element at a time.
     *
     * TODO: This can be made slightly more efficient, but only after the tests
     *	are written.
     */
    for (i = 0; i < sz - pos; ++i) {
      SetFuncs[s->itype][s->itype](1, (char*)(s->ija) + (sz+n-1-i)*ITYPE_SIZES[s->itype], 0, (char*)(s->ija) + (sz-1-i)*ITYPE_SIZES[s->itype], 0);

      if (!struct_only) {
        SetFuncs[s->dtype    ][s->dtype      ](1, (char*)(s->a)   + (sz+n-1-i)*DTYPE_SIZES[s->dtype      ], 0, (char*)(s->a)   + (sz-1-i)*DTYPE_SIZES[s->dtype      ], 0);
      }
    }
  }

  // Now insert the new values.
  if (struct_only) {
  	yale_vector_replace_j(s, pos, j);

  } else {
  	yale_storage_vector_replace(s, pos, j, val, n);
  }

  return 'i';
}

/*
 * Clear out the D portion of the A vector (clearing the diagonal and setting
 * the zero value).
 *
 * Note: This sets a literal 0 value. If your dtype is RUBYOBJ (a Ruby object),
 * it'll actually be INT2FIX(0) instead of a string of NULLs.
 */
static void clear_diagonal_and_zero(YALE_STORAGE* s) {
  y_size_t i;
  // Clear out the diagonal + one extra entry
  if (s->dtype == RUBYOBJ) {
    for (i = 0; i < YALE_IA_SIZE(s)+1; ++i) // insert Ruby zeros
      *(VALUE*)( (char*)(s->a) + i*DTYPE_SIZES[s->dtype] ) = INT2FIX(0);
  } else { // just insert regular zeros
    memset(s->a, 0, DTYPE_SIZES[s->dtype] * (YALE_IA_SIZE(s)+1));
  }
}

/*
 * If we add n items to row i, we need to increment ija[i+1] and onward.
 *
 * TODO: Add function pointer array for AddFuncs, the same way we have SetFuncs.
 */
static void yale_storage_increment_ia_after(YALE_STORAGE* s, y_size_t ija_size, y_size_t i, y_size_t n) {
  y_size_t val;

  ++i;
  for (; i <= ija_size; ++i) {
    YaleGetIJA(val, s, i);
    val += n;
    YaleSetIJA(i, s, val);
  }
}

/*
 * Binary search for returning insertion points.
 */
static y_size_t yale_storage_insert_search(YALE_STORAGE* s, y_size_t left, y_size_t right, y_size_t key, bool* found) {
  y_size_t mid = (left + right)/2, mid_j;

  if (left > right) {
    *found = false;
    return left;
  }

  YaleGetIJA(mid_j, s, mid);
  if (mid_j == key) {
    *found = true;
    return mid;

  } else if (mid_j > key) {
  	return yale_storage_insert_search(s, left, mid-1, key, found);

  } else {
  	return yale_storage_insert_search(s, mid+1, right, key, found);
  }
}

/*
 * Just like yale_storage_vector_replace, but doesn't replace the contents of the cell,
 * only the column index.
 */
static char yale_storage_vector_replace_j(YALE_STORAGE* s, y_size_t pos, y_size_t* j) {
  SetFuncs[s->itype][Y_SIZE_T](1,
                                     pos*ITYPE_SIZES[s->itype] + (char*)(s->ija),
                                     ITYPE_SIZES[s->itype],
                                     j,
                                     sizeof(y_size_t));

  return 'r';
}

/////////////////////////
// Copying and Casting //
/////////////////////////

/*
 * Copy constructor for changing dtypes.
 */
STORAGE* yale_storage_cast_copy(const STORAGE* rhs, dtype_t new_dtype) {
  LR_DTYPE_TEMPLATE_TABLE(yale_storage_cast_copy_template, YALE_STORAGE*, const YALE_STORAGE*, dtype_t);

  return (STORAGE*)ttable[new_dtype][rhs->dtype]((YALE_STORAGE*)rhs, new_dtype);
}


/*
 * Templated copy constructor for changing dtypes.
 */
template <typename LDType, typename RDType>
YALE_STORAGE* yale_storage_cast_copy_template(const YALE_STORAGE* rhs, dtype_t new_dtype) {

  // Allocate a new structure
  size_t size = yale_storage_get_size(rhs);
  YALE_STORAGE* lhs = yale_storage_copy_alloc_struct(rhs, new_dtype, rhs->capacity, size);

  if (rhs->dtype == new_dtype) {  // FIXME: Test if this condition is actually faster; second condition should work just as well.

    memcpy(lhs->a, rhs->a, size * DTYPE_SIZES[new_dtype]);

  } else {

    DType*    la = reinterpret_cast<LDType*>(lhs->a);
    NewDType* ra = reinterpret_cast<RDType*>(rhs->a);

    for (size_t index = 0; index < size; ++index) {
      la[index] = static_cast<LDType>(ra[index]);
    }

  }

  return lhs;
}


/*
 * Returns size of Yale storage as a size_t (no matter what the itype is).
 */
static inline size_t yale_storage_get_size(const YALE_STORAGE* storage) {
  ITYPE_TEMPLATE_TABLE(yale_storage_get_size_template, size_t, const YALE_STORAGE*);

  return ttable[storage->itype](storage);
}


/*
 * Template access for getting the size of Yale storage.
 */
template <typename IType>
static inline size_t yale_storage_get_size_template(const YALE_STORAGE* storage) {
  return static_cast<size_t>(reinterpret_cast<IType*>(storage->ija)[ storage->shape[0] ]);
}



/*
 * Copy constructor.
 */
YALE_STORAGE* yale_storage_copy(YALE_STORAGE* rhs) {
  y_size_t size;
  YALE_STORAGE* lhs;

  YaleGetSize(size, rhs);
  lhs = yale_storage_copy_alloc_struct(rhs, rhs->dtype, rhs->capacity, size);

  // Now copy the contents -- but only within the boundaries set by the size. Leave
  // the rest uninitialized.
  memcpy(lhs->a, rhs->a, size * DTYPE_SIZES[lhs->dtype]);

  return lhs;
}

/*
 * Allocate for a copy or copy-cast operation, and copy the IJA portion of the
 * matrix (the structure).
 */
static YALE_STORAGE* yale_storage_copy_alloc_struct(const YALE_STORAGE* rhs, const dtype_t new_dtype, const y_size_t new_capacity, const y_size_t new_size) {
  YALE_STORAGE* lhs = ALLOC( YALE_STORAGE );
  lhs->rank         = rhs->rank;
  lhs->shape        = ALLOC_N( size_t, lhs->rank );
  memcpy(lhs->shape, rhs->shape, lhs->rank * sizeof(size_t));
  lhs->itype        = rhs->itype;
  lhs->capacity     = new_capacity;
  lhs->dtype        = new_dtype;
  lhs->ndnz         = rhs->ndnz;

  lhs->ija          = ALLOC_N( char, ITYPE_SIZES[lhs->itype] * lhs->capacity );
  lhs->a            = ALLOC_N( char, DTYPE_SIZES[lhs->dtype] * lhs->capacity );

  // Now copy the contents -- but only within the boundaries set by the size. Leave
  // the rest uninitialized.
  memcpy(lhs->ija, rhs->ija, yale_storage_get_size(rhs) * ITYPE_SIZES[lhs->itype]); // indices

  return lhs;
}


STORAGE* yale_storage_matrix_multiply(STORAGE_PAIR casted_storage, size_t* resulting_shape, bool vector) {
  LI_DTYPE_TEMPLATE_TABLE(yale_storage_matrix_multiply_template, NMATRIX*, STORAGE_PAIR, size_t*, bool);

  YALE_STORAGE* storage_access = reinterpret_cast<YALE_STORAGE*>(casted_storage.left);

  return ttable[storage_access->dtype][storage_access->itype](casted_storage, resulting_shape, vector);
}


template <typename DType, typename IType>
static STORAGE* yale_storage_matrix_multiply_template(STORAGE_PAIR casted_storage, size_t* resulting_shape, bool vector) {
  YALE_STORAGE *left  = (YALE_STORAGE*)(casted_storage.left),
               *right = (YALE_STORAGE*)(casted_storage.right);

  // We can safely get dtype from the casted matrices; post-condition of binary_storage_cast_alloc is that dtype is the
  // same for left and right.
  // int8_t dtype = left->dtype;

  // Create result storage.
  YALE_STORAGE* result = create_yale_storage(left->dtype, resulting_shape, 2, left->capacity + right->capacity);
  init_yale_storage(result);

  NAMED_ITYPE_TEMPLATE_TABLE(symbmm_table, symbmm, void, const unsigned int, const unsigned int, const IType*, const IType*, const bool, const IType*, const IType*, const bool, IType*, const bool);
  NAMED_LI_DTYPE_TEMPLATE_TABLE(numbmm_table, numbmm, void, const unsigned int, const unsigned int, const IType*, const IType*, const DType*, const bool, const IType* ib, const IType* jb, const DType* b, const bool diagb, IType* ic, IType* jc, DType* c, const bool);
  NAMED_LI_DTYPE_TEMPLATE_TABLE(smmp_sort_columns_table, smmp_sort_columns, const unsigned int n, const IType* ia, IType* ja, DType* a);

  // Do the multiplication
  symbmm_table[left->itype](result->shape[0], result->shape[1], left->ija, left->ija, true, right->ija, right->ija, true, result->ija, true);
  nubmm_table[left->dtype][left->itype](result->shape[0], result->shape[1], left->ija, left->ija, left->a, true, right->ija, right->ija, right->a, true, result->ija, result->ija, result->a, true);

  // Sort the columns
  smmp_sort_columns_table[left->dtype][left->itype](result->shape[0], result->ija, result->ija, result->a);

  return reinterpret_cast<STORAGE*>(result);
}