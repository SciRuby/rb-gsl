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
// == yale.h
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

#ifndef YALE_H
#define YALE_H

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

#include "types.h"

#include "data/data.h"

#include "common.h"

/*
 * Macros
 */

#define YALE_GROWTH_CONSTANT 1.5

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
#define YaleGetIJA(victim,s,i)              //(SetFuncs[Y_SIZE_T][(s)->itype](1, &(victim), 0, YALE_IJA((s), DTYPE_SIZES[s->itype], (i)), 0))
#define YaleSetIJA(i,s,from)                //(SetFuncs[s->itype][Y_SIZE_T](1, YALE_IJA((s), DTYPE_SIZES[s->itype], (i)), 0, &(from), 0))
#define YaleGetSize(sz,s)                   //(SetFuncs[Y_SIZE_T][((YALE_STORAGE*)s)->itype](1, &sz, 0, (YALE_SIZE_PTR(((YALE_STORAGE*)s), DTYPE_SIZES[((YALE_STORAGE*)s)->itype])), 0))
//#define YALE_FIRST_NZ_ROW_ENTRY(sptr,elem_size,i)

/*
 * Types
 */

typedef struct {
	// Common elements found in all storage types.  Must not be re-arranged.
	dtype_t	dtype;
	size_t	rank;
	size_t*	shape;
	size_t*	offset;
	
	// Yale storage specific elements.
	void* a;

	// Strictly non-diagonal non-zero count!
	size_t ndnz;

	size_t	capacity;
	itype_t	itype;
	void*		ija;
} YALE_STORAGE;

/*
 * Data
 */
 

/*
 * Functions
 */

///////////////
// Lifecycle //
///////////////

YALE_STORAGE* yale_storage_create(dtype_t dtype, size_t* shape, size_t rank, size_t init_capacity);
YALE_STORAGE* yale_storage_create_from_old_yale(dtype_t dtype, size_t* shape, char* ia, char* ja, char* a, dtype_t from_dtype, itype_t from_itype);
YALE_STORAGE*	yale_storage_create_merged(const YALE_STORAGE* merge_template, const YALE_STORAGE* other);
void					yale_storage_delete(STORAGE* s);
void					yale_storage_init(YALE_STORAGE* s);
void					yale_storage_mark(void*);

///////////////
// Accessors //
///////////////

void* yale_storage_get(STORAGE* s, SLICE* slice);
void*	yale_storage_ref(STORAGE* s, SLICE* slice);
char	yale_storage_set(STORAGE* s, SLICE* slice, void* v);

///////////
// Tests //
///////////

bool yale_storage_eqeq(const YALE_STORAGE* left, const YALE_STORAGE* right);

/////////////
// Utility //
/////////////

dtype_t	yale_storage_itype(YALE_STORAGE* s);
void		yale_storage_print_vectors(YALE_STORAGE* s);

int yale_storage_binary_search(YALE_STORAGE* s, y_size_t left, y_size_t right, y_size_t key);

char yale_storage_set_diagonal(YALE_STORAGE* s, y_size_t i, void* v);

char yale_storage_vector_replace(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n);
char yale_storage_vector_insert_resize(YALE_STORAGE* s, y_size_t current_size, y_size_t pos, y_size_t* j, void* val, y_size_t n, bool struct_only);
char yale_storage_vector_insert(YALE_STORAGE* s, y_size_t pos, y_size_t* j, void* val, y_size_t n, bool struct_only);

/////////////////////////
// Copying and Casting //
/////////////////////////

STORAGE* yale_storage_cast_copy(const STORAGE* rhs, dtype_t new_dtype);
YALE_STORAGE* yale_storage_copy(YALE_STORAGE* rhs);

#endif
