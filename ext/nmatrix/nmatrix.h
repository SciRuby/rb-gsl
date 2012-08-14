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
// == nmatrix.h
//

#ifndef NMATRIX_H
#define NMATRIX_H

/*
 * Standard Includes
 */

#include <math.h>
#include <ruby.h>
#include <string.h>

#ifdef BENCHMARK
	// SOURCE: http://stackoverflow.com/questions/2349776/how-can-i-benchmark-a-c-program-easily
	
	#include <sys/time.h>
	#include <sys/resource.h>
#endif

/*
 * Project Includes
 */

#include "types.h"

#include "data/data.h"

#include "math.h"

#include "storage/storage.h"

/*
 * Macros
 */

#define RUBY_ZERO INT2FIX(0)


#ifndef SIZEOF_INT
 #error SIZEOF_INT undefined
#else
 #if SIZEOF_INT == 8
	#define DEFAULT_DTYPE  INT64
	#define SIZE_T         INT64
 #else
	#if SIZEOF_INT == 4
		#define DEFAULT_DTYPE INT32
		#define SIZE_T        INT32
	#else
	 #if SIZEOF_INT == 2
		#define DEFAULT_DTYPE INT16
		#define SIZE_T        INT16
	 #else
	  #error Unhandled SIZEOF_INT -- please #define SIZE_T and DEFAULT_DTYPE manually.
	 #endif
	#endif
 #endif
#endif

// The definition above is partially overridden by the definition below when size_t is undefined.
#ifndef HAVE_SIZE_T
 typedef u_int64_t  size_t;
 #define SIZE_T     INT64
#endif

#define NM_MAX_RANK 15

#define UnwrapNMatrix(obj,var)  Data_Get_Struct(obj, NMATRIX, var)

#define NM_STORAGE(val)         (NM_STRUCT(val)->storage)
#define NM_LIST_STORAGE(val)    ((LIST_STORAGE*)(NM_STORAGE(val)))
#define NM_YALE_STORAGE(val)    ((YALE_STORAGE*)(NM_STORAGE(val)))
#define NM_DENSE_STORAGE(val)   ((DENSE_STORAGE*)(NM_STORAGE(val)))


//#define NM_PTR(a, p)            ((a)->ptr+(p)*nm_sizeof[(a)->type])
//#define NM_PTR_TYPE(val,type)   (type)(((struct numeric_matrix*)DATA_PTR(val))->ptr)

#define NM_DENSE_SRC(val)       (NM_DENSE_STORAGE(val)->src)
#define NM_STRUCT(val)          (reinterpret_cast<NMATRIX*>(DATA_PTR(val)))
#define NM_RANK(val)            (NM_STORAGE(val)->rank)
#define NM_DTYPE(val)           (NM_STORAGE(val)->dtype)
#define NM_ITYPE(val)           (NM_YALE_STORAGE(val)->itype)
#define NM_STYPE(val)           (NM_STRUCT(val)->stype)
#define NM_SHAPE(val,i)         (NM_STORAGE(val)->shape[(i)])
#define NM_SHAPE0(val)          (NM_STORAGE(val)->shape[0])
#define NM_SHAPE1(val)          (NM_STORAGE(val)->shape[1])

#define NM_DENSE_COUNT(val)     (storage_count_max_elements(NM_DENSE_STORAGE(val)))
#define NM_SIZEOF_DTYPE(val)    (nm_sizeof[NM_DTYPE(val)])
#define NM_REF(val,slice)      (RefFuncs[NM_STYPE(val)]( NM_STORAGE(val), slice, NM_SIZEOF_DTYPE(val) ))
    
#define NM_MAX(a,b) (((a)>(b))?(a):(b))
#define NM_MIN(a,b) (((a)>(b))?(b):(a))
#define NM_SWAP(a,b,tmp) {(tmp)=(a);(a)=(b);(b)=(tmp);}

#define NM_CHECK_ALLOC(x) if (!x) rb_raise(rb_eNoMemError, "insufficient memory");

#define CheckNMatrixType(v)   if (TYPE(v) != T_DATA || (RDATA(v)->dfree != (RUBY_DATA_FUNC)nm_delete && RDATA(v)->dfree != (RUBY_DATA_FUNC)nm_delete_ref)) rb_raise(rb_eTypeError, "expected NMatrix on left-hand side of operation");

#define NM_IsNMatrix(obj) \
  (rb_obj_is_kind_of(obj, cNMatrix) == Qtrue)

#define NM_IsNVector(obj) \
  (rb_obj_is_kind_of(obj, cNVector) == Qtrue)

// FIXME: What should this actually be?
//#define NM_INDEX_TYPES  NM_FLOAT32
#define NM_INDEX_TYPES  5

/*
 * Types
 */

// Two vectors and a capacity.
typedef struct {
  void*  ija;
  void*  a;
  size_t capacity;
} VECTOR;

typedef struct NMATRIX {
	// Method of storage (csc, dense, etc).
	stype_t		stype;
	// Pointer to storage struct.
	STORAGE*	storage;
} NMATRIX;

/*
 * Data
 */

#ifndef NMATRIX_C
	extern VALUE cNMatrix;
#endif

//extern uint8_t (*MathHomOps_b[5])(const uint8_t, const uint8_t);
//extern int64_t (*MathHomOps_i64[5])(const int64_t, const int64_t);
//extern int32_t (*MathHomOps_i32[5])(const int32_t, const int32_t);
//extern int16_t (*MathHomOps_i16[5])(const int16_t, const int16_t);
//extern int8_t (*MathHomOps_i8[5])(const int8_t, const int8_t);
//extern float (*MathHomOps_f32[5])(const float, const float);
//extern double (*MathHomOps_f64[5])(const double, const double);
//extern complex64 (*MathHomOps_c64[5])(const complex64, const complex64);
//extern complex128 (*MathHomOps_c128[5])(const complex128, const complex128);
//extern rational32 (*MathHomOps_r32[5])(rational32, rational32);
//extern rational64 (*MathHomOps_r64[5])(rational64, rational64);
//extern rational128 (*MathHomOps_r128[5])(rational128, rational128);
//extern VALUE (*MathHomOps_v[5])(const VALUE, const VALUE);
//extern void (*SmmpSortColumns[15][7])(const unsigned int, const void *, void *, void *);
//extern void (*Transp[15][7])(const unsigned int, const unsigned int, const void *, const void *, const void *, const bool, void *, void *, void *, const bool);
//extern void (*DetExact[15])(const int, const void *, const int, void *);

//extern int (*EwDenseHom[15])(const void *, const void *, void *, const int, enum MathHomOps);
//extern int (*EwDenseBool[15])(const void *, const void *, void *, const int, const enum MathBoolOps);
//extern int (*EwDenseBit[15])(const void *, const void *, void *, const int, const enum MathBitOps);
//extern int (*EwYaleHom[15][7])(const void *, const void *, void *, const void *, const void *, void *, const unsigned int, const unsigned int, const enum MathHomOps);
//extern int (*EwYaleBool[15][7])(const void *, const void *, void *, const void *, const void *, void *, const unsigned int, const unsigned int, const enum MathBoolOps);
//extern int (*EwYaleBit[15][7])(const void *, const void *, void *, const void *, const void *, void *, const unsigned int, const unsigned int, const enum MathBitOps);

/*
 * Functions
 */

// FIXME: Does this belong here?
void transp(size_t n, size_t m, void* ia, void* ja, bool diaga, void* a, void* ib, void* jb, void* b, bool move, int8_t itype, dtype_t dtype);

void cast_copy_value_single(void* to, const void* from, dtype_t l_dtype, dtype_t r_dtype);

extern "C" {
	void Init_nmatrix();
	
	// External API
	VALUE rb_nmatrix_dense_create(dtype_t dtype, size_t* shape, size_t rank, void* elements, size_t length);
	VALUE rb_nvector_dense_create(dtype_t dtype, void* elements, size_t length);
}

#endif // NMATRIX_H
