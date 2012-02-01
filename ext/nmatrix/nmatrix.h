/* Derived mostly from NArray, with modifications */

#ifndef NMATRIX_H
#define NMATRIX_H

#include "nmatrix_config.h"

// Matrix multiplication
#include <cblas.h>

#include <math.h>

#ifdef HAVE_STDBOOL_H
# include <stdbool.h>
#else
typedef char    bool;
# define true    1;
# define false   0;
#endif

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef HAVE_STDINT_H
# include <stdint.h>
#endif

#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>
#endif


/*
  Data types used in NArray / NMatrix :
  Please modify these types if your system has any different type.
*/


/* NM_BYTE : unsigned 8-bit integer */
#ifndef HAVE_U_INT8_T
# ifdef HAVE_UINT8_T
typedef uint8_t			u_int8_t;
# else
typedef unsigned char		u_int8_t;
# endif
#endif

//#ifndef HAVE_INT8_T
//typedef char                   int8_t;
//#endif

#ifndef HAVE_INT16_T
# if SIZEOF_SHORT == 2
typedef short                  int16_t;
# else
---->> Please define int16_t manually because sizeof(short) != 2. <<----
# endif
#endif /* HAVE_INT16_T */

#ifndef HAVE_INT32_T
# if SIZEOF_LONG == 4
typedef long                   int32_t;
# else
#  if SIZEOF_INT == 4
typedef int                    int32_t;
#  else
---->> Please define int32_t manually because sizeof(long) != 4. <<----
#  endif
# endif
#endif /* HAVE_INT32_T */

/* unsigned 32-bit integer */
#ifndef HAVE_U_INT32_T
# ifdef HAVE_UINT32_T
typedef uint32_t			u_int32_t;
# else
#  if SIZEOF_LONG == 4
typedef unsigned long                   u_int32_t;
#  else
#   if SIZEOF_INT == 4
typedef unsigned int                    u_int32_t;
#   else
---->> Please define u_int32_t manually because sizeof(long) != 4. <<----
#   endif
#  endif
# endif
#endif /* HAVE_U_INT32_T */

#ifndef HAVE_INT64_T
# if SIZEOF_QUAD == 8
typedef quad                   int64_t;
# else
#  if SIZEOF_LONG == 8
typedef long                   int64_t;
#  else
---->> Please define int64_t manually because sizeof(quad) != 8. <<----
#  endif
# endif
#endif /* HAVE_INT64_T */

/* unsigned 64-bit integer */
#ifndef HAVE_U_INT64_T
# ifdef HAVE_UINT64_T
typedef uint64_t            u_int64_t;
# else
#  if SIZEOF_QUAD == 8
typedef unsigned quad       u_int64_t;
#  else
#   if SIZEOF_LONG == 8
typedef unsigned long       u_int64_t;
#   else
---->> Please define u_int64_t manually because sizeof(quad) != 8. <<----
#   endif
#  endif
# endif
#endif /* HAVE_U_INT64_T */


typedef struct { float r,i; } complex64;
typedef struct { double r,i; } complex128;
typedef struct { int16_t n,d; } rational32;
typedef struct { int32_t n,d; } rational64;
typedef struct { int64_t n,d; } rational128;

#include "dtypes.h"

#if SIZEOF_INT == 8
# define DEFAULT_DTYPE  NM_INT64
#else
# if SIZEOF_INT == 4
#  define DEFAULT_DTYPE NM_INT32
# else
#  define DEFAULT_DTYPE NM_INT16
# endif
#endif


enum NMatrix_STypes {
  S_DENSE,
  S_LIST,
  S_COMPRESSED,
  S_TYPES
};



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


// two vectors and a capacity
typedef struct y_vector {
  void*  ija;
  void*  a;
  size_t capacity;
} VECTOR;


typedef struct common_s { // Common elements found in all _s types.
  int8_t    dtype;
  size_t    rank;
  size_t*   shape;
} STORAGE;


typedef struct list_s {
  int8_t    dtype;
  size_t    rank;
  size_t*   shape;
  void*     default_val;
  LIST*     rows;
} LIST_STORAGE;


typedef struct dense_s {
  int8_t    dtype;
  size_t    rank;
  size_t*   shape;
  void*     elements;
} DENSE_STORAGE;


typedef struct yale_s {
  int8_t    dtype;
  size_t    rank;
  size_t*   shape;
  size_t    ndnz; // strictly non-diagonal non-zero count!
  size_t    capacity;
  int8_t    index_dtype;
  void*     ija;
  void*     a;
} YALE_STORAGE;


typedef struct numeric_matrix {
  //int8_t   dtype;             /* data type (int8, int32, rational128, etc) */
  int8_t   stype;             /* method of storage (csc, dense, etc) */
  STORAGE* storage;           /* pointer to storage struct */
  // VALUE    ref;            /* NMatrix object wrapping this structure */
} NMATRIX;




#ifndef NMATRIX_C
extern VALUE cNMatrix;

extern const int nm_sizeof[NM_TYPES+1];
#endif


#define NM_MAX_RANK 15

#define UnwrapNMatrix(obj,var)  Data_Get_Struct(obj, struct numeric_matrix, var)
#define IsNMatrix(obj)          (rb_obj_is_kind_of(obj, CNMatrix)==Qtrue)

#define NM_STORAGE(val)         (((struct numeric_matrix*)DATA_PTR(val))->storage)
//#define NM_PTR(a, p)            ((a)->ptr+(p)*nm_sizeof[(a)->type])
#define NM_STRUCT(val)          ((struct numeric_matrix*)DATA_PTR(val))
//#define NM_PTR_TYPE(val,type)   (type)(((struct numeric_matrix*)DATA_PTR(val))->ptr)
#define NM_RANK(val)            (((STORAGE*)(NM_STORAGE(val)))->rank)
#define NM_DTYPE(val)           (((STORAGE*)(NM_STORAGE(val)))->dtype)
#define NM_STYPE(val)           (((struct numeric_matrix*)DATA_PTR(val))->stype)
#define NM_SHAPE(val,i)         (((STORAGE*)(NM_STORAGE(val)))->shape[(i)])
#define NM_SHAPE0(val)          (((struct numeric_matrix*)DATA_PTR(val))->shape[0])
#define NM_SHAPE1(val)          (((struct numeric_matrix*)DATA_PTR(val))->shape[1])
#define NM_SIZEOF_DTYPE(val)    (nm_sizeof[NM_DTYPE(val)])
#define NM_REF(val,coords)      (RefFuncs[NM_STYPE(val)]( NM_STORAGE(val), coords, NM_SIZEOF_DTYPE(val) ))

#define NM_IsNMatrix(obj) (rb_obj_is_kind_of(obj, cNMatrix)==Qtrue)
#define NM_IsArray(obj)   (TYPE(obj)==T_ARRAY || rb_obj_is_kind_of(obj,cNMatrix)==Qtrue)
#define NM_IsROBJ(d) ((d)->dtype==NM_ROBJ)
#define NM_IsINTEGER(a) \
    (NM_DTYPE(a)==NM_BYTE || NM_DTYPE(a)==NM_INT8 || NM_DTYPE(a)==NM_INT16 || NM_DTYPE(a)==NM_INT32 || NM_DTYPE(a)==NM_INT64)
#define NM_IsCOMPLEX(a) \
    (NM_DTYPE(a)==NM_COMPLEX32 || NM_DTYPE(a)==NM_COMPLEX64)
#define NM_MAX(a,b) (((a)>(b))?(a):(b))
#define NM_SWAP(a,b,tmp) {(tmp)=(a);(a)=(b);(b)=(tmp);}

#define NUM2REAL(v) NUM2DBL( rb_funcall((v),nm_id_real,0) )
#define NUM2IMAG(v) NUM2DBL( rb_funcall((v),nm_id_imag,0) )

#define NUM2NUMER(v) NUM2INT( rb_funcall((v), nm_id_numer,0) )
#define NUM2DENOM(v) NUM2INT( rb_funcall((v), nm_id_denom,0) )

#define IS_NUMERIC(v)   (FIXNUM_P(v) || TYPE(v) == T_FLOAT || TYPE(v) == T_COMPLEX || TYPE(v) == T_RATIONAL)
#define IS_STRING(v)    (TYPE(v) == T_STRING)

/* Local */

typedef union {
  u_int8_t b[2];
  int16_t s;
} nm_size16_t;

typedef union {
  u_int8_t b[4];
  int32_t  i;
  float    f;
} nm_size32_t;

typedef union {
  u_int8_t b[8];
  int64_t  q;
  float    f[2];
  double   d;
} nm_size64_t;

typedef union {
  u_int8_t b[16];
  double   d[2];
} nm_size128_t;

#if !defined RSTRING_LEN
#define RSTRING_LEN(a) RSTRING(a)->len
#endif
#if !defined RSTRING_PTR
#define RSTRING_PTR(a) RSTRING(a)->ptr
#endif
#if !defined RARRAY_LEN
#define RARRAY_LEN(a) RARRAY(a)->len
#endif
#if !defined RARRAY_PTR
#define RARRAY_PTR(a) RARRAY(a)->ptr
#endif


typedef void     (*nm_setfunc_t[NM_TYPES][NM_TYPES])();
typedef VALUE    (*nm_stype_ref_t[S_TYPES])();
typedef VALUE    (*nm_create_t[S_TYPES])();
typedef void     (*nm_delete_t[S_TYPES])();
typedef STORAGE* (*nm_copy_s_t[S_TYPES])();
typedef void     (*nm_gemm_t[NM_TYPES])(); // TODO: Write for non-dense?
//typedef void (*nm_setsf_t[S_TYPES][S_TYPES])();
//typedef void (*nm_setdf_t[NM_DTYPES][NM_DTYPES])();

extern nm_setfunc_t SetFuncs;
extern ID nm_id_real, nm_id_imag;
extern ID nm_id_denom, nm_id_numer;


/* dense.c */
DENSE_STORAGE*  create_dense_storage(int8_t dtype, size_t* shape, size_t rank);
void            delete_dense_storage(DENSE_STORAGE* s);
DENSE_STORAGE*  copy_dense_storage(DENSE_STORAGE* rhs);

size_t          count_dense_storage_elements(DENSE_STORAGE* s);

size_t          dense_storage_pos(DENSE_STORAGE* s, size_t* coords);
void*           dense_storage_get(DENSE_STORAGE* s, size_t* coords);
void            dense_storage_set(DENSE_STORAGE* s, size_t* coords, void* val);

/* list.c */
LIST_STORAGE*   create_list_storage(int8_t dtype, size_t* shape, size_t rank, void* init_val);
void            delete_list_storage(LIST_STORAGE* s);
LIST_STORAGE*   copy_list_storage(LIST_STORAGE* rhs);

void*           list_storage_get(LIST_STORAGE* s, size_t* coords);
void*           list_storage_insert(LIST_STORAGE* s, size_t* coords, void* val);
void*           list_storage_remove(LIST_STORAGE* s, size_t* coords);

LIST*           create_list();
void            delete_list(LIST* list, size_t recursions);
void            copy_list_contents(LIST* lhs, LIST* rhs, int8_t dtype, size_t recursions);
NODE*           list_find(LIST* list, size_t key);
NODE*           list_find_preceding_from(NODE* prev, size_t key);
NODE*           list_find_nearest_from(NODE* prev, size_t key);
void*           list_remove(LIST* list, size_t key);
NODE*           list_insert(LIST* list, bool replace, size_t key, void* val);
NODE*           list_insert_after(NODE* node, size_t key, void* val);
void            list_print_int(LIST* list);


/* nmatrix.c */
int8_t nm_dtypestring_to_dtype(VALUE str);
int8_t nm_dtypesymbol_to_dtype(VALUE sym);
int8_t nm_stypestring_to_stype(VALUE str);
int8_t nm_stypesymbol_to_stype(VALUE sym);
int8_t nm_guess_dtype(VALUE v);
size_t* nm_interpret_shape_arg(VALUE arg, size_t* rank);
VALUE nm_dense_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, VALUE self);
VALUE nm_list_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, VALUE self);
VALUE nm_init(int argc, VALUE* argv, VALUE self);
// VALUE nm_init_copy(VALUE copy, VALUE original);
NMATRIX* nm_create(int8_t stype, void* storage);
void Init_nmatrix();

#endif
