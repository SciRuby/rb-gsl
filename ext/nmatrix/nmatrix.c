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
// == nmatrix.c
//
#ifndef NMATRIX_C
# define NMATRIX_C

#include <ruby.h>

#include "nmatrix.h"
//#include "types.h"

VALUE cNMatrix, cNVector;
VALUE nm_eDataTypeError, nm_eStorageTypeError;

ID nm_id_real, nm_id_imag;
ID nm_id_numer, nm_id_denom;
ID nm_id_transpose, nm_id_no_transpose, nm_id_complex_conjugate; // cblas
ID nm_id_list, nm_id_dense;
ID nm_id_mult, nm_id_multeq;
ID nm_id_add;

#include "dtypes.c"

#ifdef BENCHMARK
double get_time() {
  struct timeval t;
  struct timezone tzp;
  gettimeofday(&t, &tzp);
  return t.tv_sec + t.tv_usec*1e-6;
}
#endif


const char *nm_stypestring[] = {
  "dense",
  "list",
  "yale",
  "stypes"
};


nm_delete_t DeleteFuncs = {
  delete_dense_storage,
  delete_list_storage,
  delete_yale_storage
};


nm_mark_t MarkFuncs = {
  mark_dense_storage,
  mark_list_storage,
  mark_yale_storage
};


nm_gemv_t GemvFuncs = {
  NULL,
  cblas_bgemv_,
  cblas_i8gemv_,
  cblas_i16gemv_,
  cblas_i32gemv_,
  cblas_i64gemv_,
  cblas_sgemv_,
  cblas_dgemv_,
  cblas_cgemv_,
  cblas_zgemv_,
  cblas_r32gemv_,
  cblas_r64gemv_,
  cblas_r128gemv_,
  NULL
};


nm_gemm_t GemmFuncs = { // by NM_TYPES
  NULL,
  cblas_bgemm_,
  cblas_i8gemm_,
  cblas_i16gemm_,
  cblas_i32gemm_,
  cblas_i64gemm_,
  cblas_sgemm_,
  cblas_dgemm_,
  cblas_cgemm_,
  cblas_zgemm_,
  cblas_r32gemm_,
  cblas_r64gemm_,
  cblas_r128gemm_,
  cblas_vgemm_
};

static void TransposeTypeErr(y_size_t n, y_size_t m, YALE_PARAM A, YALE_PARAM B, bool move) {
  rb_raise(nm_eDataTypeError, "illegal operation with this matrix type");
}


// First dimension is dtype, second dimension is index dtype (so lots of nulls)
nm_smmp_transpose_t SparseTransposeFuncs = {
  {TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_NONE
  {TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_BYTE
  {TransposeTypeErr, TransposeTypeErr, i8_i8_transp_, i16_i8_transp_, i32_i8_transp_, i64_i8_transp_, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_INT8
  {TransposeTypeErr, TransposeTypeErr, i8_i16_transp_, i16_i16_transp_, i32_i16_transp_, i64_i16_transp_, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_INT16
  {TransposeTypeErr, TransposeTypeErr, i8_i32_transp_, i16_i32_transp_, i32_i32_transp_, i64_i32_transp_, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_INT32
  {TransposeTypeErr, TransposeTypeErr, i8_i64_transp_, i16_i64_transp_, i32_i64_transp_, i64_i64_transp_, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_INT64
  {TransposeTypeErr, TransposeTypeErr, i8_f32_transp_, i16_f32_transp_, i32_f32_transp_, i64_f32_transp_, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_FLOAT32
  {TransposeTypeErr, TransposeTypeErr, i8_f64_transp_, i16_f64_transp_, i32_f64_transp_, i64_f64_transp_, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_FLOAT64
  {TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_COMPLEX64
  {TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_COMPLEX128
  {TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_RATIONAL32
  {TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_RATIONAL64
  {TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}, // NM_RATIONAL128
  {TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr, TransposeTypeErr}  // NM_ROBJ
};

static void SmmpTypeErr(y_size_t n, y_size_t m, YALE_PARAM A, YALE_PARAM B, YALE_PARAM C) {
  rb_raise(nm_eDataTypeError, "illegal operation with this matrix type");
}

// First dimension is dtype, second dimension is index dtype (so lots of nulls)
nm_smmp_t SmmpFuncs = {
  {SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_NONE
  {SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_BYTE
  {SmmpTypeErr, SmmpTypeErr, i8_i8_smmp, i16_i8_smmp, i32_i8_smmp, i64_i8_smmp, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_INT8
  {SmmpTypeErr, SmmpTypeErr, i8_i16_smmp, i16_i16_smmp, i32_i16_smmp, i64_i16_smmp, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_INT16
  {SmmpTypeErr, SmmpTypeErr, i8_i32_smmp, i16_i32_smmp, i32_i32_smmp, i64_i32_smmp, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_INT32
  {SmmpTypeErr, SmmpTypeErr, i8_i64_smmp, i16_i64_smmp, i32_i64_smmp, i64_i64_smmp, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_INT64
  {SmmpTypeErr, SmmpTypeErr, i8_f32_smmp, i16_f32_smmp, i32_f32_smmp, i64_f32_smmp, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_FLOAT32
  {SmmpTypeErr, SmmpTypeErr, i8_f64_smmp, i16_f64_smmp, i32_f64_smmp, i64_f64_smmp, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_FLOAT64
  {SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_COMPLEX64
  {SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_COMPLEX128
  {SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_RATIONAL32
  {SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_RATIONAL64
  {SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}, // NM_RATIONAL128
  {SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr, SmmpTypeErr}  // NM_ROBJ
};


static inline DENSE_PARAM cblas_params_for_multiply(const DENSE_STORAGE* left, const DENSE_STORAGE* right, const DENSE_STORAGE* result, bool vector) {
  DENSE_PARAM p;

  p.A = left->elements;
  p.B = right->elements;    // for vector, this is actually x
  p.C = result->elements;   // vector Y

  p.M = left->shape[0];
  p.lda = left->shape[1];

  if (vector) {
    p.N = left->shape[1];

    p.ldb = 1;  // incX
    p.ldc = 1;  // incY
  } else {
    p.N = right->shape[1];
    p.K = left->shape[1];

    p.ldb = right->shape[1];
    p.ldc = result->shape[1];
  }

  switch(left->dtype) {
  case NM_FLOAT32:
  case NM_FLOAT64:
    p.alpha.d[0] = 1.0;
    p.beta.d[0] = 0.0;
    break;

  case NM_COMPLEX64:
    p.alpha.c[0].r = 1.0;
    p.alpha.c[0].i = 0.0;
    p.beta.c[0].r = 0.0;
    p.beta.c[0].i = 0.0;
    break;

  case NM_COMPLEX128:
    p.alpha.z.r = 1.0;
    p.alpha.z.i = 0.0;
    p.beta.z.r = 0.0;
    p.beta.z.i = 0.0;
    break;

  case NM_BYTE:
    p.alpha.b[0] = 1;
    p.beta.b[0] = 0;
    break;

  case NM_INT8:
  case NM_INT16:
  case NM_INT32:
  case NM_INT64:
    p.alpha.i[0] = 1;
    p.beta.i[0]  = 0;
    break;

  case NM_RATIONAL32:
    p.alpha.r[0].n = 1;
    p.alpha.r[0].d = 1;
    p.beta.r[0].n = 0;
    p.beta.r[0].d = 1;
    break;

  case NM_RATIONAL64:
    p.alpha.ra[0].n = 1;
    p.alpha.ra[0].d = 1;
    p.beta.ra[0].n = 0;
    p.beta.ra[0].d = 1;
    break;

  case NM_RATIONAL128:
    p.alpha.rat.n = 1;
    p.alpha.rat.d = 1;
    p.beta.rat.n = 0;
    p.beta.rat.d = 1;
    break;

  case NM_ROBJ:
    p.alpha.v[0] = RUBY_ZERO;
    p.beta.v[0]  = INT2FIX(1);
    break;

  default:
    rb_raise(nm_eDataTypeError, "unexpected dtype");

  }

  return p;
}


static NMATRIX* multiply_matrix_dense_casted(STORAGE_PAIR casted_storage, size_t* resulting_shape, bool vector) {
  DENSE_STORAGE *left  = (DENSE_STORAGE*)(casted_storage.left),
                *right = (DENSE_STORAGE*)(casted_storage.right),
                *result;

  // We can safely get dtype from the casted matrices; post-condition of binary_storage_cast_alloc is that dtype is the
  // same for left and right.
  int8_t dtype = left->dtype;

  // Create result storage.
  result = create_dense_storage(dtype, resulting_shape, 2, NULL, 0);

  // Do the multiplication
  if (vector) GemvFuncs[dtype](CblasRowMajor, CblasNoTrans, cblas_params_for_multiply(left, right, result, true));
  else        GemmFuncs[dtype](CblasRowMajor, CblasNoTrans, CblasNoTrans, cblas_params_for_multiply(left, right, result, false));

  return nm_create(S_DENSE, result);
}


static NMATRIX* multiply_matrix_yale_casted(STORAGE_PAIR casted_storage, size_t* resulting_shape, bool vector) {
  YALE_STORAGE *left  = (YALE_STORAGE*)(casted_storage.left),
               *right = (YALE_STORAGE*)(casted_storage.right),
               *result;
  YALE_PARAM A, B, C;

  // We can safely get dtype from the casted matrices; post-condition of binary_storage_cast_alloc is that dtype is the
  // same for left and right.
  int8_t dtype = left->dtype;

  // Create result storage.
  result = create_yale_storage(dtype, resulting_shape, 2, left->capacity + right->capacity);
  init_yale_storage(result);

  // Set multiplication parameters
  A.ia = A.ja = left->ija;
  A.a  = left->a;
  B.ia = B.ja = right->ija;
  B.a  = right->a;
  C.ia = C.ja = result->ija;
  C.a  = result->a;

  A.diag = B.diag = C.diag = true;

  // Do the multiplication
  SmmpFuncs[dtype][left->index_dtype](result->shape[0], result->shape[1], A, B, C);

  return nm_create(S_YALE, result);
}


static NMATRIX* multiply_matrix_list_casted(STORAGE_PAIR casted_storage, size_t* resulting_shape) {
  rb_raise(rb_eNotImpError, "multiplication not implemented for list-of-list matrices");
  free(resulting_shape);
  return NULL;
}


nm_matrix_multiply_op_t CastedMultiplyFuncs = {
  multiply_matrix_dense_casted,
  multiply_matrix_list_casted,
  multiply_matrix_yale_casted
};


static NMATRIX* elementwise_dense_casted(STORAGE_PAIR casted_storage, char op) {
  DENSE_STORAGE *left  = (DENSE_STORAGE*)(casted_storage.left),
                *right = (DENSE_STORAGE*)(casted_storage.right),
                *result;

  // We can safely get dtype from the casted matrices; post-condition of binary_storage_cast_alloc is that dtype is the
  // same for left and right.
  size_t i;
  int8_t dtype = left->dtype;

  // Setup matrix shape for result
  size_t* shape = ALLOC_N(size_t, left->rank);
  for (i = 0; i < left->rank; ++i) shape[i] = left->shape[i];

  // Create result storage.
  result = create_dense_storage(dtype, shape, left->rank, NULL, 0);

  // Do the operation
  // TODO: Operation

  return nm_create(S_DENSE, result);
}


static NMATRIX* elementwise_list_casted(STORAGE_PAIR casted_storage, char op) {
  rb_raise(rb_eNotImpError, "elementwise operations not implemented for list-of-list matrices");
  return NULL;
}


mm_elementwise_binary_op_t CastedElementwiseFuncs = {
  elementwise_dense_casted,
  elementwise_list_casted,
  elementwise_yale_casted
};


static void nm_delete(NMATRIX* mat) {
  DeleteFuncs[mat->stype](mat->storage);
}


static STORAGE* nm_dense_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, size_t init_val_len, VALUE self) {
  return (STORAGE*)(create_dense_storage(dtype, shape, rank, init_val, init_val_len));
}

static STORAGE* nm_list_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, size_t init_val_len, VALUE self) {
  if (init_val_len > 1) {
    rb_raise(rb_eArgError, "list storage needs initial size, not initial value\n");
    return NULL;
  }
  return (STORAGE*)(create_list_storage(dtype, shape, rank, init_val));
}


static STORAGE* nm_yale_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, size_t init_val_len, VALUE self) {
  YALE_STORAGE* s;

  if (init_val_len > 1) {
    rb_raise(rb_eArgError, "list storage needs initial size, not initial value\n");
    return NULL;
  }

  s = create_yale_storage(dtype, shape, rank, *(size_t*)init_val);
  init_yale_storage(s);
  free(init_val);

  if (!s) rb_raise(rb_eNoMemError, "Yale allocation failed");

  return (STORAGE*)(s);
  //return Data_Wrap_Struct(self, NULL, nm_delete, matrix);
}


nm_create_storage_t CreateFuncs = {
  nm_dense_new,
  nm_list_new,
  nm_yale_new
};


nm_cast_copy_storage_t CastCopyFuncs = {
  cast_copy_dense_storage,
  cast_copy_list_storage,
  cast_copy_yale_storage
};


nm_stype_ref_t RefFuncs = {
  dense_storage_get,
  list_storage_get,
  yale_storage_ref
};


VALUE nm_dense_set(STORAGE* s, size_t* coords, VALUE val) {
  void* v = ALLOCA_N(char, nm_sizeof[s->dtype]);
  SetFuncs[s->dtype][NM_ROBJ](1, v, 0, &val, 0);
  dense_storage_set( (DENSE_STORAGE*)s, coords, v );
  return val;
}


// Should work exactly the same as nm_dense_set.
VALUE nm_yale_set(STORAGE* s, size_t* coords, VALUE val) {
  void* v = ALLOCA_N(char, nm_sizeof[s->dtype]);
  SetFuncs[s->dtype][NM_ROBJ](1, v, 0, &val, 0);
  yale_storage_set( (YALE_STORAGE*)s, coords, v );
  return val;
}


// TODO: Why can't you be more like your brothers, nm_dense_set and nm_yale_set?
VALUE nm_list_set(STORAGE* s, size_t* coords, VALUE val) {
  void *v = ALLOC_N(char, nm_sizeof[s->dtype]), *rm;
  LIST_STORAGE* ls = (LIST_STORAGE*)s;

  //fprintf(stderr, "    create_val: %p\n", v);

  SetFuncs[s->dtype][NM_ROBJ](1, v, 0, &val, 0);

  if (!memcmp(ls->default_val, v, nm_sizeof[s->dtype])) {
    // User asked to insert default_value, which is actually node *removal*.
    // So let's do that instead.

    rm = list_storage_remove( ls, coords );

    //if (rm) fprintf(stderr, "    remove_val: %p\n", rm);

    if (rm) free(rm);
    return val;

  } else if (list_storage_insert( ls, coords, v ))    return val;
  return Qnil;
  // No need to free; the list keeps v.
}



nm_stype_ins_t InsFuncs = {
  nm_dense_set,
  nm_list_set,
  nm_yale_set,
};



// Converts a typestring to a typecode for storage. Only looks at the first three characters.
int8_t nm_stypestring_to_stype(VALUE str) {
  int8_t i;
  for (i = 0; i < S_TYPES; ++i)
    if ( !strncmp(RSTRING_PTR(str), nm_stypestring[i], 3) ) return i;
  return S_DENSE;
}

int8_t nm_stypesymbol_to_stype(VALUE sym) {
  int8_t i;
  for (i = 0; i < S_TYPES; ++i)
    if (SYM2ID(sym) == rb_intern(nm_stypestring[i])) return i;
  return S_DENSE;
}


int8_t nm_dtypestring_to_dtype(VALUE str) {
  int8_t i;
  for (i = 0; i < NM_TYPES; ++i)
    if ( !strncmp(RSTRING_PTR(str), nm_dtypestring[i], RSTRING_LEN(str)) ) return i;
  return NM_NONE;
}

int8_t nm_dtypesymbol_to_dtype(VALUE sym) {
  int8_t i;
  for (i = 0; i < NM_TYPES; ++i)
    if (SYM2ID(sym) == rb_intern(nm_dtypestring[i])) return i;
  return NM_NONE;
}


// TODO: Probably needs some work for Bignum.
int8_t nm_guess_dtype(VALUE v) {
  switch(TYPE(v)) {
  case T_TRUE:
  case T_FALSE:
    return NM_BYTE;
  case T_STRING:
    if (RSTRING_LEN(v) == 1) return NM_BYTE;
    else                     return NM_NONE;

#if SIZEOF_INT == 8
  case T_FIXNUM:
    return NM_INT64;
  case T_RATIONAL:
    return NM_RATIONAL128;
#else
# if SIZEOF_INT == 4
  case T_FIXNUM:
    return NM_INT32;
  case T_RATIONAL:
    return NM_RATIONAL64;
#else
  case T_FIXNUM:
    return NM_INT16;
  case T_RATIONAL:
    return NM_RATIONAL32;
# endif
#endif

  case T_BIGNUM:
    return NM_INT64;

#if SIZEOF_FLOAT == 4
  case T_COMPLEX:
    return NM_COMPLEX128;
  case T_FLOAT:
    return NM_FLOAT64;
#else
# if SIZEOF_FLOAT == 2
  case T_COMPLEX:
    return NM_COMPLEX64;
  case T_FLOAT:
    return NM_FLOAT32;
# endif
#endif

  case T_ARRAY: // may be passed for dense -- for now, just look at the first element.
    return nm_guess_dtype(RARRAY_PTR(v)[0]);
    // TODO: Look at entire array for most specific type.

  case T_NIL:
  default:
    return NM_NONE;
  }
}


// Read the shape argument to NMatrix.new, which may be either an array or a single number.
// Second argument is where the shape is stored at the end of the function; returns the rank.
// You are responsible for freeing shape!
size_t* nm_interpret_shape_arg(VALUE arg, size_t* rank) {
  size_t i;
  size_t* shape;

  if (TYPE(arg) == T_ARRAY) {
    *rank = RARRAY_LEN(arg);
    shape = ALLOC_N(size_t, *rank);
    for (i = 0; i < *rank; ++i)
      shape[i]  = (size_t)(FIX2UINT(RARRAY_PTR(arg)[i]));
  } else if (FIXNUM_P(arg)) {
    *rank = 2;
    shape = ALLOC_N(size_t, *rank);
    for (i = 0; i < *rank; ++i)
      shape[i]  = (size_t)(FIX2UINT(arg));
  } else {
    *rank  = 0;
    shape = NULL;
    rb_raise(rb_eArgError, "Expected an array of numbers or a single fixnum for matrix shape");
  }

  return shape;
}


// argv will be either 1 or 2 elements. If 1, could be either initial or dtype. If 2, is initial and dtype.
// This function returns the dtype.
int8_t nm_interpret_dtype(int argc, VALUE* argv, int8_t stype) {
  int offset = 0; // if argc == 1
  if (argc == 2) offset = 1;
  else if (argc != 1) rb_raise(rb_eArgError, "Need an initial value or a dtype");

  if (SYMBOL_P(argv[offset])) return nm_dtypesymbol_to_dtype(argv[offset]);
  else if (IS_STRING(argv[offset])) return nm_dtypestring_to_dtype(StringValue(argv[offset]));
  else if (stype == S_YALE) rb_raise(rb_eArgError, "yale requires dtype");
  else return nm_guess_dtype(argv[0]);

  return NM_NONE;
}

int8_t nm_interpret_stype(VALUE arg) {
  if (SYMBOL_P(arg)) return nm_stypesymbol_to_stype(arg);
  else if (IS_STRING(arg)) return nm_stypestring_to_stype(StringValue(arg));
  else rb_raise(rb_eArgError, "Expected storage type");
  return S_DENSE;
}


void* nm_interpret_initial_value(VALUE arg, int8_t dtype) {
  void* init_val;

  if (TYPE(arg) == T_ARRAY) { // array
    init_val = ALLOC_N(char, nm_sizeof[dtype] * RARRAY_LEN(arg));
    SetFuncs[dtype][NM_ROBJ](RARRAY_LEN(arg), init_val, nm_sizeof[dtype], RARRAY_PTR(arg), nm_sizeof[NM_ROBJ]);
  } else { // single value
    init_val = ALLOC_N(char, nm_sizeof[dtype]);
    SetFuncs[dtype][NM_ROBJ](1, init_val, 0, &arg, 0);
  }

  return init_val;
}


size_t* nm_interpret_initial_capacity(VALUE arg) {
  size_t* init_cap = ALLOC(size_t);
  *init_cap = FIX2UINT(arg);
  return init_cap;
}


static VALUE nm_init(int argc, VALUE* argv, VALUE nm) {
  char    ZERO = 0;
  VALUE   QNIL = Qnil;
  int8_t  dtype, stype, offset = 0;
  size_t  rank;
  size_t* shape;
  size_t  init_val_len = 0;
  void*   init_val = NULL;
  NMATRIX* nmatrix;

  // READ ARGUMENTS

  //fprintf(stderr, "Called nmatrix new with %d arguments\n", argc);

  if (argc < 2 || argc > 4) { rb_raise(rb_eArgError, "Expected 2, 3, or 4 arguments"); return Qnil; }

  if (!SYMBOL_P(argv[0]) && !IS_STRING(argv[0])) {
    stype  = S_DENSE;
  } else {
    stype  = nm_interpret_stype(argv[0]);                        // 0: String or Symbol
    offset = 1;
  }
  shape    = nm_interpret_shape_arg(argv[offset], &rank);        // 1: String or Symbol
  dtype    = nm_interpret_dtype(argc-1-offset, argv+offset+1, stype); // 2-3: dtype

  if (IS_NUMERIC(argv[1+offset]) || TYPE(argv[1+offset]) == T_ARRAY) { // initial value provided (could also be initial capacity, if yale)
    if (stype == S_YALE) {
      init_val     = nm_interpret_initial_capacity(argv[1+offset]);
      init_val_len = 1;
    } else {
      init_val = nm_interpret_initial_value(argv[1+offset], dtype);// 4: initial value / dtype
      if (TYPE(argv[1+offset]) == T_ARRAY) init_val_len = RARRAY_LEN(argv[1+offset]);
      else                                 init_val_len = 1;
    }
  } else {
    if (stype == S_DENSE) { // no need to initialize dense with any kind of default value unless it's an NM_ROBJ matrix
      if (dtype == NM_ROBJ) { // pretend [nil] was passed for ROBJ.
        init_val = ALLOC(VALUE);
        SetFuncs[NM_ROBJ][NM_ROBJ](1, init_val, 0, &QNIL, 0);
        init_val_len = 1;
      } else init_val = NULL;
    } else if (stype == S_YALE) { // if it's a list or compressed, we want to assume default of 0 even if none provided
      init_val = ALLOC(size_t);
      *(size_t*)init_val = 0;
    } else {
      init_val = ALLOC_N(char, nm_sizeof[dtype]);
      //memset(init_val, 0, nm_sizeof[dtype]); // TODO: See if this works instead of the next line (with NM_ROBJ matrix). Cleaner.
      SetFuncs[dtype][NM_BYTE](1, init_val, 0, &ZERO, 0);
    }
  }


  // TODO: Update to allow an array as the initial value.

  if (dtype == NM_NONE) {
    rb_raise(rb_eArgError, "Could not recognize dtype");
    free(init_val);
    free(shape);
    return nm;
  }

  if (stype < S_TYPES) {
    UnwrapNMatrix( nm, nmatrix );

    nmatrix->stype   = stype;
    nmatrix->storage = CreateFuncs[stype](shape, rank, dtype, init_val, init_val_len, nm);

    return nm;
  } else
    rb_raise(rb_eNotImpError, "Unrecognized storage type");

  free(shape);
  free(init_val);
  return nm;
}


static VALUE nm_alloc(VALUE klass) {
  NMATRIX* mat = ALLOC(NMATRIX);
  mat->storage = NULL;
  mat->stype   = S_TYPES;
  return Data_Wrap_Struct(klass, MarkFuncs[mat->stype], nm_delete, mat);
}


// This is the "back-door initializer," for when Ruby needs to create the object in an atypical way.
//
// Note that objects created this way will have NULL storage.
/*static VALUE nm_initialize(VALUE self, VALUE stype, VALUE dtype) {
  NMATRIX* matrix;
  UnwrapNMatrix(self, matrix);

  matrix->stype   = nm_interpret_stype(stype);
  matrix->dtype   = nm_interpret_dtype(1, &dtype, stype);
  matrix->storage = NULL;

  return self;
}*/


static VALUE nm_init_copy(VALUE copy, VALUE original) {
  NMATRIX *lhs, *rhs;

  CheckNMatrixType(original);

  if (copy == original) return copy;

  UnwrapNMatrix( original, rhs );
  UnwrapNMatrix( copy,     lhs );

  lhs->stype = rhs->stype;

  // Copy the storage
  lhs->storage = CastCopyFuncs[rhs->stype](rhs->storage, rhs->storage->dtype);

  return copy;
}


static VALUE nm_init_cast_copy(VALUE copy, VALUE original, VALUE new_dtype_symbol) {
  NMATRIX *lhs, *rhs;
  int8_t new_dtype = nm_dtypesymbol_to_dtype(new_dtype_symbol);
  //fprintf(stderr,"In copy constructor\n");

  CheckNMatrixType(original);

  if (copy == original) return copy;

  UnwrapNMatrix( original, rhs );
  UnwrapNMatrix( copy,     lhs );

  lhs->stype = rhs->stype;

  // Copy the storage
  lhs->storage = CastCopyFuncs[rhs->stype](rhs->storage, new_dtype);

  return copy;
}


static VALUE nm_cast_copy(VALUE self, VALUE new_dtype_symbol) {
  NMATRIX *original, *copy;
  int8_t new_dtype = nm_dtypesymbol_to_dtype(new_dtype_symbol);

  CheckNMatrixType(self);

  UnwrapNMatrix(self, original);

  copy = ALLOC(NMATRIX);
  copy->stype = original->stype;
  copy->storage = CastCopyFuncs[original->stype](original->storage, new_dtype);

  return Data_Wrap_Struct(cNMatrix, MarkFuncs[copy->stype], nm_delete, copy);
}


// Cast a single matrix to a new dtype (unless it's already casted, then just return it). Helper for binary_storage_cast_alloc.
static inline STORAGE* storage_cast_alloc(NMATRIX* matrix, int8_t new_dtype) {
  if (matrix->storage->dtype == new_dtype) return matrix->storage;
  else                                     return CastCopyFuncs[matrix->stype](matrix->storage, new_dtype);
}


// Cast a pair of matrices for a binary operation to a new dtype (which this function determines). Technically, only
// does an actual cast on matrices that are the wrong dtype; otherwise returns a reference to the original. Bear this in
// mind when freeing memory!
static inline STORAGE_PAIR binary_storage_cast_alloc(NMATRIX* left_matrix, NMATRIX* right_matrix) {
  STORAGE_PAIR casted;
  int8_t new_dtype = Upcast[left_matrix->storage->dtype][right_matrix->storage->dtype];

  casted.left  = storage_cast_alloc(left_matrix, new_dtype);
  casted.right = storage_cast_alloc(right_matrix, new_dtype);

  return casted;
}


// Preconditions: NMatrices with SAME stype but possibly different dtype. This function does NOT verify stype.
// Returns: result NMatrix with upcasted dtype and same stype as left and right.
static VALUE multiply_matrix(NMATRIX* left, NMATRIX* right) {
  ///TODO: multiplication for non-dense and/or non-decimal matrices
  size_t*  resulting_shape   = ALLOC_N(size_t, 2);
  NMATRIX* result;
  bool     vector = false;

  // Make sure both of our matrices are of the correct type.
  STORAGE_PAIR casted = binary_storage_cast_alloc(left, right);

  resulting_shape[0] = left->storage->shape[0];
  resulting_shape[1] = right->storage->shape[1];

  // Sometimes we only need to use matrix-vector multiplication (e.g., GEMM versus GEMV). Find out.
  if (resulting_shape[1] == 1) vector = true;

  result = CastedMultiplyFuncs[left->stype](casted, resulting_shape, vector);

  // Free any casted-storage we created for the multiplication.
  // TODO: This is probably an indication that we should use the Ruby garbage collector instead.
  // If we did that, we night not have to re-create these every time, right? Or wrong? Need to do
  // more research.
  if (left->storage != casted.left)   DeleteFuncs[left->stype](casted.left);
  if (right->storage != casted.right) DeleteFuncs[left->stype](casted.right);

  if (result) return Data_Wrap_Struct(cNMatrix, MarkFuncs[result->stype], nm_delete, result);
  return Qnil; // Only if we try to multiply list matrices should we return Qnil.
}


static VALUE multiply_scalar(NMATRIX* left, VALUE scalar) {
  rb_raise(rb_eNotImpError, "matrix-scalar multiplication not implemented yet");
  return Qnil;
}


static VALUE nm_multiply(VALUE left_v, VALUE right_v) {
  NMATRIX *left, *right;

  // left has to be of type NMatrix.
  CheckNMatrixType(left_v);

  UnwrapNMatrix( left_v, left );

  if (IS_NUMERIC(right_v))
    return multiply_scalar(left, right_v);

  else if (TYPE(right_v) == T_ARRAY)
    rb_raise(rb_eNotImpError, "for matrix-vector multiplication, please use an NVector instead of an Array for now");

  //if (RDATA(right_v)->dfree != (RUBY_DATA_FUNC)nm_delete) {
  else if (TYPE(right_v) == T_DATA && RDATA(right_v)->dfree == (RUBY_DATA_FUNC)nm_delete) { // both are matrices
    UnwrapNMatrix( right_v, right );

    if (left->storage->shape[1] != right->storage->shape[0])
      rb_raise(rb_eArgError, "incompatible dimensions");

    if (left->stype != right->stype)
      rb_raise(rb_eNotImpError, "matrices must have same stype");

    return multiply_matrix(left, right);

  } else rb_raise(rb_eTypeError, "expected right operand to be NMatrix, NVector, or single numeric value");

  return Qnil;
}


static VALUE elementwise(NMATRIX* left, NMATRIX* right, char op) {
  ///TODO: multiplication for non-dense and/or non-decimal matrices
  NMATRIX* result;

  // Make sure both of our matrices are of the correct type.
  STORAGE_PAIR casted = binary_storage_cast_alloc(left, right);

  result = CastedElementwiseFuncs[left->stype](casted, op);

  // Free up temporary casted matrices
  if (left->storage != casted.left)   DeleteFuncs[left->stype](casted.left);
  if (right->storage != casted.right) DeleteFuncs[left->stype](casted.right);

  if (result) return Data_Wrap_Struct(cNMatrix, MarkFuncs[result->stype], nm_delete, result);
  return Qnil; // Only if we try to multiply list matrices should we return Qnil.
}


static VALUE nm_add(VALUE leftv, VALUE rightv) {
  NMATRIX *left, *right;

  // left has to be of type NMatrix.
  CheckNMatrixType(leftv);
  CheckNMatrixType(rightv);

  UnwrapNMatrix(rightv, right);
  UnwrapNMatrix(leftv, left);





  return Qnil;
}


// Borrowed this function from NArray. Handles 'each' iteration on a dense matrix.
//
// Additionally, handles separately matrices containing VALUEs and matrices containing
// other types of data.
static VALUE nm_dense_each(VALUE nmatrix) {
  DENSE_STORAGE* s = NM_STORAGE(nmatrix);
  VALUE v;
  size_t i;

  void (*copy)();

  if (NM_DTYPE(nmatrix) == NM_ROBJ) {

    // matrix of Ruby objects -- yield directly
    for (i = 0; i < count_dense_storage_elements(s); ++i)
      rb_yield( *((VALUE*)((char*)(s->elements) + i*nm_sizeof[NM_DTYPE(nmatrix)])) );

  } else {
    // We're going to copy the matrix element into a Ruby VALUE and then operate on it.
    copy = SetFuncs[NM_ROBJ][NM_DTYPE(nmatrix)];

    for (i = 0; i < count_dense_storage_elements(s); ++i) {
      (*copy)(1, &v, 0, (char*)(s->elements) + i*nm_sizeof[NM_DTYPE(nmatrix)], 0);
      rb_yield(v); // yield to the copy we made
    }
  }

  return nmatrix;
}


// iterate through contents of an NMatrix using `each`
static VALUE nm_each(VALUE nmatrix) {
  volatile VALUE nm = nmatrix; // not sure why we do this, but it gets done in ruby's array.c.

  switch(NM_STYPE(nm)) {
  case S_DENSE:
    return nm_dense_each(nm);
  default:
    rb_raise(rb_eNotImpError, "only dense each works right now");
  }
}


// Does not create storage, but does destroy it.
NMATRIX* nm_create(int8_t stype, void* storage) {
  NMATRIX* mat = ALLOC(NMATRIX);

  mat->stype   = stype;
  mat->storage = storage;

  return mat;
}


size_t* convert_coords(size_t rank, VALUE* c, VALUE self) {
  size_t r;
  size_t* coords = ALLOC_N(size_t,rank);

  for (r = 0; r < rank; ++r) {
    coords[r] = FIX2UINT(c[r]);
    if (coords[r] >= NM_SHAPE(self,r)) rb_raise(rb_eArgError, "out of range");
  }

  return coords;
}


VALUE nm_mref(int argc, VALUE* argv, VALUE self) {
  VALUE v;

  if (NM_RANK(self) == (size_t)(argc)) {

    SetFuncs[NM_ROBJ][NM_DTYPE(self)](1, &v, 0,
      RefFuncs[NM_STYPE(self)](NM_STORAGE(self),
                               convert_coords((size_t)(argc), argv, self)
                              ), 0);
    return v;

  } else if (NM_RANK(self) < (size_t)(argc)) {
    rb_raise(rb_eArgError, "Coordinates given exceed matrix rank");
  } else {
    rb_raise(rb_eNotImpError, "Slicing not supported yet");
  }
  return Qnil;
}


VALUE nm_mset(int argc, VALUE* argv, VALUE self) {
  size_t rank = argc - 1; // last arg is the value

  if (argc <= 1) {
    rb_raise(rb_eArgError, "Expected coordinates and r-value");

  } else if (NM_RANK(self) == rank) {
    return (*(InsFuncs[NM_STYPE(self)]))( NM_STORAGE(self),
                                         convert_coords(rank, argv, self),
                                         argv[rank] );

  } else if (NM_RANK(self) < rank) {
    rb_raise(rb_eArgError, "Coordinates given exceed matrix rank");
  } else {
    rb_raise(rb_eNotImpError, "Slicing not supported yet");
  }
  return Qnil;
}


VALUE nm_rank(VALUE self) {
  VALUE ret;
  SetFuncs[NM_ROBJ][NM_INT64]( 1, &ret, 0, &(NM_STORAGE(self)->rank), 0 );
  return ret;
}


VALUE nm_shape(VALUE self) {
  STORAGE* s   = NM_STORAGE(self);

  // Copy elements into a VALUE array and then use those to create a Ruby array with rb_ary_new4.
  VALUE* shape = ALLOCA_N(VALUE, s->rank);
  SetFuncs[NM_ROBJ][NM_SIZE_T]( s->rank, shape, sizeof(VALUE), s->shape, sizeof(size_t));

  return rb_ary_new4(s->rank, shape);
}


static VALUE nm_stype(VALUE self) {
  ID stype = rb_intern(nm_stypestring[NM_STYPE(self)]);
  return ID2SYM(stype);
}


static VALUE nm_dtype(VALUE self) {
  ID dtype = rb_intern(nm_dtypestring[NM_DTYPE(self)]);
  return ID2SYM(dtype);
}


// Interprets cblas argument which could be any of false/:no_transpose, :transpose, or :complex_conjugate,
// into an enum recognized by cblas.
//
// Called by nm_cblas_gemm -- basically inline.
static char gemm_op_sym(VALUE op) {
  if (op == false || rb_to_id(op) == nm_id_no_transpose) return CblasNoTrans;
  else if (rb_to_id(op) == nm_id_transpose) return CblasTrans;
  else if (rb_to_id(op) == nm_id_complex_conjugate) return CblasConjTrans;
  else rb_raise(rb_eArgError, "Expected false, :transpose, or :complex_conjugate");
  return CblasNoTrans;
}


// Directly call cblas_xgemm as quickly as possible.
//
//     C: = alpha*op(A)*op(B) + beta*C
//
// where op(X) is one of:
//
//     op(X) = X   or   op(X) = X**T
//
// == Arguments
// See: http://www.netlib.org/blas/dgemm.f
//
// You probably don't want to call this function. Instead of __cblas_gemm__, why don't you try cblas_gemm, which
// is more flexible with its arguments?
//
// This function does almost no type checking. Seriously, be really careful when you call it! There's no exception
// handling!
static VALUE nm_cblas_gemm(VALUE self,
                           VALUE trans_a, VALUE trans_b,
                           VALUE m, VALUE n, VALUE k,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE b, VALUE ldb,
                           VALUE beta,
                           VALUE c, VALUE ldc)
{
  struct cblas_param_t p = cblas_params_for_multiply(((DENSE_STORAGE*)(NM_STORAGE(a))), ((DENSE_STORAGE*)(NM_STORAGE(b))), ((DENSE_STORAGE*)(NM_STORAGE(c))), false);
  p.M = FIX2INT(m);
  p.N = FIX2INT(n);
  p.K = FIX2INT(k);

  p.lda = FIX2INT(lda);
  p.ldb = FIX2INT(ldb);
  p.ldc = FIX2INT(ldc);

  switch(NM_DTYPE(c)) {
  case NM_FLOAT32:
  case NM_FLOAT64:
    p.alpha.d[0] = NUM2DBL(alpha);
    p.beta.d[0] = NUM2DBL(beta);
    break;

  case NM_COMPLEX64:
    p.alpha.c[0].r = REAL2DBL(alpha);
    p.alpha.c[0].i = IMAG2DBL(alpha);
    p.beta.c[0].r = REAL2DBL(beta);
    p.beta.c[0].i = IMAG2DBL(beta);
    break;

  case NM_COMPLEX128:
    p.alpha.z.r = REAL2DBL(alpha);
    p.alpha.z.i = IMAG2DBL(alpha);
    p.beta.z.r = REAL2DBL(beta);
    p.beta.z.i = IMAG2DBL(beta);
    break;

  case NM_BYTE:
    p.alpha.b[0] = FIX2INT(alpha);
    p.beta.b[0] = FIX2INT(beta);
    break;

  case NM_INT8:
  case NM_INT16:
  case NM_INT32:
  case NM_INT64:
    p.alpha.i[0] = FIX2INT(alpha);
    p.beta.i[0]  = FIX2INT(beta);
    break;

  case NM_RATIONAL32:
    p.alpha.r[0].n = NUMER2INT(alpha);
    p.alpha.r[0].d = DENOM2INT(alpha);
    p.beta.r[0].n = NUMER2INT(beta);
    p.beta.r[0].d = DENOM2INT(beta);
    break;

  case NM_RATIONAL64:
    p.alpha.ra[0].n = NUMER2INT(alpha);
    p.alpha.ra[0].d = DENOM2INT(alpha);
    p.beta.ra[0].n = NUMER2INT(beta);
    p.beta.ra[0].d = DENOM2INT(beta);
    break;

  case NM_RATIONAL128:
    p.alpha.rat.n = NUMER2INT(alpha);
    p.alpha.rat.d = DENOM2INT(alpha);
    p.beta.rat.n = NUMER2INT(beta);
    p.beta.rat.d = DENOM2INT(beta);
    break;

  case NM_ROBJ:
    p.alpha.v[0] = alpha;
    p.beta.v[0]  = beta;
    break;

  default:
    rb_raise(nm_eDataTypeError, "unexpected dtype");

  }

  /* fprintf(stderr, "cblas_gemm: %d %d %d %d %d %f %d %d %f %d\n", trans_a_, trans_b_,
         m_, n_, k_, alpha_, lda_, ldb_, beta_, ldc_); */

  GemmFuncs[NM_DTYPE(c)](CblasRowMajor, gemm_op_sym(trans_a), gemm_op_sym(trans_b), p);

  return Qtrue;
}


static VALUE nm_cblas_gemv(VALUE self,
                           VALUE trans_a,
                           VALUE m, VALUE n,
                           VALUE alpha,
                           VALUE a, VALUE lda,
                           VALUE x, VALUE incx,
                           VALUE beta,
                           VALUE y, VALUE incy) {

  struct cblas_param_t p;
  p.M = FIX2INT(m);
  p.N = FIX2INT(n);
  p.A = ((DENSE_STORAGE*)(NM_STORAGE(a)))->elements;
  p.B = ((DENSE_STORAGE*)(NM_STORAGE(x)))->elements;
  p.C = ((DENSE_STORAGE*)(NM_STORAGE(y)))->elements;
  p.lda = FIX2INT(lda);
  p.ldb = FIX2INT(incx);
  p.ldc = FIX2INT(incy);

  switch(NM_DTYPE(y)) {
  case NM_FLOAT32:
  case NM_FLOAT64:
    p.alpha.d[0] = REAL2DBL(alpha);
    p.beta.d[0]  = REAL2DBL(beta);
    break;
  case NM_COMPLEX64:
    p.alpha.c[0].r = REAL2DBL(alpha);
    p.alpha.c[0].i = IMAG2DBL(alpha);
    p.beta.c[0].r = REAL2DBL(beta);
    p.beta.c[0].i = IMAG2DBL(beta);
    break;
  case NM_COMPLEX128:
    p.alpha.z.r = REAL2DBL(alpha);
    p.alpha.z.i = IMAG2DBL(alpha);
    p.beta.z.r = REAL2DBL(beta);
    p.beta.z.i = IMAG2DBL(beta);
    break;
  }

  /* fprintf(stderr, "cblas_gemm: %d %d %d %d %d %f %d %d %f %d\n", trans_a_, trans_b_,
         m_, n_, k_, alpha_, lda_, ldb_, beta_, ldc_); */

  GemvFuncs[NM_DTYPE(y)](CblasRowMajor, gemm_op_sym(trans_a), p);

  return Qtrue;
}


static VALUE nm_capacity(VALUE self) {
  VALUE cap;
  size_t r, tmp=1;

  switch(NM_STYPE(self)) {
  case S_YALE:

    cap = UINT2NUM(((YALE_STORAGE*)(NM_STORAGE(self)))->capacity);
    break;

  case S_DENSE:

    for (r = 0; r < NM_STORAGE(self)->rank; ++r)
      tmp *= NM_STORAGE(self)->shape[r];

    SetFuncs[NM_ROBJ][NM_INT64](1, &cap, 0, tmp, 0);
    break;

  case S_LIST:
  default:
    rb_raise(rb_eNotImpError, "TODO: implement capacity/size on other storage types");
  }

  return cap;
}



static VALUE nm_yale_size(VALUE self) {
  VALUE sz;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  SetFuncs[NM_ROBJ][s->index_dtype](1, &sz, 0, (YALE_SIZE_PTR((s), nm_sizeof[s->index_dtype])), 0);
  return sz;
}


static VALUE nm_yale_a(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*sz);

  SetFuncs[NM_ROBJ][s->dtype](sz, vals, nm_sizeof[NM_ROBJ], s->a, nm_sizeof[s->dtype]);
  ary = rb_ary_new4(sz, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}


static VALUE nm_yale_d(VALUE self) {
  y_size_t sz;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*s->shape[0]);

  SetFuncs[NM_ROBJ][s->dtype](s->shape[0], vals, nm_sizeof[NM_ROBJ], s->a, nm_sizeof[s->dtype]);
  ary = rb_ary_new4(s->shape[0], vals);

  return ary;
}


static VALUE nm_yale_lu(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*(s->capacity - s->shape[0]));

  SetFuncs[NM_ROBJ][s->dtype](sz - s->shape[0] - 1, vals, nm_sizeof[NM_ROBJ], (char*)(s->a) + (s->shape[0] + 1)*nm_sizeof[s->dtype], nm_sizeof[s->dtype]);
  ary = rb_ary_new4(sz - s->shape[0] - 1, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}


// Only works for i8/f64
static VALUE nm_yale_print_vectors(VALUE self) {
  if (NM_STYPE(self) != S_YALE || NM_DTYPE(self) != NM_FLOAT64) rb_raise(nm_eDataTypeError, "must be yale float64 matrix");

  print_vectors((YALE_STORAGE*)(NM_STORAGE(self)));

  return Qnil;
}


static VALUE nm_yale_ia(VALUE self) {
  y_size_t sz;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*(s->shape[0]+1));

  SetFuncs[NM_ROBJ][s->index_dtype](s->shape[0]+1, vals, nm_sizeof[NM_ROBJ], s->ija, nm_sizeof[s->index_dtype]);
  ary = rb_ary_new4(s->shape[0]+1, vals);

  return ary;
}


static VALUE nm_yale_ja(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*(s->capacity - s->shape[0]));

  SetFuncs[NM_ROBJ][s->index_dtype](sz - s->shape[0] - 1, vals, nm_sizeof[NM_ROBJ], (char*)(s->ija) + (s->shape[0] + 1)*nm_sizeof[s->index_dtype], nm_sizeof[s->index_dtype]);
  ary = rb_ary_new4(sz - s->shape[0] - 1, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}


static VALUE nm_yale_ija(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = (YALE_STORAGE*)NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(nm_eStorageTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*s->capacity);

  SetFuncs[NM_ROBJ][s->index_dtype](sz, vals, nm_sizeof[NM_ROBJ], s->ija, nm_sizeof[s->index_dtype]);
  ary = rb_ary_new4(sz, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}


static VALUE nm_transpose_new(VALUE self) {
  NMATRIX *self_m, *result, *result2;
  size_t sz;
  size_t* shape   = ALLOC_N(size_t, 2);
  YALE_PARAM A, B;
#ifdef BENCHMARK
  double t1, t2;
#endif

  UnwrapNMatrix( self, self_m );

  // switch the dimensions
  shape[1] = self_m->storage->shape[0];
  shape[0] = self_m->storage->shape[1];

  switch(self_m->stype) {
  case S_DENSE:
    //result   = nm_create(S_DENSE, create_dense_storage(left->storage->dtype, shape, 2));

    // call CBLAS xgemm (type-specific general matrix multiplication)
    // good explanation: http://www.umbc.edu/hpcf/resources-tara/how-to-BLAS.html
    /*GemmFuncs[left->storage->dtype](
    //cblas_sgemm(
          CblasRowMajor,
          CblasNoTrans,
          CblasNoTrans,
          shape[0],                // M = number of rows in left
          shape[1],                // N = number of columns in right and result
          left->storage->shape[1], // K = number of columns in left
          1.0,
          ((DENSE_STORAGE*)(left->storage))->elements,
              left->storage->shape[1], // * nm_sizeof[left->dtype],
          ((DENSE_STORAGE*)(right->storage))->elements,
              right->storage->shape[1], // * nm_sizeof[right->dtype],
          0.0,
          ((DENSE_STORAGE*)(result->storage))->elements,
              shape[1] // * nm_sizeof[result->dtype]
          ); */
    rb_raise(rb_eNotImpError, "need dense transpose function");
    break;
  case S_YALE:
    YaleGetSize(sz, (YALE_STORAGE*)(self_m->storage)); // size of new matrix is going to be size of old matrix
    result = nm_create(S_YALE, create_yale_storage(self_m->storage->dtype, shape, 2, sz));

    // TODO: Do we really need to initialize the whole thing? Or just the A portion?
    init_yale_storage((YALE_STORAGE*)(result->storage));

    result2 = nm_create(S_YALE, create_yale_storage(self_m->storage->dtype, shape, 2, sz));
    init_yale_storage((YALE_STORAGE*)(result2->storage));

    A.ia = A.ja = ((YALE_STORAGE*)(self_m->storage))->ija;
    B.ia = B.ja = ((YALE_STORAGE*)(result->storage))->ija;
    A.a  = ((YALE_STORAGE*)(self_m->storage))->a;
    B.a  = ((YALE_STORAGE*)(result->storage))->a;
    A.diag = true;

#ifdef BENCHMARK
    t1 = get_time();
#endif

    // call the appropriate function pointer
    SparseTransposeFuncs[ self_m->storage->dtype ][ ((YALE_STORAGE*)(self_m->storage))->index_dtype ](shape[0], shape[1], A, B, true);
#ifdef BENCHMARK
    t1 = get_time() - t1;
/*
    t2 = get_time();
    transp(
          shape[0],
          shape[1],
          ((YALE_STORAGE*)(self_m->storage))->ija,
          ((YALE_STORAGE*)(self_m->storage))->ija,
          true,
          ((YALE_STORAGE*)(self_m->storage))->a,
          ((YALE_STORAGE*)(result2->storage))->ija,
          ((YALE_STORAGE*)(result2->storage))->ija,
          ((YALE_STORAGE*)(result2->storage))->a,
          true, // move
          ((YALE_STORAGE*)(self_m->storage))->index_dtype,
          self_m->storage->dtype
          );

    t2 = get_time() - t2;
    fprintf(stderr, "t1: %f\nt2: %f\n", t1, t2);
*/
#endif

    break;
  default:
    rb_raise(rb_eNotImpError, "transpose for this type not implemented yet");
  }

  return Data_Wrap_Struct(cNMatrix, MarkFuncs[result->stype], nm_delete, result);
}

//static VALUE nm_transpose_auto(VALUE self) {
//
//}


void Init_nmatrix() {
    /* Require Complex class */
    //rb_require("complex");
    //cComplex = rb_const_get( rb_cObject, rb_intern("Complex") );

    /* Define NMatrix class */
    cNMatrix = rb_define_class("NMatrix", rb_cObject);

    /* class methods */
    rb_define_singleton_method(cNMatrix, "__cblas_gemm__", nm_cblas_gemm, 13);
    rb_define_singleton_method(cNMatrix, "__cblas_gemv__", nm_cblas_gemv, 11);

    rb_define_alloc_func(cNMatrix, nm_alloc);
    rb_define_method(cNMatrix, "initialize", nm_init, -1);
    // rb_define_singleton_method(cNMatrix, "new", nm_init, -1);


    rb_define_method(cNMatrix, "initialize_copy", nm_init_copy, 1);
    rb_define_method(cNMatrix, "initialize_cast_copy", nm_init_cast_copy, 2);
    rb_define_method(cNMatrix, "as_dtype", nm_cast_copy, 1);

    /* methods */
    rb_define_method(cNMatrix, "dtype", nm_dtype, 0);
    rb_define_method(cNMatrix, "stype", nm_stype, 0);

    rb_define_method(cNMatrix, "[]", nm_mref, -1);
    rb_define_method(cNMatrix, "[]=", nm_mset, -1);
    rb_define_method(cNMatrix, "rank", nm_rank, 0);
    rb_define_alias(cNMatrix, "dim", "rank");
    rb_define_method(cNMatrix, "shape", nm_shape, 0);
    rb_define_method(cNMatrix, "transpose", nm_transpose_new, 0);
    //rb_define_method(cNMatrix, "transpose!", nm_transpose_auto, 0);

    rb_define_method(cNMatrix, "each", nm_each, 0);

    rb_define_method(cNMatrix, "*", nm_multiply, 1);
    //rb_define_method(cNMatrix, "+", nm_add, 1);
    rb_define_alias(cNMatrix, "dot", "*");


    rb_define_method(cNMatrix, "capacity", nm_capacity, 0);

    rb_define_method(cNMatrix, "__yale_print__", nm_yale_print_vectors, 0);
    rb_define_method(cNMatrix, "__yale_ija__", nm_yale_ija, 0);
    rb_define_method(cNMatrix, "__yale_a__", nm_yale_a, 0);
    rb_define_method(cNMatrix, "__yale_size__", nm_yale_size, 0);
    rb_define_method(cNMatrix, "__yale_ia__", nm_yale_ia, 0);
    rb_define_method(cNMatrix, "__yale_ja__", nm_yale_ja, 0);
    rb_define_method(cNMatrix, "__yale_d__", nm_yale_d, 0);
    rb_define_method(cNMatrix, "__yale_lu__", nm_yale_lu, 0);
    rb_define_const(cNMatrix, "YALE_GROWTH_CONSTANT", rb_float_new(YALE_GROWTH_CONSTANT));


    cNVector = rb_define_class("NVector", cNMatrix);

    // Special exceptions
    nm_eDataTypeError    = rb_define_class("DataTypeError", rb_eStandardError);
    nm_eStorageTypeError = rb_define_class("StorageTypeError", rb_eStandardError);

    nm_id_real  = rb_intern("real");
    nm_id_imag  = rb_intern("imag");
    nm_id_numer = rb_intern("numerator");
    nm_id_denom = rb_intern("denominator");
    nm_id_mult  = rb_intern("*");
    nm_id_add   = rb_intern("+");
    nm_id_multeq= rb_intern("*=");

    nm_id_transpose = rb_intern("transpose");
    nm_id_no_transpose = rb_intern("no_transpose");
    nm_id_complex_conjugate = rb_intern("complex_conjugate");

    nm_id_dense = rb_intern("dense");
    nm_id_list = rb_intern("list");

}

#endif