/* Loosely based upon NArray */
#ifndef NMATRIX_C
# define NMATRIX_C

#include <ruby.h>

#include "nmatrix.h"

VALUE cNMatrix;

ID nm_id_real, nm_id_imag;
ID nm_id_numer, nm_id_denom;
ID nm_id_transpose, nm_id_no_transpose, nm_id_complex_conjugate; // cblas
ID nm_id_list, nm_id_dense;

#include "dtypes.c"


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

nm_gemm_t GemmFuncs = { // by NM_TYPES
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  cblas_sgemm,
  cblas_dgemm,
  cblas_cgemm,
  cblas_zgemm,
  NULL,
  NULL,
  NULL,
  NULL
};


static void nm_delete(NMATRIX* mat) {
  DeleteFuncs[mat->stype](mat->storage);
}



VALUE nm_dense_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, VALUE self) {
  DENSE_STORAGE* s = create_dense_storage(dtype, shape, rank);
  NMATRIX* matrix  = nm_create(S_DENSE, s);
  size_t i, n;

  // User provides single initialization value sometimes. This accommodates List, which
  // really doesn't make sense without a user-specified default (initialization) value
  // (but can be done).
  //
  // If the user instead wants to provide an array or something, that's done through a
  // different mechanism. He or she should then NOT provide a default value for dense, only
  // because that requires an unnecessary pass through the array of elements.
  if (init_val) {
    n = count_dense_storage_elements(s);
    for (i = 0; i < n; ++i)
      memcpy((char*)(s->elements) + i*nm_sizeof[dtype], init_val, nm_sizeof[dtype]);
    free(init_val);
  }
  return Data_Wrap_Struct(self, NULL, nm_delete, matrix);
}

VALUE nm_list_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, VALUE self) {
  NMATRIX* matrix = nm_create(S_LIST, create_list_storage(dtype, shape, rank, init_val));
  return Data_Wrap_Struct(self, NULL, nm_delete, matrix);
}


VALUE nm_yale_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, VALUE self) {
  NMATRIX* matrix;
  YALE_STORAGE* s = create_yale_storage(dtype, shape, rank, *(size_t*)init_val);
                    init_yale_storage(s);

  if (!s) rb_raise(rb_eNoMemError, "Yale allocation failed");

  matrix = nm_create(S_YALE, s);

  free(init_val);
  return Data_Wrap_Struct(self, NULL, nm_delete, matrix);
}


nm_create_t CreateFuncs = {
  nm_dense_new,
  nm_list_new,
  nm_yale_new
};


nm_copy_s_t CopyFuncs = {
  copy_dense_storage,
  copy_list_storage,
  copy_yale_storage
};


VALUE nm_dense_get(STORAGE* s, size_t* coords) {
  VALUE v;
  SetFuncs[NM_ROBJ][s->dtype](1, &v, 0, dense_storage_get((DENSE_STORAGE*)s, coords), 0);
  return v;
}

VALUE nm_list_get(STORAGE* s, size_t* coords) {
  VALUE v;
  SetFuncs[NM_ROBJ][s->dtype](1, &v, 0, list_storage_get((LIST_STORAGE*)s, coords), 0);
  return v;
}


VALUE nm_yale_get(STORAGE* s, size_t* coords) {
  VALUE v;
  SetFuncs[NM_ROBJ][s->dtype](1, &v, 0, yale_storage_ref((YALE_STORAGE*)s, coords), 0);
  return v;
}


nm_stype_ref_t RefFuncs = {
  nm_dense_get,
  nm_list_get,
  nm_yale_get
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


VALUE nm_list_set(STORAGE* s, size_t* coords, VALUE val) {
  void *v = malloc(nm_sizeof[s->dtype]), *rm;
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



nm_stype_ref_t InsFuncs = {
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
    shape = malloc( sizeof(size_t) * (*rank) );
    for (i = 0; i < *rank; ++i)
      shape[i]  = (size_t)(FIX2UINT(RARRAY_PTR(arg)[i]));
  } else if (FIXNUM_P(arg)) {
    *rank = 2;
    shape = malloc( sizeof(size_t) * (*rank) );
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
  if (argc == 2) {
    if (SYMBOL_P(argv[1])) return nm_dtypesymbol_to_dtype(argv[1]);
    else if (IS_STRING(argv[1])) return nm_dtypestring_to_dtype(StringValue(argv[1]));
    else if (stype == S_YALE) rb_raise(rb_eArgError, "yale requires dtype");
    else return nm_guess_dtype(argv[0]);
  } else if (argc == 1) {
    if (SYMBOL_P(argv[0])) return nm_dtypesymbol_to_dtype(argv[0]);
    else if (IS_STRING(argv[0])) return nm_dtypestring_to_dtype(StringValue(argv[0]));
    else if (stype == S_YALE) rb_raise(rb_eArgError, "yale requires dtype");
    else return nm_guess_dtype(argv[0]);
  } else rb_raise(rb_eArgError, "Need an initial value or a dtype");

  return NM_NONE;
}

int8_t nm_interpret_stype(VALUE arg) {
  if (SYMBOL_P(arg)) return nm_stypesymbol_to_stype(arg);
  else if (IS_STRING(arg)) return nm_stypestring_to_stype(StringValue(arg));
  else rb_raise(rb_eArgError, "Expected storage type");
  return S_DENSE;
}

void* nm_interpret_initial_value(VALUE arg, int8_t dtype) {
  void* init_val = malloc(nm_sizeof[dtype]);
  SetFuncs[dtype][NM_ROBJ](1, init_val, 0, &arg, 0);
  return init_val;
}


size_t* nm_interpret_initial_capacity(VALUE arg) {
  size_t* init_cap = malloc(sizeof(size_t));
  *init_cap = FIX2UINT(arg);
  return init_cap;
}


VALUE nm_init(int argc, VALUE* argv, VALUE self) {
  char    ZERO = 0;
  int8_t  dtype, stype, offset = 0;
  size_t  rank;
  size_t* shape;
  void*   init_val = NULL;

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

  if (IS_NUMERIC(argv[1+offset])) { // initial value provided (could also be initial capacity, if yale)
    if (stype == S_YALE)      init_val = nm_interpret_initial_capacity(argv[1+offset]);
    else                      init_val = nm_interpret_initial_value(argv[1+offset], dtype);// 4: initial value / dtype
  } else if (stype == S_DENSE) {
    init_val = NULL; // no need to initialize dense with any kind of default value.
  } else { // if it's a list or compressed, we want to assume default of 0 even if none provided
    if (stype == S_YALE) {
      init_val = malloc(sizeof(size_t));
      *(size_t*)init_val = 0;
    } else {
      init_val = malloc(nm_sizeof[dtype]);
      SetFuncs[dtype][NM_BYTE](1, init_val, 0, &ZERO, 0);
    }
  }


  // TODO: Update to allow an array as the initial value.

  if (dtype == NM_NONE) {
    rb_raise(rb_eArgError, "Could not recognize dtype");
    free(init_val);
    free(shape);
    return Qnil;
  }

  if (stype < S_TYPES)
    return CreateFuncs[stype](shape, rank, dtype, init_val, self);
  else
    rb_raise(rb_eNotImpError, "Unrecognized storage type");

  free(shape);
  free(init_val);
  return Qnil;

}


static VALUE nm_alloc(VALUE klass) {
  NMATRIX* mat = ALLOC(NMATRIX);
  mat->storage = NULL;
  return Data_Wrap_Struct(klass, 0, nm_delete, mat);
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

  //fprintf(stderr,"In copy constructor\n");

  if (copy == original) return copy;

  if (TYPE(original) != T_DATA || RDATA(original)->dfree != (RUBY_DATA_FUNC)nm_delete)
    rb_raise(rb_eTypeError, "wrong argument type");

  UnwrapNMatrix( original, rhs );
  UnwrapNMatrix( copy,     lhs );

  lhs->stype = rhs->stype;
  //lhs->dtype = rhs->dtype;

  // Copy the storage
  lhs->storage = CopyFuncs[rhs->stype](rhs->storage, nm_sizeof[rhs->storage->dtype]);

  return copy;
}


static VALUE nm_multiply_matrix(NMATRIX* left, NMATRIX* right) {
  ///TODO: multiplication for non-dense and/or non-decimal matrices
  size_t* shape   = malloc(sizeof(size_t)*2);
  NMATRIX* result;

  shape[0] = left->storage->shape[0];
  shape[1] = right->storage->shape[1];

  result   = nm_create(S_DENSE, create_dense_storage(left->storage->dtype, shape, 2));

  //fprintf(stderr, "M=%d, N=%d, K=%d\n", shape[0], shape[1], left->storage->shape[1]);

  // call CBLAS xgemm (type-specific general matrix multiplication)
  // good explanation: http://www.umbc.edu/hpcf/resources-tara/how-to-BLAS.html
  GemmFuncs[left->storage->dtype](
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
        );

  return Data_Wrap_Struct(cNMatrix, 0, nm_delete, result);
}


static VALUE nm_multiply_scalar(NMATRIX* left, VALUE scalar) {
  return Qnil;
}



static VALUE nm_multiply(VALUE left_v, VALUE right_v) {
  NMATRIX *left, *right;

  // left has to be of type NMatrix.
  if (TYPE(left_v) != T_DATA || RDATA(left_v)->dfree != (RUBY_DATA_FUNC)nm_delete)
    rb_raise(rb_eTypeError, "wrong argument type");

  UnwrapNMatrix( left_v, left );

  //if (RDATA(right_v)->dfree != (RUBY_DATA_FUNC)nm_delete) {
  UnwrapNMatrix( right_v, right );

  if (left->storage->shape[1] != right->storage->shape[0])
    rb_raise(rb_eArgError, "incompatible dimensions");

  if (left->stype != S_DENSE && right->stype != S_DENSE)
    rb_raise(rb_eNotImpError, "dense matrices expected");

  if (left->storage->dtype != right->storage->dtype)
    rb_raise(rb_eNotImpError, "dtype mismatch");

  return nm_multiply_matrix(left, right);
  //} else {
  //  rb_raise(rb_eNotImpError, "scalar multiplication not supported yet");
  //  return Qnil;
 // }
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
  if (NM_RANK(self) == (size_t)(argc)) {
    return (*(RefFuncs[NM_STYPE(self)]))( NM_STORAGE(self),
                                          convert_coords((size_t)(argc), argv, self) );

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
  SetFuncs[NM_ROBJ][NM_INT64]( s->rank, shape, sizeof(VALUE), s->shape, sizeof(size_t));

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
                           VALUE c, VALUE ldc) {

  int m_ = FIX2INT(m), n_ = FIX2INT(n), k_ = FIX2INT(k), lda_ = FIX2INT(lda), ldb_ = FIX2INT(ldb), ldc_ = FIX2INT(ldc);
  char trans_a_ = gemm_op_sym(trans_a), trans_b_ = gemm_op_sym(trans_b);
  double alpha_ = NUM2DBL(alpha), beta_ = NUM2DBL(beta);

  /* fprintf(stderr, "cblas_gemm: %d %d %d %d %d %f %d %d %f %d\n", trans_a_, trans_b_,
         m_, n_, k_, alpha_, lda_, ldb_, beta_, ldc_); */

  GemmFuncs[NM_DTYPE(c)](CblasRowMajor, trans_a_, trans_b_, m_, n_, k_,
        alpha_,  ((DENSE_STORAGE*)(NM_STORAGE(a)))->elements, lda_,
                 ((DENSE_STORAGE*)(NM_STORAGE(b)))->elements, ldb_,
        beta_,   ((DENSE_STORAGE*)(NM_STORAGE(c)))->elements, ldc_);

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
  YALE_STORAGE* s = NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(rb_eTypeError, "wrong storage type");

  SetFuncs[NM_ROBJ][s->index_dtype](1, &sz, 0, (YALE_SIZE_PTR((s), nm_sizeof[s->index_dtype])), 0);
  return sz;
}


static VALUE nm_yale_a(VALUE self) {
  y_size_t sz, i;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(rb_eTypeError, "wrong storage type");

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
  YALE_STORAGE* s = NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(rb_eTypeError, "wrong storage type");

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
  YALE_STORAGE* s = NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(rb_eTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*(s->capacity - s->shape[0]));

  SetFuncs[NM_ROBJ][s->dtype](sz - s->shape[0] - 1, vals, nm_sizeof[NM_ROBJ], (char*)(s->a) + (s->shape[0] + 1)*nm_sizeof[s->dtype], nm_sizeof[s->dtype]);
  ary = rb_ary_new4(sz - s->shape[0] - 1, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}


static VALUE nm_yale_print_vectors(VALUE self) {
  if (NM_STYPE(self) != S_YALE) rb_raise(rb_eTypeError, "must be yale matrix");

  print_vectors(NM_STORAGE(self));

  return Qnil;
}


static VALUE nm_yale_ia(VALUE self) {
  y_size_t sz;
  void* vals;
  VALUE ary;
  YALE_STORAGE* s = NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(rb_eTypeError, "wrong storage type");

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
  YALE_STORAGE* s = NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(rb_eTypeError, "wrong storage type");

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
  YALE_STORAGE* s = NM_STORAGE(self);

  if (NM_STYPE(self) != S_YALE) rb_raise(rb_eTypeError, "wrong storage type");

  YaleGetSize(sz, s);
  vals = ALLOC_N(char, nm_sizeof[NM_ROBJ]*s->capacity);

  SetFuncs[NM_ROBJ][s->index_dtype](sz, vals, nm_sizeof[NM_ROBJ], s->ija, nm_sizeof[s->index_dtype]);
  ary = rb_ary_new4(sz, vals);

  for (i = sz; i < s->capacity; ++i)
    rb_ary_push(ary, Qnil);

  return ary;
}


void Init_nmatrix() {
    /* Require Complex class */
    //rb_require("complex");
    //cComplex = rb_const_get( rb_cObject, rb_intern("Complex") );

    /* Define NMatrix class */
    cNMatrix = rb_define_class("NMatrix", rb_cObject);
    rb_define_alloc_func(cNMatrix, nm_alloc);

    /* class methods */
    rb_define_singleton_method(cNMatrix, "__cblas_gemm__", nm_cblas_gemm, 13);
    rb_define_method(cNMatrix, "initialize", nm_init, 2);
    rb_define_singleton_method(cNMatrix, "new", nm_init, -1);


    rb_define_method(cNMatrix, "initialize_copy", nm_init_copy, 1);

    /* methods */
    rb_define_method(cNMatrix, "dtype", nm_dtype, 0);
    rb_define_method(cNMatrix, "stype", nm_stype, 0);
    rb_define_method(cNMatrix, "[]", nm_mref, -1);
    rb_define_method(cNMatrix, "[]=", nm_mset, -1);
    rb_define_method(cNMatrix, "rank", nm_rank, 0);
    rb_define_alias(cNMatrix, "dim", "rank");
    rb_define_method(cNMatrix, "shape", nm_shape, 0);

    rb_define_method(cNMatrix, "*", nm_multiply, 1);
    rb_define_alias(cNMatrix, "multiply", "*");


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

    nm_id_real  = rb_intern("real");
    nm_id_imag  = rb_intern("imag");
    nm_id_numer = rb_intern("numerator");
    nm_id_denom = rb_intern("denominator");

    nm_id_transpose = rb_intern("transpose");
    nm_id_no_transpose = rb_intern("no_transpose");
    nm_id_complex_conjugate = rb_intern("complex_conjugate");

    nm_id_dense = rb_intern("dense");
    nm_id_list = rb_intern("list");

}

#endif