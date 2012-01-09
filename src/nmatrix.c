/* Loosely based upon NArray */
#ifndef NMATRIX_C
# define NMATRIX_C

#include <ruby.h>

#include "nmatrix.h"

VALUE cNMatrix;


const int nm_sizeof[NM_TYPES+1] = {
  0, /* NM_NONE */
  sizeof(u_int8_t), /* NM_BYTE */
  sizeof(int8_t),   /* NM_INT8 */
  sizeof(int16_t),  /* NM_INT16 */
  sizeof(int32_t),  /* NM_INT32 */
  sizeof(int64_t),  /* NM_INT64 */
  sizeof(float),    /* NM_FLOAT32 */
  sizeof(double),   /* NM_FLOAT64 */
  sizeof(complex64), /* NM_COMPLEX64 */
  sizeof(complex128),
  sizeof(rational32),
  sizeof(rational64),
  sizeof(rational128),
  sizeof(VALUE),     /* NM_ROBJ */
  0                  /* NM_TYPES */
};

const char *nm_dtypestring[] = {
  "none",
  "byte",	/* 1 */
  "int8",   /* 2 */
  "int16",	/* 3 */
  "int32",	/* 4 */
  "int64",  /* 5 */
  "float32",	/* 6 */
  "float64",	/* 7 */
  "complex64",	/* 8 */
  "complex128",	/* 9 */
  "rational32", /* 11 */
  "rational64", /* 12 */
  "rational128",  /* 13 */
  "object",	/* 14 */
  "ntypes"	/* 15 */
};


const char *nm_stypestring[] = {
  "dense",
  "list",
  "compressed",
  "stypes"
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


VALUE nm_dense_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, VALUE self) {
  NMATRIX* matrix = nm_create(dtype, S_DENSE, create_dense_storage(nm_sizeof[dtype], shape, rank, init_val));
  return Data_Wrap_Struct(self, NULL, nm_delete, matrix);
}

/* static VALUE nm_list_new(int argc, VALUE* argv, int8_t dtype, VALUE klass) {

} */


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
int8_t nm_interpret_dtype(int argc, VALUE* argv) {
  if (argc == 2) {
    if (SYMBOL_P(argv[1])) return nm_dtypesymbol_to_dtype(argv[1]);
    else if (IS_STRING(argv[1])) return nm_dtypestring_to_dtype(StringValue(argv[1]));
    else return nm_guess_dtype(argv[0]);
  } else if (argc == 1) {
    if (SYMBOL_P(argv[0])) return nm_dtypesymbol_to_dtype(argv[0]);
    else if (IS_STRING(argv[0])) return nm_dtypestring_to_dtype(StringValue(argv[0]));
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
  switch(dtype) {
  case NM_BYTE:
  case NM_INT8:
    *((int8_t*)(init_val)) = NUM2CHR(arg);
    return init_val;
  case NM_INT16:
    *((int16_t*)(init_val)) = NUM2INT(arg);
    return init_val;
  case NM_INT32:
    *((int32_t*)(init_val)) = NUM2INT(arg);
    return init_val;
  case NM_INT64:
    *((int64_t*)(init_val)) = NUM2INT(arg);
    return init_val;
  case NM_FLOAT32:
    *((float*)(init_val)) = NUM2FLT(arg);
    return init_val;
  case NM_FLOAT64:
    *((double*)(init_val)) = NUM2DBL(arg);
    return init_val;
  case NM_COMPLEX64:
  case NM_COMPLEX128:
  case NM_RATIONAL32:
  case NM_RATIONAL64:
  case NM_RATIONAL128:
  case NM_ROBJ:
  case NM_TYPES:
  case NM_NONE:
    rb_raise(rb_eNotImpError, "Type not implemented");
  }
  return init_val;
}


VALUE nm_new(int argc, VALUE* argv, VALUE self) {
  int8_t  dtype, stype, argp = 0;
  size_t  i, rank, elem_size = 1;
  size_t* shape;
  VALUE*  init_val_arg = NULL;
  void*   init_val = NULL;

  // READ ARGUMENTS
  if (!SYMBOL_P(argv[0]) && !IS_STRING(argv[0])) {
    stype    = S_DENSE;                                       // Dense by default.
    argp--;
  } else
    stype    = nm_interpret_stype(argv[0]);                   // 1: String or Symbol

  shape      = nm_interpret_shape_arg(argv[argp+1], &rank);   // 2: Either Fixnum or Array

  dtype      = nm_interpret_dtype(argc-2-argp, argv+2+argp);  // 3-4: dtype

  if (argc == 4 || (argc == 3 && IS_NUMERIC(argv[argp+2])))
    init_val = nm_interpret_initial_value(argv[argp+2], dtype);// 3: initial value / dtype
  else
    init_val = malloc(nm_sizeof[dtype]);

  if (argc > 4) {
    rb_raise(rb_eArgError, "Too many arguments");
    return Qnil;
  }

  // TODO: Update to allow an array as the initial value.

  if (dtype == NM_NONE) {
    rb_raise(rb_eArgError, "Could not recognize dtype");
    return Qnil;
  }


  if (stype == S_DENSE) // nm_dense_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, VALUE klass)
    return nm_dense_new(shape, rank, dtype, init_val, self);
  //else if (stype == S_LIST)
  //  nm_list_new(argc-1, argv+1, dtype, self);
  else
    rb_raise(rb_eNotImpError, "Only dense and list currently implemented");

  return Qnil;

}


// Does not create storage, but does destroy it.
NMATRIX* nm_create(int8_t dtype, int8_t stype, void* storage) {
  NMATRIX* mat = ALLOC(NMATRIX);

  mat->dtype   = dtype;
  mat->stype   = stype;
  mat->storage = storage;

  return mat;
}

static void nm_delete(NMATRIX* mat) {
  if (mat->stype == S_DENSE) delete_dense_storage(mat->storage);
  else if (mat->stype == S_LIST) delete_list_storage(mat->storage);
  else
    rb_raise(rb_eNotImpError, "Only dense and list deletion are implemented");
}



void Init_nmatrix() {
    /* Require Complex class */
    //rb_require("complex");
    //cComplex = rb_const_get( rb_cObject, rb_intern("Complex") );

    /* Define NMatrix class */
    cNMatrix = rb_define_class("NMatrix", rb_cObject);

    /* class methods */
    rb_define_singleton_method(cNMatrix, "new", nm_new, -1);

    /* methods */
    //rb_define_method(cNMatrix, "[]", nm_mref, -1);
    //rb_define_method(cNMatrix, "[]=", nm_mset, -1);
}

#endif