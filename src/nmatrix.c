/* Loosely based upon NArray */
#ifndef NMATRIX_C
# define NMATRIX_C

#include <ruby.h>

#include "nmatrix.h"

static void nm_delete(NMATRIX* mat);

VALUE cNMatrix;

ID nm_id_real, nm_id_imag;
ID nm_id_numer, nm_id_denom;


#include "dtypes.c"


const char *nm_stypestring[] = {
  "dense",
  "list",
  "compressed",
  "stypes"
};


VALUE nm_dense_get(STORAGE* s, size_t* coords, int8_t dtype) {
  VALUE v;
  SetFuncs[NM_ROBJ][dtype](1, &v, sizeof(VALUE), dense_storage_get((DENSE_STORAGE*)s, coords, nm_sizeof[dtype]), nm_sizeof[dtype]);
  return v;
}

VALUE nm_list_get(STORAGE* s, size_t* coords, int8_t dtype) {
  VALUE v;
  SetFuncs[NM_ROBJ][dtype](1, &v, sizeof(VALUE), list_storage_get((LIST_STORAGE*)s, coords), nm_sizeof[dtype]);
  return v;
}


nm_stype_ref_t RefFuncs = {
  nm_dense_get,
  nm_list_get,
  NULL,
  NULL
};


VALUE nm_dense_set(STORAGE* s, size_t* coords, VALUE val, int8_t dtype) {
  void* v = malloc(nm_sizeof[dtype]);

  SetFuncs[dtype][NM_ROBJ](1, v, nm_sizeof[dtype], &val, sizeof(VALUE));

  dense_storage_set( (DENSE_STORAGE*)s, coords, v, nm_sizeof[dtype] );
  free(v); // dense makes a copy, so free it

  return val;
}


VALUE nm_list_set(STORAGE* s, size_t* coords, VALUE val, int8_t dtype) {
  void* v = malloc(nm_sizeof[dtype]);

  SetFuncs[dtype][NM_ROBJ](1, v, nm_sizeof[dtype], &val, sizeof(VALUE));

  if (list_storage_insert( (LIST_STORAGE*)s, coords, v ))    return val;
  else                                                       return Qnil;
  // No need to free; the list keeps v.
}


nm_stype_ref_t InsFuncs = {
  nm_dense_set,
  nm_list_set,
  NULL,
  NULL
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

VALUE nm_list_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, VALUE self) {
  NMATRIX* matrix = nm_create(dtype, S_LIST, create_list_storage(nm_sizeof[dtype], shape, rank, init_val));
  return Data_Wrap_Struct(self, NULL, nm_delete, matrix);
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
  SetFuncs[dtype][NM_ROBJ](1, init_val, nm_sizeof[dtype], &arg, sizeof(VALUE));
  return init_val;
}


VALUE nm_new(int argc, VALUE* argv, VALUE self) {
  char    ZERO = 0;
  int8_t  dtype, stype, offset = 0;
  size_t  rank;
  size_t* shape;
  void*   init_val = NULL;

  // READ ARGUMENTS

  if (argc < 2 || argc > 4) { rb_raise(rb_eArgError, "Expected 2, 3, or 4 arguments"); return Qnil; }

  if (!SYMBOL_P(argv[0]) && !IS_STRING(argv[0])) {
    stype  = S_DENSE;
  } else {
    stype  = nm_interpret_stype(argv[0]);                        // 0: String or Symbol
    offset = 1;
  }
  shape    = nm_interpret_shape_arg(argv[offset], &rank);        // 1: String or Symbol
  dtype    = nm_interpret_dtype(argc-1-offset, argv+offset+1);   // 2-3: dtype

  if (IS_NUMERIC(argv[1+offset])) { // initial provided
    init_val = nm_interpret_initial_value(argv[1+offset], dtype);// 4: initial value / dtype
  } else { // dtype was provided, use 0
    init_val = malloc(nm_sizeof[dtype]);
    SetFuncs[dtype][NM_BYTE](1, init_val, nm_sizeof[dtype], &ZERO, sizeof(char));
  }


  // TODO: Update to allow an array as the initial value.

  if (dtype == NM_NONE) {
    rb_raise(rb_eArgError, "Could not recognize dtype");
    return Qnil;
  }


  if (stype == S_DENSE) // nm_dense_new(size_t* shape, size_t rank, int8_t dtype, void* init_val, VALUE klass)
    return nm_dense_new(shape, rank, dtype, init_val, self);
  else if (stype == S_LIST)
    return nm_list_new(shape, rank, dtype, init_val, self);
  else
    rb_raise(rb_eNotImpError, "Only dense and list currently implemented");

  return Qnil;

}


static void nm_delete(NMATRIX* mat) {
  if (mat->stype == S_DENSE)     delete_dense_storage(mat->storage);
  else if (mat->stype == S_LIST) delete_list_storage(mat->storage);
  else                           rb_raise(rb_eNotImpError, "Only dense and list deletion are implemented");
}


// Does not create storage, but does destroy it.
NMATRIX* nm_create(int8_t dtype, int8_t stype, void* storage) {
  NMATRIX* mat = ALLOC(NMATRIX);

  mat->dtype   = dtype;
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
                                         convert_coords((size_t)(argc), argv, self),
                                         NM_DTYPE(self) );

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
                                         argv[rank], // value to set it to
                                         NM_DTYPE(self) );

  } else if (NM_RANK(self) < rank) {
    rb_raise(rb_eArgError, "Coordinates given exceed matrix rank");
  } else {
    rb_raise(rb_eNotImpError, "Slicing not supported yet");
  }
  return Qnil;
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
    rb_define_method(cNMatrix, "[]", nm_mref, -1);
    rb_define_method(cNMatrix, "[]=", nm_mset, -1);

    nm_id_real  = rb_intern("real");
    nm_id_imag  = rb_intern("imag");
    nm_id_numer = rb_intern("numerator");
    nm_id_denom = rb_intern("denominator");

}

#endif