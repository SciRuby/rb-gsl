/*
  gsl_nmatrix.c

  Written by John Woods     Jul/2012
 */

#include "rb_gsl_config.h"
#ifdef HAVE_NMATRIX_H
#include "rb_gsl_array.h"
#include "nmatrix.h"
#include "rb_gsl_with_nmatrix.h"

// TODO: Figure out how to import these two enums from nmatrix
typedef enum {
  NM_DENSE_STORE,
  NM_LIST_STORE,
  NM_YALE_STORE
} nm_stype_t;

typedef enum {
	NM_BYTE				  =  0, // unsigned char
	NM_INT8				  =  1, // char
	NM_INT16				=  2, // short
	NM_INT32				=  3, // int
	NM_INT64				=  4, // long
	NM_FLOAT32			=  5, // float
	NM_FLOAT64			=  6, // double
	NM_COMPLEX64		=  7, // Complex64 class
	NM_COMPLEX128	  =  8, // Complex128 class
	NM_RATIONAL32	  =  9, // Rational32 class
	NM_RATIONAL64	  = 10, // Rational64 class
	NM_RATIONAL128	= 11, // Rational128 class
	NM_RUBYOBJ			= 12  // Ruby VALUE type
} nm_dtype_t;

/* GSL::Vector -> NMatrix */

static VALUE rb_gsl_vector_to_nvector(VALUE obj, VALUE klass) {
  gsl_vector *v = NULL;
  Data_Get_Struct(obj, gsl_vector, v);

  return rb_nvector_dense_create(NM_FLOAT64, v->data, v->size);
}


static VALUE rb_gsl_vector_complex_to_nvector(VALUE obj, VALUE klass) {
  gsl_vector *v = NULL;
  Data_Get_Struct(obj, gsl_vector, v);

  return rb_nvector_dense_create(NM_COMPLEX128, v->data, v->size);
}


static VALUE rb_gsl_vector_int_to_nvector(VALUE obj, VALUE klass) {
  gsl_vector *v = NULL;
  Data_Get_Struct(obj, gsl_vector, v);

  return rb_nvector_dense_create(NM_INT64, v->data, v->size);
}


static VALUE rb_gsl_matrix_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_matrix *m = NULL;
  Data_Get_Struct(obj, gsl_matrix, m);

  return rb_nmatrix_dense_create(NM_FLOAT64, &(m->size1), 2, m->data, m->size1 * m->size2);
}


static VALUE rb_gsl_matrix_complex_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_matrix *m = NULL;
  Data_Get_Struct(obj, gsl_matrix, m);

  return rb_nmatrix_dense_create(NM_COMPLEX128, &(m->size1), 2, m->data, m->size1 * m->size2);
}


static VALUE rb_gsl_matrix_int_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_matrix *m = NULL;
  Data_Get_Struct(obj, gsl_matrix, m);

  return rb_nmatrix_dense_create(NM_INT64, &(m->size1), 2, m->data, m->size1 * m->size2);
}


gsl_vector* nv_to_gv(VALUE nm) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)(nm->storage);
  gsl_vector* v = gsl_vector_alloc( s->count );

  if (s->dtype != FLOAT64) {
    rb_raise(nm_eDataTypeError, "requires dtype of :float64 to convert to a GSL vector");
  }

  memcpy(v->data, s->elements, v->size*sizeof(double));

  return v;
}


gsl_vector_complex* nv_to_gv_complex(VALUE nm) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)(nm->storage);
  gsl_vector_complex* v = gsl_vector_complex_alloc( s->count );

  if (s->dtype != COMPLEX128) {
    rb_raise(nm_eDataTypeError, "requires dtype of :complex128 to convert to a GSL complex vector");
  }

  memcpy(v->data, s->elements, v->size*sizeof(double)*2);

  return v;
}


gsl_vector_int* nv_to_gv_int(VALUE nm) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)(nm->storage);
  gsl_vector_int* v = gsl_vector_int_alloc( s->count );

  if (s->dtype != INT64) {
    rb_raise(nm_eDataTypeError, "requires dtype of :int64 to convert to a GSL int vector");
  }

  memcpy(v->data, s->elements, v->size*sizeof(int));

  return v;
}


gsl_matrix* nm_to_gm(VALUE nm) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)(nm->storage);
  gsl_matrix* m = gsl_matrix_alloc( s->shape[0], s->shape[1] );

  if (s->dtype != FLOAT64) {
    rb_raise(nm_eDataTypeError, "requires dtype of :float64 to convert from a GSL double vector");
  }

  memcpy(m->data, s->elements, s->count);
  return m;
}

gsl_matrix_complex* nm_to_gm_complex(VALUE nm) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)(nm->storage);
  gsl_matrix_complex* m = gsl_matrix_complex_alloc( s->shape[0], s->shape[1] );

  if (s->dtype != COMPLEX128) {
    rb_raise(nm_eDataTypeError, "requires dtype of :complex128 to convert from a GSL complex vector");
  }

  memcpy(m->data, s->elements, s->count);
  return m;
}


gsl_matrix_int* nm_to_gm_int(VALUE nm) {
  DENSE_STORAGE* s = (DENSE_STORAGE*)(nm->storage);
  gsl_matrix_int* m = gsl_matrix_int_alloc( s->shape[0], s->shape[1] );

  if (s->dtype != INT64) {
    rb_raise(nm_eDataTypeError, "requires dtype of :int64 to convert from a GSL int vector");
  }

  memcpy(m->data, s->elements, s->count);
  return m;
}


static VALUE rb_gsl_nv_to_gsl_vector(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, nv_to_gv(n));
}

static VALUE rb_gsl_nv_to_gsl_vector_complex(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, nv_to_gv_complex(n));
}

static VALUE rb_gsl_nv_to_gsl_vector_int(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, nm_to_gv_int(n));
}


static VALUE rb_gsl_nm_to_gsl_matrix(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, nm_to_gm(n));
}

static VALUE rb_gsl_nm_to_gsl_matrix_complex(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, nm_to_gm_complex(n));
}

static VALUE rb_gsl_nm_to_gsl_matrix_int(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_matrix_int, 0, gsl_matrix_int_free, nm_to_gm_int(n));
}





void Init_gsl_nmatrix(VALUE module) {
  rb_define_method(cgsl_vector, "to_nv", rb_gsl_vector_to_nvector, 0);
  rb_define_alias(cgsl_vector, "to_nm", "to_nv");

  rb_define_singleton_method(cgsl_vector, "to_gslv",    rb_gsl_nv_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector, "to_gv",      rb_gsl_nv_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector, "nv_to_gslv", rb_gsl_nv_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector, "nv_to_gv",   rb_gsl_nv_to_gsl_vector, 1);

  rb_define_method(cgsl_vector_complex, "to_nv", rb_gsl_vector_to_nvector_complex, 0);
  rb_define_alias(cgsl_vector_complex, "to_nm", "to_nv");

  rb_define_singleton_method(cgsl_vector_complex, "to_gslv",    rb_gsl_nv_to_gsl_vector_complex, 1);
  rb_define_singleton_method(cgsl_vector_complex, "to_gv",      rb_gsl_nv_to_gsl_vector_complex, 1);
  rb_define_singleton_method(cgsl_vector_complex, "nv_to_gslv", rb_gsl_nv_to_gsl_vector_complex, 1);
  rb_define_singleton_method(cgsl_vector_complex, "nv_to_gv",   rb_gsl_nv_to_gsl_vector_complex, 1);

  rb_define_method(cgsl_vector_int, "to_nv", rb_gsl_vector_to_nvector_int, 0);
  rb_define_alias(cgsl_vector_int, "to_nm", "to_nv");

  rb_define_singleton_method(cgsl_vector_int, "to_gslv",    rb_gsl_nv_to_gsl_vector_int, 1);
  rb_define_singleton_method(cgsl_vector_int, "to_gv",      rb_gsl_nv_to_gsl_vector_int, 1);
  rb_define_singleton_method(cgsl_vector_int, "nv_to_gslv", rb_gsl_nv_to_gsl_vector_int, 1);
  rb_define_singleton_method(cgsl_vector_int, "nv_to_gv",   rb_gsl_nv_to_gsl_vector_int, 1);

  rb_define_method(cgsl_matrix, "to_nm", rb_gsl_matrix_to_nmatrix, 0);
  rb_define_singleton_method(cgsl_matrix, "nm_to_gslm",  rb_gsl_nm_to_gsl_matrix, 1);

  rb_define_method(cgsl_matrix_complex, "to_nm", rb_gsl_matrix_complex_to_nmatrix, 0);
  rb_define_singleton_method(cgsl_matrix_complex, "nm_to_gslm",  rb_gsl_nm_to_gsl_matrix_complex, 1);

  rb_define_method(cgsl_matrix_int, "to_nm", rb_gsl_matrix_int_to_nmatrix, 0);
  rb_define_singleton_method(cgsl_matrix_int, "nm_to_gslm",  rb_gsl_nm_to_gsl_matrix_int, 1);
}

#endif // HAVE_NMATRIX_H
