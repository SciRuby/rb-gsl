/*
  gsl_nmatrix.c

  Written by Sameer Deshmukh (@v0dro)      feb/2016
 */
#ifdef HAVE_NMATRIX_H

#include "include/rb_gsl_with_nmatrix.h"

// nmatrix external API
extern VALUE rb_nmatrix_dense_create(nm_dtype_t dtype, size_t* shape, size_t rank, void* elements, size_t length);
extern VALUE rb_nvector_dense_create(nm_dtype_t dtype, void* elements, size_t length);

/* GSL::Vector -> NMatrix */
static VALUE rb_gsl_vector_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_vector *v = NULL;
  Data_Get_Struct(obj, gsl_vector, v);

  return rb_nvector_dense_create(FLOAT64, v->data, v->size);
}

static VALUE rb_gsl_vector_int_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_vector_int *v = NULL;
  Data_Get_Struct(obj, gsl_vector_int, v);

  return rb_nvector_dense_create(INT32, v->data, v->size);
}

static VALUE rb_gsl_vector_complex_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_vector_complex *v = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);

  return rb_nvector_dense_create(COMPLEX128, v->data, v->size);
}

static VALUE rb_gsl_matrix_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_matrix *m = NULL;
  Data_Get_Struct(obj, gsl_matrix, m);

  return rb_nmatrix_dense_create(FLOAT64, &(m->size1), 2, m->data, m->size1 * m->size2);
}

static VALUE rb_gsl_matrix_int_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_matrix_int *m = NULL;
  Data_Get_Struct(obj, gsl_matrix_int, m);

  return rb_nmatrix_dense_create(INT64, &(m->size1), 2, m->data, m->size1 * m->size2);
}

static VALUE rb_gsl_matrix_complex_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_matrix_complex *m = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);

  return rb_nmatrix_dense_create(COMPLEX128, &(m->size1), 2, m->data, m->size1 * m->size2);
}

gsl_vector* nm_to_gv(VALUE nm) {
  NM_DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);
  // TODO: change this once nmatrix is fixed
  gsl_vector* v = gsl_vector_alloc( FIX2INT(rb_funcall(nm, rb_intern("count"), 0, NULL)) );

  if (s->dtype != FLOAT64) {
    rb_raise(nm_eDataTypeError, "requires dtype of :float64 to convert to a GSL vector");
  }

  memcpy(v->data, s->elements, v->size*sizeof(double));

  return v;
}

gsl_vector_int* nm_to_gv_int(VALUE nm) {
  NM_DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);
  gsl_vector_int* v = gsl_vector_int_alloc( FIX2INT(rb_funcall(nm, rb_intern("count"), 0, NULL)) );

  if (s->dtype != INT32) {
    rb_raise(nm_eDataTypeError, "requires dtype of :int32 to convert to a GSL int vector");
  }

  memcpy(v->data, s->elements, v->size*sizeof(int32_t));

  return v;
}

gsl_vector_complex* nm_to_gv_complex(VALUE nm) {
  NM_DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);
  gsl_vector_complex* v = gsl_vector_complex_alloc( FIX2INT(rb_funcall(nm, rb_intern("count"), 0, NULL)) );

  if (s->dtype != COMPLEX128) {
    rb_raise(nm_eDataTypeError, "requires dtype of :complex128 to convert to a GSL complex vector");
  }

  memcpy(v->data, s->elements, v->size*sizeof(double)*2);

  return v;
}

gsl_matrix* nm_to_gm(VALUE nm) {
  NM_DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);
  gsl_matrix* m = gsl_matrix_alloc( s->shape[0], s->shape[1] );

  if (s->dtype != FLOAT64) {
    rb_raise(nm_eDataTypeError, "requires dtype of :float64 to convert from a GSL double vector");
  }

  memcpy(m->data, s->elements, s->count);
  return m;
}

gsl_matrix_int* nm_to_gm_int(VALUE nm) {
  NM_DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);
  gsl_matrix_int* m = gsl_matrix_int_alloc( s->shape[0], s->shape[1] );

  if (s->dtype != INT32) {
    rb_raise(nm_eDataTypeError, "requires dtype of :int32 to convert from a GSL int vector");
  }

  memcpy(m->data, s->elements, s->count);
  return m;
}

gsl_matrix_complex* nm_to_gm_complex(VALUE nm) {
  NM_DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);
  gsl_matrix_complex* m = gsl_matrix_complex_alloc( s->shape[0], s->shape[1] );

  if (s->dtype != COMPLEX128) {
    rb_raise(nm_eDataTypeError, "requires dtype of :complex128 to convert from a GSL complex vector");
  }

  memcpy(m->data, s->elements, s->count);
  return m;
}

VALUE rb_gsl_nm_to_gsl_vector(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, nm_to_gv(n));
}

VALUE rb_gsl_nm_to_gsl_vector_int(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, nm_to_gv_int(n));
}

VALUE rb_gsl_nm_to_gsl_vector_complex(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, nm_to_gv_complex(n));
}

VALUE rb_gsl_nm_to_gsl_matrix(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, nm_to_gm(n));
}

VALUE rb_gsl_nm_to_gsl_matrix_int(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_matrix_int, 0, gsl_matrix_int_free, nm_to_gm_int(n));
}

VALUE rb_gsl_nm_to_gsl_matrix_complex(VALUE obj, VALUE n) {
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, nm_to_gm_complex(n));
}

VALUE rb_gsl_nm_to_gsl_vector_method(VALUE nm) {
  VALUE v;

  if (NM_DTYPE(nm) == COMPLEX64 || NM_DTYPE(nm) == COMPLEX128) {
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, nm_to_gv_complex(nm));
  } 
  else if (NM_DTYPE(nm) == INT32) {
    return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, nm_to_gv_int(nm));
  }
  else {
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, nm_to_gv(nm));
  }

  return v;
}

Init_gsl_nmatrix(VALUE module)
{
  rb_define_method(cgsl_vector, "to_nm", rb_gsl_vector_to_nmatrix, 0);
  rb_define_method(cgsl_vector_int, "to_nm", rb_gsl_vector_int_to_nmatrix, 0);
  rb_define_method(cgsl_vector_complex, "to_nm", rb_gsl_vector_complex_to_nmatrix, 0);

  rb_define_singleton_method(cgsl_vector, "nm_to_gslv", rb_gsl_nm_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector_int, "nm_to_gslv", rb_gsl_nm_to_gsl_vector_int, 1);
  rb_define_singleton_method(cgsl_vector_complex, "nm_to_gslv", rb_gsl_nm_to_gsl_vector_complex, 1);
  
  rb_define_method(cgsl_matrix, "to_nm", rb_gsl_matrix_to_nmatrix, 0);
  rb_define_method(cgsl_matrix_int, "to_nm", rb_gsl_matrix_int_to_nmatrix, 0);
  rb_define_method(cgsl_matrix_complex, "to_nm", rb_gsl_matrix_complex_to_nmatrix, 0);

  rb_define_singleton_method(cgsl_matrix, "nm_to_gslm",  rb_gsl_nm_to_gsl_matrix, 1);
  rb_define_singleton_method(cgsl_matrix_int, "nm_to_gslm",  rb_gsl_nm_to_gsl_matrix_int, 1);
  rb_define_singleton_method(cgsl_matrix_complex, "nm_to_gslm",  rb_gsl_nm_to_gsl_matrix_complex, 1);

  rb_define_method(cNMatrix, "to_gslv", rb_gsl_nm_to_gsl_vector_method, 0);
}

#endif