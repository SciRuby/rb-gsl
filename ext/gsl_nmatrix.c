/*
  gsl_nmatrix.c

  Written by John Woods     Jul/2012
 */

#include "rb_gsl_config.h"
#ifdef HAVE_NMATRIX_H
#include "rb_gsl_array.h"
#include "nmatrix.h"
#include "rb_gsl_with_nmatrix.h"


/* GSL::Vector -> NMatrix */

static VALUE rb_gsl_vector_to_nmatrix(VALUE obj, VALUE klass) {
  gsl_vector *v = NULL;
  NMATRIX* nm;
  size_t* shape;
  size_t elements_length;
  double* elements;
  unsigned int i;
  Data_Get_Struct(obj, gsl_vector, v);
  shape = ALLOC_N(size_t, 1);
  shape[0] = v->size;

  if (v->stride == 1) {
    memcpy(elements, v->data, shape[0] * sizeof(double));
  } else {
    for (i = 0; i < v->size; ++i) {
      elements[i] = gsl_vector_get(v, i);
    }
  }

  nm = nm_create(DENSE_STORE, dense_storage_create(FLOAT64, shape, 1, elements, v->size));

  return Data_Wrap_Struct(cNMatrix, mark_dense_storage, nm_delete, nm);
}


static VALUE rb_gsl_vector_complex_to_nmatrix(VALUE obj, VALUE klass) {
  // Mostly need to copy from above once it works.
  rb_raise(rb_eNotImpError, "needs complex GSL to nmatrix conversion code");
  return Qnil;
}


static VALUE rb_gsl_vector_to_nm(VALUE obj) {
  if (VECTOR_P(obj))
    return rb_gsl_vector_to_nmatrix(obj, cNMatrix);
  else if (VECTOR_COMPLEX_P(obj))
    return rb_gsl_vector_complex_to_nmatrix(obj, CNMatrix);
  else
    rb_raise(rb_eRuntimeError, "unexpected type '%s'", rb_obj_classname(obj));

  return Qnil;
}


static VALUE rb_gsl_vector_complex_to_nm(VALUE obj) {
  return rb_gsl_vector_complex_to_nmatrix(obj, cNMatrix);
}


static VALUE rb_gsl_vector_to_nmatrix(VALUE obj) {
  return rb_gsl_vector_to_nmatrix(obj, cNMatrix);
}

static VALUE rb_gsl_vector_complex_to_nmatrix(VALUE obj) {
  return rb_gsl_vector_to_nmatrix(obj, cNMatrix);
}


gsl_vector* nm_to_gv(VALUE nm) {
  gsl_vector* v = NULL;
  VALUE nmat = nm;
  DENSE_STORAGE* s = (DENSE_STORAGE*)(nm->storage);
  v = gsl_vector_alloc( s->count );
  if (nm->s->dtype != FLOAT64) {
    rb_raise(rb_eNotImpError, "to do: change dtype code here");
  }
  memcpy(v->data, s->elements, v->size*sizeof(double));
  return v;
}


void Init_gsl_nmatrix(VALUE module) {
  rb_define_method(cgsl_vector, "to_nm", rb_gsl_vector_to_nm, 0);
  rb_define_alias(cgsl_vector, "to_nmatrix", "to_nm");

  rb_define_singleton_method(cgsl_vector, "to_gslv",    rb_gsl_nm_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector, "to_gv",      rb_gsl_nm_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector, "nm_to_gslv", rb_gsl_nm_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector, "nm_to_gv",   rb_gsl_nm_to_gsl_vector, 1);

  rb_define_method(csgl_vector_complex, "to_nm", rb_gsl_vector_complex_to_nm, 0);
  rb_define_alias(csgl_vector_complex, "to_nmatrix", "to_nm");

  /*rb_define_singleton_method(cgsl_vector_complex, "to_gslv",    rb_gsl_nm_to_gsl_vector_complex, 1);
  rb_define_singleton_method(cgsl_vector_complex, "to_gv",      rb_gsl_nm_to_gsl_vector_complex, 1);
  rb_define_singleton_method(cgsl_vector_complex, "nm_to_gslv", rb_gsl_nm_to_gsl_vector_complex, 1);
  rb_define_singleton_method(cgsl_vector_complex, "nm_to_gv",   rb_gsl_nm_to_gsl_vector_complex, 1); */
}

#endif