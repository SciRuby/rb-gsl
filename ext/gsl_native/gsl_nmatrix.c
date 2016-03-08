/*
  gsl_nmatrix.c

  Written by Sameer Deshmukh (@v0dro)      feb/2016
 */
#ifdef HAVE_NMATRIX_H

#include "include/rb_gsl_with_nmatrix.h"

// functions to convert GSL::Vectors to 1D NMatrix
static VALUE rb_gsl_vector_to_nmatrix(VALUE obj);
static VALUE rb_gsl_vector_int_to_nmatrix(VALUE obj);
static VALUE rb_gsl_vector_complex_to_nmatrix(VALUE obj);

// functions to convert GSL::Matrix to 2D NMatrix
static VALUE rb_gsl_matrix_to_nmatrix(VALUE obj);
static VALUE rb_gsl_matrix_int_to_nmatrix(VALUE obj);
static VALUE rb_gsl_matrix_complex_to_nmatrix(VALUE obj);

/* GSL::Vector -> NMatrix */
static VALUE rb_gsl_vector_to_nmatrix(VALUE obj) {
  gsl_vector *v = NULL;
  Data_Get_Struct(obj, gsl_vector, v);

  return rb_nvector_dense_create(FLOAT64, v->data, v->size);
}

static VALUE rb_gsl_vector_int_to_nmatrix(VALUE obj) {
  gsl_vector_int *v = NULL;
  Data_Get_Struct(obj, gsl_vector_int, v);

  return rb_nvector_dense_create(INT32, v->data, v->size);
}

static VALUE rb_gsl_vector_complex_to_nmatrix(VALUE obj) {
  gsl_vector_complex *v = NULL;
  Data_Get_Struct(obj, gsl_vector_complex, v);

  return rb_nvector_dense_create(COMPLEX128, v->data, v->size);
}

static VALUE rb_gsl_matrix_to_nmatrix(VALUE obj) {
  gsl_matrix *m = NULL;
  Data_Get_Struct(obj, gsl_matrix, m);

  return rb_nmatrix_dense_create(FLOAT64, &(m->size1), 2, m->data, m->size1 * m->size2);
}

static VALUE rb_gsl_matrix_int_to_nmatrix(VALUE obj) {
  gsl_matrix_int *m = NULL;
  Data_Get_Struct(obj, gsl_matrix_int, m);

  return rb_nmatrix_dense_create(INT32, &(m->size1), 2, m->data, m->size1 * m->size2);
}

static VALUE rb_gsl_matrix_complex_to_nmatrix(VALUE obj) {
  gsl_matrix_complex *m = NULL;
  Data_Get_Struct(obj, gsl_matrix_complex, m);

  return rb_nmatrix_dense_create(COMPLEX128, &(m->size1), 2, m->data, m->size1 * m->size2);
}

// functions called on NMatrix objects. 'nm' is of type NMatrix.
gsl_vector* rb_gsl_nmatrix_to_gv(VALUE nm) {
  NM_DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);
  // TODO: change this once nmatrix is fixed
  gsl_vector* v = gsl_vector_alloc( FIX2INT(rb_funcall(nm, rb_intern("count"), 0, NULL)) );

  if (s->dtype != FLOAT64) {
    rb_raise(rb_eStandardError, "requires dtype of :float64 to convert to a GSL vector");
  }

  memcpy(v->data, s->elements, v->size*sizeof(double));

  return v;
}

gsl_vector_int* rb_gsl_nmatrix_to_gv_int(VALUE nm) {
  NM_DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);
  gsl_vector_int* v = gsl_vector_int_alloc(
    FIX2INT(rb_funcall(nm, rb_intern("count"), 0, NULL)));

  if (s->dtype != INT32) {
    rb_raise(rb_eStandardError, "requires dtype of :int32 to convert to a GSL int vector");
  }

  memcpy(v->data, s->elements, v->size*sizeof(int32_t));

  return v;
}

gsl_vector_complex* rb_gsl_nmatrix_to_gv_complex(VALUE nm) {
  NM_DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);
  gsl_vector_complex* v = gsl_vector_complex_alloc(
    FIX2INT(rb_funcall(nm, rb_intern("count"), 0, NULL)) );

  if (s->dtype != COMPLEX128) {
    rb_raise(rb_eStandardError, "requires dtype of :complex128 to convert to a GSL complex vector");
  }

  memcpy(v->data, s->elements, v->size*sizeof(double)*2);

  return v;
}

// NMatrix function to convert NMatrix to GSL::Vector object depending on dtype
VALUE rb_gsl_nmatrix_to_gsl_vector_method(VALUE nm) {
  if (NM_DTYPE(nm) == COMPLEX64 || NM_DTYPE(nm) == COMPLEX128) {
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, 
      rb_gsl_nmatrix_to_gv_complex(nm));
  } 
  else if (NM_DTYPE(nm) == INT32) {
    return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free,
      rb_gsl_nmatrix_to_gv_int(nm));
  }
  else if (NM_DTYPE(nm) == FLOAT64) {
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free,
      rb_gsl_nmatrix_to_gv(nm));
  }
  else {
    rb_raise(rb_eStandardError, 
      "NMatrix should be either :complex64, :complex128, :int32 or :float64 type."); 
  }
}

gsl_matrix* rb_gsl_nmatrix_to_gm(VALUE nm) {
  NM_DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);
  gsl_matrix* m = gsl_matrix_alloc( s->shape[0], s->shape[1] );

  if (s->dtype != FLOAT64) {
    rb_raise(rb_eStandardError, "requires dtype of :float64 to convert from a GSL double vector");
  }

  memcpy(m->data, s->elements, s->shape[0]*s->shape[1]*sizeof(double)); // double is nm :float64
  return m;
}

gsl_matrix_int* rb_gsl_nmatrix_to_gm_int(VALUE nm) {
  NM_DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);
  gsl_matrix_int* m = gsl_matrix_int_alloc( s->shape[0], s->shape[1] );

  if (s->dtype != INT32) {
    rb_raise(rb_eStandardError, "requires dtype of :int32 to convert from a GSL int vector");
  }

  memcpy(m->data, s->elements, s->shape[0]*s->shape[1]*sizeof(int32_t)); // int32_t is nm :int32
  return m;
}

gsl_matrix_complex* rb_gsl_nmatrix_to_gm_complex(VALUE nm) {
  NM_DENSE_STORAGE* s = NM_STORAGE_DENSE(nm);
  gsl_matrix_complex* m = gsl_matrix_complex_alloc( s->shape[0], s->shape[1]);

  if (s->dtype != COMPLEX128) {
    rb_raise(rb_eStandardError, "requires dtype of :complex128 to convert from a GSL complex vector");
  }

  memcpy(m->data, s->elements, s->shape[0]*s->shape[1]*sizeof(double)*2);
  return m;
}

// NMatrix to GSL::Matrix conversion based on dtype
VALUE rb_gsl_nmatrix_to_gsl_matrix_method(VALUE nmatrix) {
  if (NM_DIM(nmatrix) > 2) {
    rb_raise(rb_eStandardError,
      "NMatrix must not have greater than 2 dimensions.");
  }

  if (NM_DTYPE(nmatrix) == COMPLEX64 || NM_DTYPE(nmatrix) == COMPLEX128) {
    return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, 
      rb_gsl_nmatrix_to_gm_complex(nmatrix));
  } 
  else if (NM_DTYPE(nmatrix) == INT32) {
    return Data_Wrap_Struct(cgsl_matrix_int, 0, gsl_matrix_int_free,
      rb_gsl_nmatrix_to_gm_int(nmatrix));
  }
  else if (NM_DTYPE(nmatrix) == FLOAT64) {
    return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free,
      rb_gsl_nmatrix_to_gm(nmatrix));
  }
  else {
    rb_raise(rb_eStandardError, 
      "NMatrix should be either :complex64, :complex128, :int32 or :float64 type."); 
  }
}

void Init_gsl_nmatrix(VALUE module)
{
  // functions to convert GSL::Vector to 1D NMatrix.
  rb_define_method(cgsl_vector, "to_nm", rb_gsl_vector_to_nmatrix, 0);
  rb_define_method(cgsl_vector_int, "to_nm", rb_gsl_vector_int_to_nmatrix, 0);
  rb_define_method(cgsl_vector_complex, "to_nm", rb_gsl_vector_complex_to_nmatrix, 0);

  // functions to convert GSL::Matrix to 2D Nmatrix.
  rb_define_method(cgsl_matrix, "to_nm", rb_gsl_matrix_to_nmatrix, 0);
  rb_define_method(cgsl_matrix_int, "to_nm", rb_gsl_matrix_int_to_nmatrix, 0);
  rb_define_method(cgsl_matrix_complex, "to_nm", rb_gsl_matrix_complex_to_nmatrix, 0);

  // functions on NMatrix to convert to GSL::Vector and GSL::Matrix
  rb_define_method(cNMatrix, "to_gslv", rb_gsl_nmatrix_to_gsl_vector_method, 0);
  rb_define_method(cNMatrix, "to_gslm", rb_gsl_nmatrix_to_gsl_matrix_method, 0);
}

#endif