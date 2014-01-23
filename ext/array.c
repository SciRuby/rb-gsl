/*
  array.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_config.h"

#include "rb_gsl_common.h"
#include "rb_gsl_array.h"
#include "rb_gsl_complex.h"
#ifdef HAVE_NARRAY_H
#include "rb_gsl_with_narray.h"
#endif

/* global variables */
VALUE cgsl_block, cgsl_block_int;
VALUE cgsl_block_uchar;
VALUE cgsl_block_complex;
VALUE cgsl_vector, cgsl_vector_view;
VALUE cgsl_vector_view_ro;
VALUE cgsl_vector_col_view_ro;
VALUE cgsl_matrix_view_ro;
VALUE cgsl_vector_complex;
VALUE cgsl_vector_complex_view;
VALUE cgsl_matrix;
VALUE cgsl_matrix_view;
VALUE cgsl_matrix_complex;
VALUE cgsl_matrix_complex_view;
VALUE cgsl_vector_complex_view_ro;
VALUE cgsl_matrix_complex_view_ro;
VALUE cgsl_permutation;
VALUE cgsl_index;
VALUE cgsl_vector_col, cgsl_vector_col_view;
VALUE cgsl_vector_complex_col, cgsl_vector_complex_col_view;

VALUE cgsl_vector_int, cgsl_vector_int_col;
VALUE cgsl_vector_int_view, cgsl_vector_int_col_view;
VALUE cgsl_vector_int_view_ro, cgsl_vector_int_col_view_ro;
VALUE cgsl_matrix_int, cgsl_matrix_int_view, cgsl_matrix_int_view_ro;

double* get_vector_ptr(VALUE ary, size_t *stride, size_t *n)
{
  gsl_vector *v = NULL;
  gsl_vector_complex *vc = NULL;
  gsl_matrix *m;
#ifdef HAVE_NARRAY_H
  VALUE ary2;
#endif
  if (VECTOR_P(ary)) {
    Data_Get_Struct(ary, gsl_vector, v);
    *stride = v->stride;
    *n = v->size;
    return v->data;
  } else if (VECTOR_COMPLEX_P(ary)) {
    Data_Get_Struct(ary, gsl_vector_complex, vc);
    *stride = vc->stride;
    *n = vc->size*2;
    return vc->data;
  } else if (MATRIX_P(ary)) {
    Data_Get_Struct(ary, gsl_matrix, m);
    *stride = 1;
    *n = m->size1*m->size2;
    return m->data;
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(ary)) {
    *n = NA_TOTAL(ary);
    *stride = 1;
    ary2 = na_change_type(ary, NA_DFLOAT);
    return NA_PTR_TYPE(ary2,double*);
#endif
  } else {
    rb_raise(rb_eTypeError,
	     "wrong argument type %s", rb_class2name(CLASS_OF(ary)));
  }
}

gsl_vector* get_cvector(VALUE obj)
{
  gsl_vector *v = NULL;
  if (rb_obj_is_kind_of(obj, cgsl_vector)) {
    Data_Get_Struct(obj, gsl_vector,  v);
#ifdef HAVE_NARRAY_H
  } else if (NA_IsArray(obj)) {
    v = make_cvector_from_rarrays(obj);
#endif
  } else {
    rb_raise(rb_eTypeError,
             "wrong argument type %s", rb_class2name(CLASS_OF(obj)));
  }
  return v; 
}

VALUE make_rarray_from_cvector(const gsl_vector *v)
{
  size_t i;
  VALUE ary;
  ary = rb_ary_new2(v->size);
  for (i = 0; i < v->size; i++) {
    rb_ary_store(ary, i, rb_float_new(gsl_vector_get(v, i)));
  }
  return ary;
}

VALUE make_rarray_from_cvector_int(const gsl_vector_int *v)
{
  size_t i;
  VALUE ary;
  ary = rb_ary_new2(v->size);
  for (i = 0; i < v->size; i++) {
    rb_ary_store(ary, i, INT2FIX(gsl_vector_int_get(v, i)));
  }
  return ary;
}

VALUE make_rarray_from_cpermutation(const gsl_permutation *v)
{
  size_t i;
  VALUE ary;
  ary = rb_ary_new2(v->size);
  for (i = 0; i < v->size; i++) {
    rb_ary_store(ary, i, rb_float_new(gsl_permutation_get(v, i)));
  }
  return ary;
}

gsl_vector* get_vector(VALUE ary)
{
  gsl_vector *v = NULL;
  if (CLASS_OF(ary) == rb_cRange) ary = rb_gsl_range2ary(ary);
  if (TYPE(ary) == T_ARRAY) {
    return make_cvector_from_rarray(ary);
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(ary)) {
    return make_cvector_from_narray(ary);
#endif
  } else if (VECTOR_P(ary)) {
    Data_Get_Struct(ary, gsl_vector, v);
    return v;
  } else {
    rb_raise(rb_eTypeError,
	     "wrong argument type %s", rb_class2name(CLASS_OF(ary)));
  }
}

gsl_vector* make_cvector_from_rarrays(VALUE ary)
{
  if (CLASS_OF(ary) == rb_cRange) ary = rb_gsl_range2ary(ary);
  if (TYPE(ary) == T_ARRAY) {
    return make_cvector_from_rarray(ary);
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(ary)) {
    return make_cvector_from_narray(ary);
#endif
  } else {
    rb_raise(rb_eTypeError,
	     "wrong argument type %s", rb_class2name(CLASS_OF(ary)));
  }
}

void cvector_set_from_carray(gsl_vector *v, const double *a)
{
  size_t i;
  for (i = 0; i < v->size; i++) gsl_vector_set(v, i, a[i]);
}

void carray_set_from_cvector(double *a, const gsl_vector *v)
{
  size_t i;
  for (i = 0; i < v->size; i++) a[i] = gsl_vector_get(v, i);
}

void carray_set_from_rarrays(double *a, VALUE ary)
{
  if (CLASS_OF(ary) == rb_cRange) ary = rb_gsl_range2ary(ary);
  if (TYPE(ary) == T_ARRAY) {
    return carray_set_from_rarray(a, ary);
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(ary)) {
    return carray_set_from_narray(a, ary);
#endif
  } else {
    rb_raise(rb_eTypeError,
	     "wrong argument type %s", rb_class2name(CLASS_OF(ary)));
  }
}

void carray_set_from_rarray(double *a, VALUE ary)
{
  size_t i, size;
  VALUE val;
  if (CLASS_OF(ary) == rb_cRange) ary = rb_gsl_range2ary(ary);
  Check_Type(ary, T_ARRAY);
  //  size = RARRAY(ary)->len;
  size = RARRAY_LEN(ary);
  if (size == 0) return;
  for (i = 0; i < size; i++) {
    val = rb_ary_entry(ary, i);
    Need_Float(val);
    a[i] = NUM2DBL(val);
    //    a[i] = RFLOAT(val)->value;
  }
}

#ifdef HAVE_NARRAY_H
/* NArray -> CArray */
void carray_set_from_narray(double *a, VALUE ary)
{
  int size;
  VALUE ary2;
  if (!NA_IsNArray(ary))
    rb_raise(rb_eTypeError,
             "wrong argument type %s", rb_class2name(CLASS_OF(ary)));
  size = NA_TOTAL(ary);
  if (size == 0) return;
  ary2 = na_change_type(ary, NA_DFLOAT);
  memcpy(a, NA_PTR_TYPE(ary2,double*), size*sizeof(double));
}

/* NArray -> GSL::Vector */
gsl_vector* make_cvector_from_narray(VALUE ary)
{
  gsl_vector *v = NULL;
  size_t size;
  VALUE ary2;
  if (!NA_IsNArray(ary))
    rb_raise(rb_eTypeError,
             "wrong argument type %s", rb_class2name(CLASS_OF(ary)));
  size = NA_TOTAL(ary);
  v = gsl_vector_alloc(size);
  if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
  ary2 = na_change_type(ary, NA_DFLOAT);
  memcpy(v->data, NA_PTR_TYPE(ary2,double*), size*sizeof(double));
  /*  cvector_set_from_narray(v, ary);*/
  return v;
}

#endif

gsl_vector_complex* make_vector_complex_clone(const gsl_vector_complex *v)
{
  gsl_vector_complex *vnew = NULL;
  vnew = gsl_vector_complex_alloc(v->size);
  if (vnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
  gsl_vector_complex_memcpy(vnew, v);
  return vnew;
}

gsl_matrix* make_matrix_clone(const gsl_matrix *m)
{
  gsl_matrix *mnew = NULL;
  mnew = gsl_matrix_alloc(m->size1, m->size2);
  if (mnew == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_alloc failed");
  gsl_matrix_memcpy(mnew, m);
  return mnew;
}

gsl_matrix_int* make_matrix_int_clone(const gsl_matrix_int *m)
{
  gsl_matrix_int *mnew = NULL;
  mnew = gsl_matrix_int_alloc(m->size1, m->size2);
  if (mnew == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_alloc failed");
  gsl_matrix_int_memcpy(mnew, m);
  return mnew;
}

VALUE make_matrix_clone2(VALUE vm)
{
  gsl_matrix *m = NULL, *mnew = NULL;
  Data_Get_Struct(vm, gsl_matrix, m);
  mnew = gsl_matrix_alloc(m->size1, m->size2);
  if (mnew == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_alloc failed");
  gsl_matrix_memcpy(mnew, m);
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);;
}

gsl_matrix_complex* make_matrix_complex_clone(const gsl_matrix_complex *m)
{
  gsl_matrix_complex *mnew = NULL;
  mnew = gsl_matrix_complex_alloc(m->size1, m->size2);
  if (mnew == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
  gsl_matrix_complex_memcpy(mnew, m);
  return mnew;
}

gsl_vector_complex* vector_to_complex(const gsl_vector *v)
{
  gsl_vector_complex *cv = NULL;
  gsl_complex z;
  size_t i;
  cv = gsl_vector_complex_alloc(v->size);
  if (cv == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
  for (i = 0; i < v->size; i++) {
    z = gsl_complex_rect(gsl_vector_get(v, i), 0.0);
    gsl_vector_complex_set(cv, i, z);
  }
  return cv;
}

gsl_matrix_complex* matrix_to_complex(const gsl_matrix *m)
{
  gsl_matrix_complex *cm = NULL;
  gsl_complex z;
  size_t i, j;
  cm = gsl_matrix_complex_alloc(m->size1, m->size2);
  if (cm == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      z = gsl_complex_rect(gsl_matrix_get(m, i, j), 0.0);
      gsl_matrix_complex_set(cm, i, j, z);
    }
  }
  return cm;
}

void gsl_matrix_complex_mul(gsl_matrix_complex *mnew, const gsl_matrix_complex *m, 
			    const gsl_matrix_complex *mb)
{
  gsl_complex a, b, c, sum;
  size_t i, j, k;
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      sum = gsl_complex_rect(0.0, 0.0);
      for (k = 0; k < m->size2; k++) {
	a = gsl_matrix_complex_get(m, j, k);
	b = gsl_matrix_complex_get(mb, k, i);
	c = gsl_complex_mul(a, b);
	sum = gsl_complex_add(sum, c);
      }
      gsl_matrix_complex_set(mnew, j, i, sum);
    }
  }
}

void gsl_matrix_mul_vector(gsl_vector *vnew, 
			   const gsl_matrix *m, const gsl_vector *v)
{
  size_t i, j;
  double val;
  for (i = 0; i < m->size1; i++) {
    val = 0;
    for (j = 0; j < m->size2; j++) 
      val += gsl_matrix_get(m, i, j)*gsl_vector_get(v, j);
    gsl_vector_set(vnew, i, val);
  }
}

void gsl_matrix_int_mul_vector(gsl_vector_int *vnew, 
			       const gsl_matrix_int *m, const gsl_vector_int *v)
{
  size_t i, j;
  int val;
  for (i = 0; i < m->size1; i++) {
    val = 0;
    for (j = 0; j < m->size2; j++) 
      val += gsl_matrix_int_get(m, i, j)*gsl_vector_int_get(v, j);
    gsl_vector_int_set(vnew, i, val);
  }
}

void gsl_matrix_complex_mul_vector(gsl_vector_complex *vnew, 
				   const gsl_matrix_complex *m, 
				   const gsl_vector_complex *v)
{
  gsl_complex a, b, c, sum;
  size_t i, j;
  for (i = 0; i < m->size1; i++) {
    sum = gsl_complex_rect(0.0, 0.0);
    for (j = 0; j < m->size2; j++) {
      a = gsl_matrix_complex_get(m, i, j);
      b = gsl_vector_complex_get(v, j);
      c = gsl_complex_mul(a, b);
      sum = gsl_complex_add(sum, c);
    }
    gsl_vector_complex_set(vnew, i, sum);
  }
  
}

/*****/
#ifndef GSL_1_12_LATER
int gsl_vector_complex_add(gsl_vector_complex *cv, const gsl_vector_complex *cv2)
{
  size_t i;
  gsl_complex a, b, c;
  for (i = 0; i < cv->size; i++) {
    a = gsl_vector_complex_get(cv, i);
    b = gsl_vector_complex_get(cv2, i);
    c = gsl_complex_add(a, b);
    gsl_vector_complex_set(cv, i, c);
  }
  return 0;
}

int gsl_vector_complex_add_constant(gsl_vector_complex *cv, gsl_complex b)
{
  size_t i;
  gsl_complex a, c;
  for (i = 0; i < cv->size; i++) {
    a = gsl_vector_complex_get(cv, i);
    c = gsl_complex_add(a, b);
    gsl_vector_complex_set(cv, i, c);
  }
  return 0;
}

int gsl_vector_complex_scale(gsl_vector_complex *cv, gsl_complex b)
{
  size_t i;
  gsl_complex a, c;
  for (i = 0; i < cv->size; i++) {
    a = gsl_vector_complex_get(cv, i);
    c = gsl_complex_mul(a, b);
    gsl_vector_complex_set(cv, i, c);
  }
  return 0;
}


int gsl_vector_complex_sub(gsl_vector_complex *cv, const gsl_vector_complex *cv2)
{
  size_t i;
  gsl_complex a, b, c;
  for (i = 0; i < cv->size; i++) {
    a = gsl_vector_complex_get(cv, i);
    b = gsl_vector_complex_get(cv2, i);
    c = gsl_complex_sub(a, b);
    gsl_vector_complex_set(cv, i, c);
  }
  return 0;
}

int gsl_vector_complex_mul(gsl_vector_complex *cv, const gsl_vector_complex *cv2)
{
  size_t i;
  gsl_complex a, b, c;
  for (i = 0; i < cv->size; i++) {
    a = gsl_vector_complex_get(cv, i);
    b = gsl_vector_complex_get(cv2, i);
    c = gsl_complex_mul(a, b);
    gsl_vector_complex_set(cv, i, c);
  }
  return 0;
}

int gsl_vector_complex_div(gsl_vector_complex *cv, const gsl_vector_complex *cv2)
{
  size_t i;
  gsl_complex a, b, c;
  for (i = 0; i < cv->size; i++) {
    a = gsl_vector_complex_get(cv, i);
    b = gsl_vector_complex_get(cv2, i);
    c = gsl_complex_div(a, b);
    gsl_vector_complex_set(cv, i, c);
  }
  return 0;
}
#endif

VALUE rb_gsl_range2ary(VALUE obj)
{
  //  double beg, en;
  //  size_t n;
  //  int step;
  VALUE ary;
  if (CLASS_OF(obj) != rb_cRange) 
    rb_raise(rb_eTypeError, "wrong argument type %s (Range expected)",
	     rb_class2name(CLASS_OF(obj)));
  ary = rb_funcall(obj, rb_gsl_id_to_a, 0);
  return ary;
}

void get_range_beg_en_n(VALUE range, double *beg, double *en, size_t *n, int *step);
VALUE rb_gsl_range2vector(VALUE obj)
{
  double beg, en;
  size_t n;
  int i, step;
  gsl_vector *v;
  if (CLASS_OF(obj) != rb_cRange) 
    rb_raise(rb_eTypeError, "wrong argument type %s (Range expected)",
	     rb_class2name(CLASS_OF(obj)));
  get_range_beg_en_n(obj, &beg, &en, &n, &step);
  v = gsl_vector_alloc(n);
  for (i = 0; i < (int) n; i++) gsl_vector_set(v, i, (double) (beg+i));
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

void get_range_int_beg_en_n(VALUE range, int *beg, int *en, size_t *n, int *step);
VALUE rb_gsl_range2vector_int(VALUE obj)
{
  int beg, en, i, step;
  size_t n;
  gsl_vector_int *v;
  if (CLASS_OF(obj) != rb_cRange) 
    rb_raise(rb_eTypeError, "wrong argument type %s (Range expected)",
	     rb_class2name(CLASS_OF(obj)));
  get_range_int_beg_en_n(obj, &beg, &en, &n, &step);
  v = gsl_vector_int_alloc(n);
  for (i = 0; i < (int) n; i++) gsl_vector_int_set(v, i, beg+i);
  return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, v);
}

gsl_vector_int_view* rb_gsl_vector_int_view_alloc(size_t n)
{
  gsl_vector_int_view *v;
  v = ALLOC(gsl_vector_int_view);
  v->vector.size = n;
  v->vector.stride = 1;
  v->vector.owner = 0;
  return v;
}

void rb_gsl_vector_int_view_free(gsl_vector_int_view *v)
{
  free((gsl_vector_int_view *) v);
}

gsl_matrix_view* rb_gsl_matrix_view_alloc()
{
  gsl_matrix_view *mv = NULL;
  mv = ALLOC(gsl_matrix_view);
  if (mv == NULL) rb_raise(rb_eRuntimeError, "malloc failed");
  return mv;
}

void rb_gsl_matrix_view_free(gsl_matrix_view * mv)
{
  free((gsl_matrix_view *) mv);
}

gsl_matrix_int_view* rb_gsl_matrix_int_view_alloc()
{
  gsl_matrix_int_view *mv = NULL;
  mv = ALLOC(gsl_matrix_int_view);
  if (mv == NULL) rb_raise(rb_eRuntimeError, "malloc failed");
  return mv;
}

void rb_gsl_matrix_int_view_free(gsl_matrix_int_view * mv)
{
  free((gsl_matrix_int_view *) mv);
}

void Init_gsl_array(VALUE module)
{
  cgsl_block = rb_define_class_under(module, "Block", 
					 cGSL_Object);
  cgsl_block_int = rb_define_class_under(cgsl_block, "Int", 
					 cGSL_Object);
  cgsl_block_uchar = rb_define_class_under(cgsl_block, "Byte", 
					   cGSL_Object);

  cgsl_block_complex = rb_define_class_under(cgsl_block, "Complex", cgsl_block);
  cgsl_vector = rb_define_class_under(module, "Vector", 
				      cGSL_Object);
  cgsl_vector_col = rb_define_class_under(cgsl_vector, "Col", 
					 cgsl_vector);
  cgsl_vector_complex = rb_define_class_under(cgsl_vector, "Complex", 
					      cGSL_Object);
  cgsl_vector_complex_col = rb_define_class_under(cgsl_vector_complex, "Col",
					      cgsl_vector_complex);
  cgsl_matrix = rb_define_class_under(module, "Matrix", cGSL_Object);
  cgsl_matrix_complex = rb_define_class_under(cgsl_matrix, "Complex", cGSL_Object);

  cgsl_vector_view = rb_define_class_under(cgsl_vector, "View", cgsl_vector);
  cgsl_vector_col_view = rb_define_class_under(cgsl_vector_col, "View", cgsl_vector_col);

  cgsl_vector_complex_view = rb_define_class_under(cgsl_vector_complex, "View", 
						   cgsl_vector_complex);
  cgsl_vector_complex_col_view = rb_define_class_under(cgsl_vector_complex_col, "View", 
						   cgsl_vector_complex_col);

  cgsl_vector_int = rb_define_class_under(cgsl_vector, "Int", cGSL_Object);
  cgsl_vector_int_col = rb_define_class_under(cgsl_vector_int, "Col", cgsl_vector_int);  
  cgsl_vector_int_view = rb_define_class_under(cgsl_vector_int, "View", cgsl_vector_int);
  cgsl_vector_int_col_view = rb_define_class_under(cgsl_vector_int_col, "View", cgsl_vector_int_col);


  /*****/

  cgsl_matrix_view = rb_define_class_under(cgsl_matrix, "View", 
					   cgsl_matrix);
  cgsl_matrix_complex_view = rb_define_class_under(cgsl_matrix_complex, "View", 
						   cgsl_matrix_complex);
  cgsl_permutation = rb_define_class_under(module, "Permutation", cGSL_Object);  
  cgsl_index = rb_define_class_under(module, "Index", cgsl_permutation);

  cgsl_vector_view_ro = rb_define_class_under(cgsl_vector_view, "ReadOnly",
					      cgsl_vector_view);
  cgsl_vector_col_view_ro = rb_define_class_under(cgsl_vector_col_view, "ReadOnly",
					      cgsl_vector_col_view);
  cgsl_vector_int_view_ro = rb_define_class_under(cgsl_vector_int_view, "ReadOnly",
					      cgsl_vector_int_view);
  cgsl_vector_int_col_view_ro = rb_define_class_under(cgsl_vector_int_col_view, "ReadOnly",
					      cgsl_vector_int_col_view);
  cgsl_matrix_view_ro = rb_define_class_under(cgsl_matrix_view, "ReadOnly",
					      cgsl_matrix_view);

  cgsl_vector_complex_view_ro = rb_define_class_under(cgsl_vector_complex_view, 
						      "ReadOnly",
						      cgsl_vector_complex_view);
  cgsl_matrix_complex_view_ro = rb_define_class_under(cgsl_matrix_complex_view, 
						      "ReadOnly",
						      cgsl_matrix_complex_view);

  /*****/
  cgsl_matrix_int = rb_define_class_under(cgsl_matrix, "Int", cGSL_Object);
  cgsl_matrix_int_view = rb_define_class_under(cgsl_matrix_int, "View", cgsl_matrix_int);
  cgsl_matrix_int_view_ro = rb_define_class_under(cgsl_matrix_int_view, "ReadOnly",
					      cgsl_matrix_int_view);
  /*****/
  Init_gsl_block_init(module);
  Init_gsl_block_int_init(module);
  Init_gsl_block_uchar_init(module);
  Init_gsl_vector(module);
  Init_gsl_vector_int(module);

  Init_gsl_vector_complex(module);
  Init_gsl_matrix(module);
  Init_gsl_matrix_int(module);
  Init_gsl_matrix_complex(module);
  Init_gsl_permutation(module);
#ifdef GSL_1_1_LATER
  Init_gsl_combination(module);
#endif
  Init_gsl_array_complex(module);
  Init_gsl_matrix_nmf();

  rb_define_method(cgsl_vector_view_ro, "set", rb_gsl_obj_read_only, -1);
  rb_define_method(cgsl_matrix_view_ro, "set", rb_gsl_obj_read_only, -1);
  rb_define_method(cgsl_vector_int_view_ro, "set", rb_gsl_obj_read_only, -1);
  rb_define_method(cgsl_matrix_int_view_ro, "set", rb_gsl_obj_read_only, -1);  
  rb_define_method(cgsl_vector_complex_view_ro, "set", rb_gsl_obj_read_only, -1);
  rb_define_method(cgsl_matrix_complex_view_ro, "set", rb_gsl_obj_read_only, -1);

}
