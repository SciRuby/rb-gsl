/*
  sort.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_config.h"
#include "rb_gsl_array.h"
#include <gsl/gsl_heapsort.h>
#include <gsl/gsl_sort.h>
EXTERN ID RBGSL_ID_call;
EXTERN VALUE cgsl_complex;

int rb_gsl_comparison_double(const void *aa, const void *bb);
int rb_gsl_comparison_complex(const void *aa, const void *bb);
int rb_gsl_comparison_double(const void *aa, const void *bb)
{
  double *a = NULL, *b = NULL;
  a = (double *) aa;
  b = (double *) bb;
  return FIX2INT(rb_funcall(rb_block_proc(), RBGSL_ID_call, 2, rb_float_new(*a), rb_float_new(*b)));
}

int rb_gsl_comparison_complex(const void *aa, const void *bb)
{
  gsl_complex *a = NULL, *b = NULL;
  a = (gsl_complex *) aa;
  b = (gsl_complex *) bb;
  return FIX2INT(rb_funcall(rb_block_proc(), RBGSL_ID_call, 2, 
			    Data_Wrap_Struct(cgsl_complex, 0, NULL, a),
			    Data_Wrap_Struct(cgsl_complex, 0, NULL, b)));
}

static VALUE rb_gsl_heapsort_vector(VALUE obj)
{
  gsl_vector *v = NULL;
  if (!rb_block_given_p()) rb_raise(rb_eRuntimeError, "Proc is not given");
  Data_Get_Struct(obj, gsl_vector, v);
  gsl_heapsort(v->data, v->size, sizeof(double), rb_gsl_comparison_double);
  return obj;
}

static VALUE rb_gsl_heapsort_vector2(VALUE obj)
{
  gsl_vector *v = NULL, *vnew = NULL;
  if (!rb_block_given_p()) rb_raise(rb_eRuntimeError, "Proc is not given");
  Data_Get_Struct(obj, gsl_vector, v);
  vnew = gsl_vector_alloc(v->size);
  gsl_vector_memcpy(vnew, v);
  gsl_heapsort(vnew->data, vnew->size, sizeof(double), rb_gsl_comparison_double);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
}

static VALUE rb_gsl_heapsort_index_vector(VALUE obj)
{
  gsl_vector *v = NULL;
  gsl_permutation *p = NULL;
  if (!rb_block_given_p()) rb_raise(rb_eRuntimeError, "Proc is not given");
  Data_Get_Struct(obj, gsl_vector, v);
  p = gsl_permutation_alloc(v->size);
  gsl_heapsort_index(p->data, v->data, v->size, sizeof(double), rb_gsl_comparison_double);
  return Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p);
}

static VALUE rb_gsl_heapsort_vector_complex(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  if (!rb_block_given_p()) rb_raise(rb_eRuntimeError, "Proc is not given");
  Data_Get_Struct(obj, gsl_vector_complex, v);
  gsl_heapsort(v->data, v->size, sizeof(gsl_complex), rb_gsl_comparison_complex);
  return obj;
}

static VALUE rb_gsl_heapsort_vector_complex2(VALUE obj)
{
  gsl_vector_complex *v = NULL, *vnew = NULL;
  if (!rb_block_given_p()) rb_raise(rb_eRuntimeError, "Proc is not given");
  Data_Get_Struct(obj, gsl_vector_complex, v);
  vnew = gsl_vector_complex_alloc(v->size);
  gsl_vector_complex_memcpy(vnew, v);
  gsl_heapsort(vnew->data, vnew->size, sizeof(gsl_complex), rb_gsl_comparison_complex);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
}

static VALUE rb_gsl_heapsort_index_vector_complex(VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_permutation *p = NULL;
  if (!rb_block_given_p()) rb_raise(rb_eRuntimeError, "Proc is not given");
  Data_Get_Struct(obj, gsl_vector_complex, v);
  p = gsl_permutation_alloc(v->size);
  gsl_heapsort_index(p->data, v->data, v->size, sizeof(gsl_complex), rb_gsl_comparison_complex);
  return Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p);
}

/* singleton */
static VALUE rb_gsl_heapsort(VALUE obj, VALUE vv)
{
  if (!rb_block_given_p()) rb_raise(rb_eRuntimeError, "Proc is not given");
  if (rb_obj_is_kind_of(vv, cgsl_vector_complex)) {
    return rb_gsl_heapsort_vector_complex(vv);
  } else if (rb_obj_is_kind_of(vv, cgsl_vector)) {
    return rb_gsl_heapsort_vector(vv);
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (Vector or Vector::Complex expected)", rb_class2name(CLASS_OF(vv)));
  }
  return vv;
}

static VALUE rb_gsl_heapsort2(VALUE obj, VALUE vv)
{
  if (!rb_block_given_p()) rb_raise(rb_eRuntimeError, "Proc is not given");
  if (rb_obj_is_kind_of(vv, cgsl_vector_complex)) {
    return rb_gsl_heapsort_vector_complex2(vv);
  } else if (rb_obj_is_kind_of(vv, cgsl_vector)) {
    return rb_gsl_heapsort_vector2(vv);
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (Vector or Vector::Complex expected)", rb_class2name(CLASS_OF(vv)));
  }
  return vv;
}

static VALUE rb_gsl_heapsort_index(VALUE obj, VALUE vv)
{
  if (!rb_block_given_p()) rb_raise(rb_eRuntimeError, "Proc is not given");
  if (rb_obj_is_kind_of(vv, cgsl_vector_complex)) {
    return rb_gsl_heapsort_index_vector_complex(vv);
  } else if (rb_obj_is_kind_of(vv, cgsl_vector)) {
    return rb_gsl_heapsort_index_vector(vv);
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (Vector or Vector::Complex expected)", rb_class2name(CLASS_OF(vv)));
  }
  return vv;
}

/*****/

#ifdef HAVE_NARRAY_H
#include "narray.h"
static VALUE rb_gsl_sort_narray(VALUE obj)
{
  struct NARRAY *na;
  size_t size, stride;
  double *ptr1, *ptr2;
  VALUE ary;
  GetNArray(obj, na);
  ptr1 = (double*) na->ptr;
  size = na->total;
  stride = 1;
  ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(obj));
  ptr2 = NA_PTR_TYPE(ary, double*);
  memcpy(ptr2, ptr1, sizeof(double)*size);
  gsl_sort(ptr2, stride, size);
  return ary;
}
static VALUE rb_gsl_sort_narray_bang(VALUE obj)
{
  struct NARRAY *na;
  size_t size, stride;
  double *ptr1;
  GetNArray(obj, na);
  ptr1 = (double*) na->ptr;
  size = na->total;
  stride = 1;
  gsl_sort(ptr1, stride, size);
  return obj;
}
static VALUE rb_gsl_sort_index_narray(VALUE obj)
{
  struct NARRAY *na;
  size_t size, stride;
  double *ptr1;
  gsl_permutation *p;
  GetNArray(obj, na);
  ptr1 = (double*) na->ptr;
  size = na->total;
  stride = 1;
  p = gsl_permutation_alloc(size);
  gsl_sort_index(p->data, ptr1, stride, size);
  return Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p);
}
#endif

void Init_gsl_sort(VALUE module)
{
  rb_define_singleton_method(module, "heapsort!", rb_gsl_heapsort, 1);
  rb_define_singleton_method(module, "heapsort", rb_gsl_heapsort2, 1);
  rb_define_singleton_method(module, "heapsort_index", rb_gsl_heapsort_index, 1);

  rb_define_method(cgsl_vector, "heapsort!", rb_gsl_heapsort_vector, 0);
  rb_define_method(cgsl_vector, "heapsort", rb_gsl_heapsort_vector2, 0);
  rb_define_method(cgsl_vector, "heapsort_index", rb_gsl_heapsort_index_vector, 0);

  rb_define_method(cgsl_vector_complex, "heapsort!", rb_gsl_heapsort_vector_complex, 0);
  rb_define_method(cgsl_vector_complex, "heapsort", rb_gsl_heapsort_vector_complex2, 0);
  rb_define_method(cgsl_vector_complex, "heapsort_index", rb_gsl_heapsort_index_vector_complex, 0);

#ifdef HAVE_NARRAY_H
  rb_define_method(cNArray, "gsl_sort", rb_gsl_sort_narray, 0);
  rb_define_method(cNArray, "gsl_sort!", rb_gsl_sort_narray_bang, 0);
  rb_define_method(cNArray, "gsl_sort_index", rb_gsl_sort_index_narray, 0);
#endif
}
