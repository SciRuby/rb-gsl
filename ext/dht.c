/*
  dht.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl.h"
#include "rb_gsl_array.h"
#include "rb_gsl_common.h"
#include <gsl/gsl_dht.h>
#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif

static VALUE rb_gsl_dht_alloc(int argc, VALUE *argv, VALUE klass)
{
  gsl_dht *t = NULL;
  switch (argc) {
  case 1:
    CHECK_FIXNUM(argv[0]);
    t = gsl_dht_alloc(FIX2INT(argv[0]));
    break;
  case 3:
    CHECK_FIXNUM(argv[0]);
    Need_Float(argv[1]); Need_Float(argv[2]);
    t = gsl_dht_new(FIX2INT(argv[0]), NUM2DBL(argv[1]), NUM2DBL(argv[2]));
   break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 3)", argc);
    break;
  }
  return Data_Wrap_Struct(klass, 0, gsl_dht_free, t);
}

static VALUE rb_gsl_dht_init(VALUE obj, VALUE nu, VALUE xmax)
{
  gsl_dht *t = NULL;
  Need_Float(nu); Need_Float(xmax);
  Data_Get_Struct(obj, gsl_dht, t);
  gsl_dht_init(t, NUM2DBL(nu), NUM2DBL(xmax));
  return obj;
}

static VALUE rb_gsl_dht_apply(int argc, VALUE *argv, VALUE obj)
{
  gsl_dht *t = NULL;
  double *ptr1, *ptr2;
  gsl_vector *vin, *vout;
  size_t size, stride;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
#endif
  VALUE ary;
  switch (argc) {
  case 2:
    Data_Get_Struct(obj, gsl_dht, t);
    ptr1 = get_vector_ptr(argv[0], &stride, &size);
    ptr2 = get_vector_ptr(argv[1], &stride, &size);
    return INT2FIX(gsl_dht_apply(t, ptr1, ptr2));
    break;
  case 1:
    Data_Get_Struct(obj, gsl_dht, t);
    if (VECTOR_P(argv[0])) {
      Data_Get_Struct(argv[0], gsl_vector, vin);
      ptr1 = vin->data;
      vout = gsl_vector_alloc(vin->size);
      ptr2 = vout->data;
      ary = Data_Wrap_Struct(VECTOR_ROW_COL(argv[0]), 0, gsl_vector_free, vout);
#ifdef HAVE_NARRAY_H
    } else if (NA_IsNArray(argv[0])) {
      GetNArray(argv[0], na);
      ptr1 = (double*)na->ptr;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv[0]));
      ptr2 = NA_PTR_TYPE(ary, double*);
#endif
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (Vector expected)",
	       rb_class2name(CLASS_OF(argv[0])));
    }
    gsl_dht_apply(t, ptr1, ptr2);
    return ary;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  return Qnil; /* never reach here */
}

static VALUE rb_gsl_dht_xk_sample(VALUE obj, VALUE n, 
			       double (*sample)(const gsl_dht*, int))
{
  gsl_dht *t = NULL;
  gsl_vector_int *vi;
  gsl_vector *v;
  size_t i, size;
  int nn;
  VALUE ary;
  double val;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  int *ptr;
  double *ptr2;
#endif
  Data_Get_Struct(obj, gsl_dht, t);
  if (CLASS_OF(n) == rb_cRange) n = rb_gsl_range2ary(n);
  switch (TYPE(n)) {
  case T_FIXNUM:
    return rb_float_new((*sample)(t, FIX2INT(n)));
    break;
  case T_ARRAY:
    //    size = RARRAY(n)->len;
    size = RARRAY_LEN(n);
    ary = rb_ary_new2(size);
    for (i = 0; i < size; i++) {
      nn = FIX2INT(rb_ary_entry(n, i));
      val = (*sample)(t, nn);
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
    if (VECTOR_INT_P(n)) {
      Data_Get_Struct(n, gsl_vector_int, vi);
      v = gsl_vector_alloc(vi->size);
      for (i = 0; i < v->size; i++) {
	nn = gsl_vector_int_get(vi, i);
	val = (*sample)(t, nn);
	gsl_vector_set(v, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
#ifdef HAVE_NARRAY_H
    } else if (NA_IsNArray(n)) {
      GetNArray(n, na);
      ptr = (int*) na->ptr;
      size = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, cNArray);
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < size; i++) {
	ptr2[i] = (*sample)(t, ptr[i]);	
      }
      return ary;
#endif
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (Vector::Int expected)",
	       rb_class2name(CLASS_OF(n)));
    }

  }
  return Qnil;
}

static VALUE rb_gsl_dht_x_sample(VALUE obj, VALUE n)
{
  return rb_gsl_dht_xk_sample(obj, n, gsl_dht_x_sample);
}

static VALUE rb_gsl_dht_k_sample(VALUE obj, VALUE n)
{
  return rb_gsl_dht_xk_sample(obj, n, gsl_dht_k_sample);
}

static VALUE rb_gsl_dht_j(VALUE obj)
{
  gsl_dht *t = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_dht, t);
  v = rb_gsl_make_vector_view(t->j, (t->size+2), 1);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, free, v);
}

static VALUE rb_gsl_dht_zero(VALUE obj)
{
  gsl_dht *t = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_dht, t);
  v = rb_gsl_make_vector_view(t->j+1, (t->size+1), 1);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, free, v);
}

static VALUE rb_gsl_dht_Jjj(VALUE obj)
{
  gsl_dht *t = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_dht, t);
  v = rb_gsl_make_vector_view(t->Jjj, t->size*(t->size+1)/2, 1);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, free, v);
}

static VALUE rb_gsl_dht_sample(int argc, VALUE *argv, VALUE obj)
{
  gsl_dht *t = NULL;
  gsl_matrix *mm = NULL;
  size_t n, m;
  double val;
  Data_Get_Struct(obj, gsl_dht, t);
  switch (argc) {
  case 0:
    mm = gsl_matrix_alloc(t->size, t->size);
    for (n = 0; n < t->size; n++) {
      for (m = 0; m < t->size; m++) {
	val =  t->j[n+1]*gsl_dht_x_sample(t, m)/t->xmax;
	gsl_matrix_set(mm, n, m, val);
      }
    }
    return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mm);
    break;
  case 2:
    n = FIX2INT(argv[0]);
    m = FIX2INT(argv[1]);
    val = t->j[n+1]*gsl_dht_x_sample(t, m)/t->xmax;
    return rb_float_new(val);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 2)", argc);
    break;
  }
  return Qnil;
}

static VALUE rb_gsl_dht_num(int argc, VALUE *argv, VALUE obj)
{
  gsl_dht *t = NULL;
  gsl_matrix *mm = NULL;
  size_t n, m;
  double val;
  Data_Get_Struct(obj, gsl_dht, t);
  switch (argc) {
  case 0:
    mm = gsl_matrix_alloc(t->size, t->size);
    for (n = 0; n < t->size; n++) {
      for (m = 0; m < t->size; m++) {
	val = gsl_sf_bessel_Jnu(t->nu, t->j[n+1]*gsl_dht_x_sample(t, m)/t->xmax);
	gsl_matrix_set(mm, n, m, val);
      }
    }
    return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mm);
    break;
  case 2:
    n = FIX2INT(argv[0]);
    m = FIX2INT(argv[1]);
    val = gsl_sf_bessel_Jnu(t->nu, t->j[n+1]*gsl_dht_x_sample(t, m)/t->xmax);
    return rb_float_new(val);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 2)", argc);
    break;
  }
  return Qnil;
}

static VALUE rb_gsl_dht_J2(VALUE obj)
{
  gsl_dht *t = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_dht, t);
  v = rb_gsl_make_vector_view(t->J2, t->size+1, 1);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, free, v);
}

static VALUE rb_gsl_dht_den(VALUE obj)
{
  gsl_dht *t = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_dht, t);
  v = rb_gsl_make_vector_view(t->J2+1, t->size, 1);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, free, v);
}

static VALUE rb_gsl_dht_size(VALUE obj)
{
  gsl_dht *t = NULL;
  Data_Get_Struct(obj, gsl_dht, t);
  return INT2FIX(t->size);
}

static VALUE rb_gsl_dht_nu(VALUE obj)
{
  gsl_dht *t = NULL;
  Data_Get_Struct(obj, gsl_dht, t);
  return rb_float_new(t->nu);
}

static VALUE rb_gsl_dht_xmax(VALUE obj)
{
  gsl_dht *t = NULL;
  Data_Get_Struct(obj, gsl_dht, t);
  return rb_float_new(t->xmax);
}

static VALUE rb_gsl_dht_kmax(VALUE obj)
{
  gsl_dht *t = NULL;
  Data_Get_Struct(obj, gsl_dht, t);
  return rb_float_new(t->kmax);
}

static VALUE rb_gsl_dht_coef(int argc, VALUE *argv, VALUE obj)
{
  gsl_dht *t = NULL;
  gsl_matrix *mm = NULL;
  size_t n, m;
  double val;
  Data_Get_Struct(obj, gsl_dht, t);
  switch (argc) {
  case 0:
    mm = gsl_matrix_alloc(t->size, t->size);
    for (n = 0; n < t->size; n++) {
      for (m = 0; m < t->size; m++) {
	val = gsl_sf_bessel_Jnu(t->nu, t->j[n+1]*gsl_dht_x_sample(t, m)/t->xmax);
	val *= (2.0/t->xmax/t->xmax)/t->J2[m+1];
	gsl_matrix_set(mm, n, m, val);
      }
    }
    return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mm);
    break;
  case 2:
    n = FIX2INT(argv[0]);
    m = FIX2INT(argv[1]);
    val = gsl_sf_bessel_Jnu(t->nu, t->j[n+1]*gsl_dht_x_sample(t, m)/t->xmax);
    val *= (2.0/t->xmax/t->xmax)/t->J2[m+1];
    return rb_float_new(val);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 2)", argc);
    break;
  }
  return Qnil;
}

void Init_gsl_dht(VALUE module)
{
  VALUE cgsl_dht;
  cgsl_dht = rb_define_class_under(module, "Dht", cGSL_Object);
  rb_define_singleton_method(cgsl_dht, "alloc", rb_gsl_dht_alloc, -1);
  rb_define_method(cgsl_dht, "init", rb_gsl_dht_init, 2);
  rb_define_method(cgsl_dht, "apply", rb_gsl_dht_apply, -1);
  rb_define_method(cgsl_dht, "x_sample", rb_gsl_dht_x_sample, 1);
  rb_define_method(cgsl_dht, "k_sample", rb_gsl_dht_k_sample, 1);

  rb_define_method(cgsl_dht, "size", rb_gsl_dht_size, 0);
  rb_define_method(cgsl_dht, "nu", rb_gsl_dht_nu, 0);
  rb_define_method(cgsl_dht, "xmax", rb_gsl_dht_xmax, 0);
  rb_define_method(cgsl_dht, "kmax", rb_gsl_dht_kmax, 0);

  rb_define_method(cgsl_dht, "j", rb_gsl_dht_j, 0);
  rb_define_method(cgsl_dht, "Jjj", rb_gsl_dht_Jjj, 0);
  rb_define_method(cgsl_dht, "J2", rb_gsl_dht_J2, 0);

  rb_define_method(cgsl_dht, "zero", rb_gsl_dht_zero, 0);
  rb_define_alias(cgsl_dht, "zeros", "zero");
  rb_define_method(cgsl_dht, "sample", rb_gsl_dht_sample, -1);
  rb_define_method(cgsl_dht, "num", rb_gsl_dht_num, -1);
  rb_define_method(cgsl_dht, "den", rb_gsl_dht_den, 0);
  rb_define_method(cgsl_dht, "coef", rb_gsl_dht_coef, -1);
}
