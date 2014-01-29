/*
  sum.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_array.h"
#include <gsl/gsl_sum.h>

static VALUE rb_gsl_sum_accel(VALUE obj)
{
  gsl_sum_levin_u_workspace *w = NULL;
  double sum, err, sum_plain, *ptr;
  size_t terms_used, n, stride;  
  ptr = get_vector_ptr(obj, &stride, &n);
  w = gsl_sum_levin_u_alloc(n);
  gsl_sum_levin_u_accel(ptr, n, w, &sum, &err);
  sum_plain = w->sum_plain;
  terms_used = w->terms_used;
  gsl_sum_levin_u_free(w);
 return rb_ary_new3(4, rb_float_new(sum), rb_float_new(err), 
		    rb_float_new(sum_plain), INT2FIX(terms_used));
}
 
static VALUE rb_gsl_utrunc_accel(VALUE obj)
{
  gsl_sum_levin_utrunc_workspace *w = NULL;
  double sum, err, sum_plain, *ptr;
  size_t terms_used, n, stride;
  ptr = get_vector_ptr(obj, &stride, &n);
  w = gsl_sum_levin_utrunc_alloc(n);
  gsl_sum_levin_utrunc_accel(ptr, n, w, &sum, &err);
  sum_plain = w->sum_plain;
  terms_used = w->terms_used;
  gsl_sum_levin_utrunc_free(w);
  return rb_ary_new3(4, rb_float_new(sum), rb_float_new(err), 
		     rb_float_new(sum_plain), INT2FIX(terms_used));
}

static VALUE rb_gsl_sum_levin_u_new(VALUE klass, VALUE nn)
{
  gsl_sum_levin_u_workspace *w = NULL;
  CHECK_FIXNUM(nn);
  w = gsl_sum_levin_u_alloc(FIX2INT(nn));
  return Data_Wrap_Struct(klass, 0, gsl_sum_levin_u_free, w);
}

static VALUE rb_gsl_sum_levin_utrunc_new(VALUE klass, VALUE nn)
{
  gsl_sum_levin_utrunc_workspace *w = NULL;
  CHECK_FIXNUM(nn);
  w = gsl_sum_levin_utrunc_alloc(FIX2INT(nn));
  return Data_Wrap_Struct(klass, 0, gsl_sum_levin_utrunc_free, w);
}

/* singleton */
static VALUE rb_gsl_sum_levin_u_accel2(VALUE obj, VALUE vv)
{
  gsl_sum_levin_u_workspace *w = NULL;
  double sum, err, sum_plain;
  size_t terms_used, n, stride;
  double *ptr;
  ptr = get_vector_ptr(vv, &stride, &n);
  w = gsl_sum_levin_u_alloc(n);
  gsl_sum_levin_u_accel(ptr, n, w, &sum, &err);
  sum_plain = w->sum_plain;
  terms_used = w->terms_used;
  gsl_sum_levin_u_free(w);
  return rb_ary_new3(4, rb_float_new(sum), rb_float_new(err), 
		     rb_float_new(sum_plain), INT2FIX(terms_used));
}

static VALUE rb_gsl_sum_levin_utrunc_accel2(VALUE obj, VALUE vv)
{
  gsl_sum_levin_utrunc_workspace *w = NULL;
  double sum, err, sum_plain;
  size_t terms_used, n, stride;
  double *ptr;
  ptr = get_vector_ptr(vv, &stride, &n);
  w = gsl_sum_levin_utrunc_alloc(n);
  gsl_sum_levin_utrunc_accel(ptr, n, w, &sum, &err);
  sum_plain = w->sum_plain;
  terms_used = w->terms_used;
  gsl_sum_levin_utrunc_free(w);
  return rb_ary_new3(4, rb_float_new(sum), rb_float_new(err), 
		     rb_float_new(sum_plain), INT2FIX(terms_used));
}

static VALUE rb_gsl_sum_levin_u_sum_plain(VALUE obj)
{
  gsl_sum_levin_u_workspace *w = NULL;
  Data_Get_Struct(obj, gsl_sum_levin_u_workspace, w);
  return rb_float_new(w->sum_plain);
}

static VALUE rb_gsl_sum_levin_u_terms_used(VALUE obj)
{
  gsl_sum_levin_u_workspace *w = NULL;
  Data_Get_Struct(obj, gsl_sum_levin_u_workspace, w);
  return INT2FIX(w->terms_used);
}

static VALUE rb_gsl_sum_levin_utrunc_sum_plain(VALUE obj)
{
  gsl_sum_levin_utrunc_workspace *w = NULL;
  Data_Get_Struct(obj, gsl_sum_levin_utrunc_workspace, w);
  return rb_float_new(w->sum_plain);
}

static VALUE rb_gsl_sum_levin_utrunc_terms_used(VALUE obj)
{
  gsl_sum_levin_utrunc_workspace *w = NULL;
  Data_Get_Struct(obj, gsl_sum_levin_utrunc_workspace, w);
  return INT2FIX(w->terms_used);
}

void Init_gsl_sum(VALUE module) 
{
  VALUE mgsl_sum;
  VALUE cgsl_sum_levin_u, cgsl_sum_levin_utrunc;

  mgsl_sum = rb_define_module_under(module, "Sum");
  cgsl_sum_levin_u = rb_define_class_under(mgsl_sum, 
					   "Levin_u", cGSL_Object);
  cgsl_sum_levin_utrunc = rb_define_class_under(mgsl_sum, 
						"Levin_utrunc", cGSL_Object);

  rb_define_singleton_method(cgsl_sum_levin_u, "new", rb_gsl_sum_levin_u_new, 1);
  rb_define_singleton_method(cgsl_sum_levin_u, "alloc", rb_gsl_sum_levin_u_new, 1);
  rb_define_singleton_method(cgsl_sum_levin_utrunc, "new", 
			     rb_gsl_sum_levin_utrunc_new, 1);
  rb_define_singleton_method(cgsl_sum_levin_utrunc, "alloc", 
			     rb_gsl_sum_levin_utrunc_new, 1);
  rb_define_singleton_method(cgsl_sum_levin_u, "accel", 
			     rb_gsl_sum_levin_u_accel2, 1);

  rb_define_singleton_method(cgsl_sum_levin_utrunc, "accel", 
			     rb_gsl_sum_levin_utrunc_accel2, 1);
  rb_define_method(cgsl_sum_levin_u, "accel", rb_gsl_sum_levin_u_accel2, 1);
  rb_define_method(cgsl_sum_levin_utrunc, "accel", 
		   rb_gsl_sum_levin_utrunc_accel2, 1);

  rb_define_method(cgsl_sum_levin_u, "sum_plain", rb_gsl_sum_levin_u_sum_plain, 0);
  rb_define_method(cgsl_sum_levin_u, "terms_used", 
		   rb_gsl_sum_levin_u_terms_used, 0);
  rb_define_method(cgsl_sum_levin_utrunc, "sum_plain", 
		   rb_gsl_sum_levin_utrunc_sum_plain, 0);
  rb_define_method(cgsl_sum_levin_utrunc, "terms_used", 
		   rb_gsl_sum_levin_utrunc_terms_used, 0);
  /***/

  rb_define_method(cgsl_vector, "accel_sum", rb_gsl_sum_accel, 0);
  rb_define_alias(cgsl_vector, "accel", "accel_sum");
  rb_define_alias(cgsl_vector, "sum_accel", "accel_sum");
  rb_define_method(cgsl_vector, "utrunc_accel", rb_gsl_utrunc_accel, 0);

#ifdef HAVE_NARRAY_H
  rb_define_method(cNArray, "accel_sum", rb_gsl_sum_accel, 0);
  rb_define_alias(cNArray, "accel", "accel_sum");
  rb_define_alias(cNArray, "sum_accel", "accel_sum");
  rb_define_method(cNArray, "utrunc_accel", rb_gsl_utrunc_accel, 0);
#endif
}
