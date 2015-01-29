/*
  fit.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_fit.h"

/* linear fit without weights: y = c0 + c1 x */
/* This returns 7 elements array */
static VALUE rb_gsl_fit_linear(int argc, VALUE *argv, VALUE obj)
{
  double *ptrx, *ptry;
  double c0, c1, cov00, cov01, cov11, sumsq;
  int status;
  size_t n, stridex, stridey;
  switch (argc) {
  case 2:
    ptrx = get_vector_ptr(argv[0], &stridex, &n);
    ptry = get_vector_ptr(argv[1], &stridey, &n);
    break;
  case 3:
    CHECK_FIXNUM(argv[2]);
    ptrx = get_vector_ptr(argv[0], &stridex, &n);
    ptry = get_vector_ptr(argv[1], &stridey, &n);
    n = FIX2INT(argv[2]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  status = gsl_fit_linear(ptrx, stridex, ptry, stridey, n, &c0, &c1, &cov00,
        &cov01, &cov11, &sumsq);
  return rb_ary_new3(7, rb_float_new(c0), rb_float_new(c1), rb_float_new(cov00),
         rb_float_new(cov01), rb_float_new(cov11), rb_float_new(sumsq),
         INT2FIX(status));
}

/* linear fit with weights: y = c0 + c1 x */
static VALUE rb_gsl_fit_wlinear(int argc, VALUE *argv, VALUE obj)
{
  double *ptrx, *ptry, *ptrw;
  double c0, c1, cov00, cov01, cov11, sumsq;
  int status;
  size_t n, stridex, stridey, stridew;
  switch (argc) {
  case 3:
    ptrx = get_vector_ptr(argv[0], &stridex, &n);
    ptrw = get_vector_ptr(argv[1], &stridew, &n);
    ptry = get_vector_ptr(argv[2], &stridey, &n);
    break;
  case 4:   
    CHECK_FIXNUM(argv[3]);
    ptrx = get_vector_ptr(argv[0], &stridex, &n);
    ptrw = get_vector_ptr(argv[1], &stridew, &n);
    ptry = get_vector_ptr(argv[2], &stridey, &n);
    n = FIX2INT(argv[3]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  status = gsl_fit_wlinear(ptrx, stridex, ptrw, stridew, ptry, stridey,
         n,
         &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
  return rb_ary_new3(7, rb_float_new(c0), rb_float_new(c1), rb_float_new(cov00),
         rb_float_new(cov01), rb_float_new(cov11), rb_float_new(sumsq),
         INT2FIX(status));
}

static VALUE rb_gsl_fit_linear_est(int argc, VALUE *argv, VALUE obj)
{
  double y, yerr, x, c0, c1, c00, c01, c11;
  int status;
  size_t i;
  switch (argc) {
  case 2:
    x = NUM2DBL(argv[0]);
    if (TYPE(argv[1]) == T_ARRAY) {
      c0 = NUM2DBL(rb_ary_entry(argv[1], 0));
      c1 = NUM2DBL(rb_ary_entry(argv[1], 1));
      c00 = NUM2DBL(rb_ary_entry(argv[1], 2));
      c01 = NUM2DBL(rb_ary_entry(argv[1], 3));
      c11 = NUM2DBL(rb_ary_entry(argv[1], 4));
    } else {
      rb_raise(rb_eTypeError, "argv[1] Array expected");
    }
    break;
  case 6:
    for (i = 0; i < 6; i++) Need_Float(argv[i]);
    x = NUM2DBL(argv[0]);
    c0 = NUM2DBL(argv[1]);
    c1 = NUM2DBL(argv[2]);
    c00 = NUM2DBL(argv[3]);
    c01 = NUM2DBL(argv[4]);
    c11 = NUM2DBL(argv[5]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 6)", argc);
  }
  status = gsl_fit_linear_est(x, c0, c1, c00, c01, c11, &y, &yerr);
  return rb_ary_new3(3, rb_float_new(y), rb_float_new(yerr), INT2FIX(status));
}

static VALUE rb_gsl_fit_mul(int argc, VALUE *argv, VALUE obj)
{
  double *ptrx, *ptry;
  double c1, cov11, sumsq;
  int status;
  size_t n, stridex, stridey;
  switch (argc) {
  case 2:
    ptrx = get_vector_ptr(argv[0], &stridex, &n);
    ptry = get_vector_ptr(argv[1], &stridey, &n);
    break;
  case 3:
    CHECK_FIXNUM(argv[2]);
    ptrx = get_vector_ptr(argv[0], &stridex, &n);
    ptry = get_vector_ptr(argv[1], &stridey, &n);
    n = FIX2INT(argv[2]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  status = gsl_fit_mul(ptrx, stridex, ptry, stridey, n, &c1, &cov11, &sumsq);
  return rb_ary_new3(4, rb_float_new(c1), 
         rb_float_new(cov11), rb_float_new(sumsq), INT2FIX(status));
}

static VALUE rb_gsl_fit_wmul(int argc, VALUE *argv, VALUE obj)
{
  double *ptrx, *ptry, *ptrw;
  double c1, cov11, sumsq;
  int status;
  size_t n, stridex, stridey, stridew;
  switch (argc) {
  case 3:
    ptrx = get_vector_ptr(argv[0], &stridex, &n);
    ptrw = get_vector_ptr(argv[1], &stridew, &n);
    ptry = get_vector_ptr(argv[2], &stridey, &n);
    break;
  case 4:   
    CHECK_FIXNUM(argv[3]);
    ptrx = get_vector_ptr(argv[0], &stridex, &n);
    ptrw = get_vector_ptr(argv[1], &stridew, &n);
    ptry = get_vector_ptr(argv[2], &stridey, &n);
    n = FIX2INT(argv[3]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  status = gsl_fit_wmul(ptrx, stridex, ptrw, stridew, ptry, stridey, 
      n, &c1, &cov11, &sumsq);
  return rb_ary_new3(4, rb_float_new(c1), 
         rb_float_new(cov11), rb_float_new(sumsq), INT2FIX(status));
}

static VALUE rb_gsl_fit_mul_est(int argc, VALUE *argv, VALUE obj)
{
  double y, yerr, x, c1, c11;
  int status;
  switch (argc) {
  case 2:
    Need_Float(argv[0]);
    if (TYPE(argv[1]) == T_ARRAY) {
      c1 = NUM2DBL(rb_ary_entry(argv[1], 0));
      c11 = NUM2DBL(rb_ary_entry(argv[1], 1));
    } else {
      rb_raise(rb_eTypeError, "argv[1]: Array expected");
    }
    x = NUM2DBL(argv[0]);
    break;
  case 3:
    Need_Float(argv[0]); Need_Float(argv[1]); Need_Float(argv[2]);
    x = NUM2DBL(argv[0]);
    c1 = NUM2DBL(argv[1]);
    c11 = NUM2DBL(argv[2]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  status = gsl_fit_mul_est(x, c1, c11, &y, &yerr);
  return rb_ary_new3(3, rb_float_new(y), rb_float_new(yerr), INT2FIX(status));
}

void Init_gsl_fit(VALUE module)
{
  VALUE mgsl_fit;
  mgsl_fit = rb_define_module_under(module, "Fit"); 
  rb_define_module_function(mgsl_fit, "linear", rb_gsl_fit_linear, -1);
  rb_define_module_function(mgsl_fit, "wlinear", rb_gsl_fit_wlinear, -1);
  rb_define_module_function(mgsl_fit, "linear_est", rb_gsl_fit_linear_est, -1);
  rb_define_module_function(mgsl_fit, "mul", rb_gsl_fit_mul, -1);
  rb_define_module_function(mgsl_fit, "wmul", rb_gsl_fit_wmul, -1);
  rb_define_module_function(mgsl_fit, "mul_est", rb_gsl_fit_mul_est, -1);
}
