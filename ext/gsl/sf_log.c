/*
  sf_log.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_sf.h"

static VALUE rb_gsl_sf_log(VALUE obj, VALUE x)
{
  if (COMPLEX_P(x) || VECTOR_COMPLEX_P(x) || MATRIX_COMPLEX_P(x)) 
    return rb_gsl_math_complex_eval(gsl_complex_log, x);
  return rb_gsl_sf_eval1(gsl_sf_log, x);
}

static VALUE rb_gsl_sf_log10(VALUE obj, VALUE x)
{
  if (COMPLEX_P(x) || VECTOR_COMPLEX_P(x) || MATRIX_COMPLEX_P(x)) 
    return rb_gsl_math_complex_eval(gsl_complex_log10, x);
  return rb_gsl_sf_eval1(log10, x);
}


static VALUE rb_gsl_sf_log_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_log_e, x);
}

static VALUE rb_gsl_sf_log_abs(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_log_abs, x);
}

static VALUE rb_gsl_sf_log_abs_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_log_abs_e, x);
}

static VALUE rb_gsl_sf_complex_log_e(int argc, VALUE *argv, VALUE obj)
{
  gsl_sf_result *rslt1 = NULL, *rslt2 = NULL;
  gsl_complex *z = NULL;
  VALUE vlnr, vtheta;
  double re, im;
  // local variable "status" was defined and set, but never used
  //int status;
  switch (argc) {
  case 1:
    CHECK_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_complex, z);
    re = GSL_REAL(*z);
    im = GSL_IMAG(*z);
    break;
  case 2:
    Need_Float(argv[0]);     Need_Float(argv[1]);
    re = NUM2DBL(argv[0]);
    im = NUM2DBL(argv[1]);
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  vlnr = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt1);
  vtheta = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt2);
  /*status =*/ gsl_sf_complex_log_e(re, im, rslt1, rslt2);
  return rb_ary_new3(2, vlnr, vtheta);
}

static VALUE rb_gsl_sf_log_1plusx(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_log_1plusx, x);
}

static VALUE rb_gsl_sf_log_1plusx_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_log_1plusx_e, x);
}

static VALUE rb_gsl_sf_log_1plusx_mx(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_log_1plusx_mx, x);
}

static VALUE rb_gsl_sf_log_1plusx_mx_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_log_1plusx_mx_e, x);
}

void Init_gsl_sf_log(VALUE module)
{
  rb_define_module_function(module, "log",  rb_gsl_sf_log, 1);
  rb_define_module_function(module, "log10",  rb_gsl_sf_log10, 1);  
  rb_define_module_function(module, "log_e",  rb_gsl_sf_log_e, 1);
  rb_define_module_function(module, "log_abs",  rb_gsl_sf_log_abs, 1);
  rb_define_module_function(module, "log_abs_e",  rb_gsl_sf_log_abs_e, 1);
  rb_define_module_function(module, "complex_log_e",  rb_gsl_sf_complex_log_e, -1);
  rb_define_module_function(module, "log_1plusx",  rb_gsl_sf_log_1plusx, 1);
  rb_define_module_function(module, "log_1plusx_e",  rb_gsl_sf_log_1plusx_e, 1);
  rb_define_module_function(module, "log_1plusx_mx",  rb_gsl_sf_log_1plusx_mx, 1);
  rb_define_module_function(module, "log_1plusx_mx_e",  rb_gsl_sf_log_1plusx_mx_e, 1);
}
