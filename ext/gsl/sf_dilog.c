/*
  sf_dilog.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "rb_gsl_sf.h"

static VALUE rb_gsl_sf_dilog(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_dilog, x);
}

static VALUE rb_gsl_sf_dilog_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_dilog_e, x);
}

static VALUE rb_gsl_sf_complex_dilog_e(VALUE obj, VALUE r, VALUE theta)
{
  gsl_sf_result *re, *im;
  VALUE vre, vim;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(r); Need_Float(theta);
  vre = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, re);
  vim = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, im);
  /*status =*/ gsl_sf_complex_dilog_e(NUM2DBL(r), NUM2DBL(theta), re, im);
  return rb_ary_new3(2, vre, vim);
}

void Init_gsl_sf_dilog(VALUE module)
{
  rb_define_module_function(module, "dilog",  rb_gsl_sf_dilog, 1);
  rb_define_module_function(module, "dilog_e",  rb_gsl_sf_dilog_e, 1);
  rb_define_module_function(module, "complex_dilog_e",  rb_gsl_sf_complex_dilog_e, 2);
}
