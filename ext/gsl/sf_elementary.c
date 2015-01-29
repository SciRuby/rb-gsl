/*
  sf_elementary.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "include/rb_gsl_sf.h"

static VALUE rb_gsl_sf_multiply_e(VALUE obj, VALUE x, VALUE y)
{
  gsl_sf_result *r;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(x); Need_Float(y);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, r);
  /*status =*/ gsl_sf_multiply_e(NUM2DBL(x), NUM2DBL(y), r);
  return v;
}

static VALUE rb_gsl_sf_multiply_err_e(VALUE obj, VALUE x, VALUE dx,
              VALUE y, VALUE dy)
{
  gsl_sf_result *r;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(x); Need_Float(y);
  Need_Float(dx); Need_Float(dy);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, r);
  /*status =*/ gsl_sf_multiply_err_e(NUM2DBL(x), NUM2DBL(dx),
         NUM2DBL(y), NUM2DBL(dy), r);
  return v;
}

void Init_gsl_sf_elementary(VALUE module)
{
  rb_define_module_function(module, "multiply_e",  rb_gsl_sf_multiply_e, 2);
  rb_define_module_function(module, "multiply_err_e",  rb_gsl_sf_multiply_err_e, 4);
}
