/*
  sf_power.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_sf.h"

VALUE rb_gsl_complex_pow(int argc, VALUE *argv, VALUE obj);
static VALUE rb_gsl_sf_pow_int(VALUE obj, VALUE x, VALUE n)
{
  VALUE argv[2];
  if (COMPLEX_P(x) || VECTOR_COMPLEX_P(x) || MATRIX_COMPLEX_P(x)) {
    argv[0] = x;
    argv[1] = n;
    return rb_gsl_complex_pow(2, argv, obj);
  }
  return rb_gsl_sf_eval_double_int(gsl_sf_pow_int, x, n);
}

static VALUE rb_gsl_sf_pow_int_e(VALUE obj, VALUE x, VALUE n)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(x);
  CHECK_FIXNUM(n);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_pow_int_e(NUM2DBL(x), FIX2INT(n), rslt);
  return v;
}

void Init_gsl_sf_power(VALUE module)
{
  VALUE mgsl_sf_pow;
  rb_define_module_function(module, "pow_int",  rb_gsl_sf_pow_int, 2);
  rb_define_module_function(module, "pow_int_e",  rb_gsl_sf_pow_int_e, 2);
  mgsl_sf_pow = rb_define_module_under(module, "Pow");
  rb_define_module_function(mgsl_sf_pow, "int",  rb_gsl_sf_pow_int, 2);
  rb_define_module_function(mgsl_sf_pow, "int_e",  rb_gsl_sf_pow_int_e, 2);
}
