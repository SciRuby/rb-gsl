/*
  sf_laguerre.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_sf.h"

static VALUE rb_gsl_sf_laguerre_X(int argc, VALUE *argv, VALUE obj,
				  double (*f)(double, double))
{
  switch (argc) {
  case 2:
    return rb_gsl_sf_eval_double_double(f, argv[0], argv[1]);
    break;
  case 1:
    return rb_gsl_sf_eval_double_double(f, INT2FIX(0), argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
}

static VALUE rb_gsl_sf_laguerre_1(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_sf_laguerre_X(argc, argv, obj, gsl_sf_laguerre_1);
}

static VALUE rb_gsl_sf_laguerre_1_e(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_laguerre_1_e, a, x);
}

static VALUE rb_gsl_sf_laguerre_2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_sf_laguerre_X(argc, argv, obj, gsl_sf_laguerre_2);
}

static VALUE rb_gsl_sf_laguerre_2_e(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_laguerre_2_e, a, x);
}

static VALUE rb_gsl_sf_laguerre_3(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_sf_laguerre_X(argc, argv, obj, gsl_sf_laguerre_3);
}

static VALUE rb_gsl_sf_laguerre_3_e(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_laguerre_3_e, a, x);
}

static VALUE rb_gsl_sf_laguerre_n(int argc, VALUE *argv, VALUE obj)
{
  switch (argc) {
  case 3:
    return rb_gsl_sf_eval_int_double_double(gsl_sf_laguerre_n, argv[0],
					    argv[1], argv[2]);
    break;
  case 2:
    return rb_gsl_sf_eval_int_double_double(gsl_sf_laguerre_n, argv[0],
					    INT2FIX(0), argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
}

static VALUE rb_gsl_sf_laguerre_n_e(VALUE obj, VALUE n, VALUE a, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  CHECK_FIXNUM(n);
  Need_Float(a); Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_laguerre_n_e(FIX2INT(n), NUM2DBL(a), NUM2DBL(x), rslt);
  return v;
}

void Init_gsl_sf_laguerre(VALUE module)
{
  VALUE mgsl_sf_laguerre;

  rb_define_module_function(module, "laguerre_1",  rb_gsl_sf_laguerre_1, -1);
  rb_define_module_function(module, "laguerre_1_e",  rb_gsl_sf_laguerre_1_e, 2);
  rb_define_module_function(module, "laguerre_2",  rb_gsl_sf_laguerre_2, -1);
  rb_define_module_function(module, "laguerre_2_e",  rb_gsl_sf_laguerre_2_e, 2);
  rb_define_module_function(module, "laguerre_3",  rb_gsl_sf_laguerre_3, -1);
  rb_define_module_function(module, "laguerre_3_e",  rb_gsl_sf_laguerre_3_e, 2);
  rb_define_module_function(module, "laguerre_n",  rb_gsl_sf_laguerre_n, -1);
  rb_define_module_function(module, "laguerre_n_e",  rb_gsl_sf_laguerre_n_e, 3);

  mgsl_sf_laguerre = rb_define_module_under(module, "Laguerre");
  rb_define_module_function(mgsl_sf_laguerre, "one",  rb_gsl_sf_laguerre_1, -1);
  rb_define_module_function(mgsl_sf_laguerre, "one_e",  rb_gsl_sf_laguerre_1_e, 2);
  rb_define_module_function(mgsl_sf_laguerre, "two",  rb_gsl_sf_laguerre_2, -1);
  rb_define_module_function(mgsl_sf_laguerre, "two_e",  rb_gsl_sf_laguerre_2_e, 2);
  rb_define_module_function(mgsl_sf_laguerre, "three_3",  rb_gsl_sf_laguerre_3, -1);
  rb_define_module_function(mgsl_sf_laguerre, "three_e",  rb_gsl_sf_laguerre_3_e, 2);
  rb_define_module_function(mgsl_sf_laguerre, "n",  rb_gsl_sf_laguerre_n, -1);
  rb_define_module_function(mgsl_sf_laguerre, "n_e",  rb_gsl_sf_laguerre_n_e, 3);

}
