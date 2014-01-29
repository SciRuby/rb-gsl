/*
  sf_gegenbauer.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_sf.h"

static VALUE rb_gsl_sf_gegenpoly_1(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_gegenpoly_1, lambda, x);
}

static VALUE rb_gsl_sf_gegenpoly_1_e(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_gegenpoly_1_e, lambda, x);
}

static VALUE rb_gsl_sf_gegenpoly_2(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_gegenpoly_2, lambda, x);
}

static VALUE rb_gsl_sf_gegenpoly_2_e(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_gegenpoly_2_e, lambda, x);
}

static VALUE rb_gsl_sf_gegenpoly_3(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_gegenpoly_3, lambda, x);
}

static VALUE rb_gsl_sf_gegenpoly_3_e(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_gegenpoly_3_e, lambda, x);
}

static VALUE rb_gsl_sf_gegenpoly_n(VALUE obj, VALUE n, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_int_double_double(gsl_sf_gegenpoly_n, n, lambda, x);
}

static VALUE rb_gsl_sf_gegenpoly_n_e(VALUE obj, VALUE n, VALUE lambda, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  CHECK_FIXNUM(n);
  Need_Float(lambda); Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_gegenpoly_n_e(FIX2INT(n), NUM2DBL(lambda), NUM2DBL(x), rslt);

  return v;
}

static VALUE rb_gsl_sf_gegenpoly_array(VALUE obj, VALUE nmax, VALUE lambda, VALUE x)
{
  gsl_vector *v = NULL;
  CHECK_FIXNUM(nmax);
  Need_Float(lambda); Need_Float(x);
  v = gsl_vector_alloc(nmax);
  gsl_sf_gegenpoly_array(FIX2INT(nmax), NUM2DBL(lambda), NUM2DBL(x), v->data);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

void Init_gsl_sf_gegenbauer(VALUE module)
{
  VALUE mgsl_sf_gegenpoly;
  rb_define_module_function(module, "gegenpoly_1",  rb_gsl_sf_gegenpoly_1, 2);
  rb_define_module_function(module, "gegenpoly_1_e",  rb_gsl_sf_gegenpoly_1_e, 2);
  rb_define_module_function(module, "gegenpoly_2",  rb_gsl_sf_gegenpoly_2, 2);
  rb_define_module_function(module, "gegenpoly_2_e",  rb_gsl_sf_gegenpoly_2_e, 2);
  rb_define_module_function(module, "gegenpoly_3",  rb_gsl_sf_gegenpoly_3, 2);
  rb_define_module_function(module, "gegenpoly_3_e",  rb_gsl_sf_gegenpoly_3_e, 2);
  rb_define_module_function(module, "gegenpoly_n",  rb_gsl_sf_gegenpoly_n, 3);
  rb_define_module_function(module, "gegenpoly_n_e",  rb_gsl_sf_gegenpoly_n_e, 3);
  rb_define_module_function(module, "gegenpoly_array",  rb_gsl_sf_gegenpoly_array, 3);

  mgsl_sf_gegenpoly = rb_define_module_under(module, "Gegenpoly");
  rb_define_module_function(mgsl_sf_gegenpoly, "one",  rb_gsl_sf_gegenpoly_1, 2);
  rb_define_module_function(mgsl_sf_gegenpoly, "one_e",  rb_gsl_sf_gegenpoly_1_e, 2);
  rb_define_module_function(mgsl_sf_gegenpoly, "two",  rb_gsl_sf_gegenpoly_2, 2);
  rb_define_module_function(mgsl_sf_gegenpoly, "two_e",  rb_gsl_sf_gegenpoly_2_e, 2);
  rb_define_module_function(mgsl_sf_gegenpoly, "three",  rb_gsl_sf_gegenpoly_3, 2);
  rb_define_module_function(mgsl_sf_gegenpoly, "three_e",  rb_gsl_sf_gegenpoly_3_e, 2);
  rb_define_module_function(mgsl_sf_gegenpoly, "n",  rb_gsl_sf_gegenpoly_n, 3);
  rb_define_module_function(mgsl_sf_gegenpoly, "n_e",  rb_gsl_sf_gegenpoly_n_e, 3);
  rb_define_module_function(mgsl_sf_gegenpoly, "array",  rb_gsl_sf_gegenpoly_array, 3);

}
