/*
  sf_exp.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "include/rb_gsl_sf.h"

static VALUE rb_gsl_sf_exp(VALUE obj, VALUE x)
{
  if (COMPLEX_P(x) || VECTOR_COMPLEX_P(x) || MATRIX_COMPLEX_P(x)) {
    return rb_gsl_math_complex_eval(gsl_complex_exp, x);
  } else {
    return rb_gsl_sf_eval1(gsl_sf_exp, x);
  }
}

static VALUE rb_gsl_sf_exp_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_exp_e, x);
}

static VALUE rb_gsl_sf_exp_e10_e(VALUE obj, VALUE x)
{
  gsl_sf_result_e10 *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result_e10, gsl_sf_result_e10, 0, free, rslt);
  /*status =*/ gsl_sf_exp_e10_e(NUM2DBL(x), rslt);
  return v;
}

static VALUE rb_gsl_sf_exp_mult(VALUE obj, VALUE x, VALUE y)
{
  Need_Float(x);   Need_Float(y);
  return rb_float_new(gsl_sf_exp_mult(NUM2DBL(x), NUM2DBL(y)));
}

static VALUE rb_gsl_sf_exp_mult_e(VALUE obj, VALUE x, VALUE y)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_exp_mult_e, x, y);
}

static VALUE rb_gsl_sf_exp_mult_e10_e(VALUE obj, VALUE x, VALUE y)
{
  gsl_sf_result_e10 *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(x);   Need_Float(y);
  v = Data_Make_Struct(cgsl_sf_result_e10, gsl_sf_result_e10, 0, free, rslt);
  /*status =*/ gsl_sf_exp_mult_e10_e(NUM2DBL(x), NUM2DBL(y), rslt);
  return v;
}

static VALUE rb_gsl_sf_expm1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_expm1, x);
}

static VALUE rb_gsl_sf_expm1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_expm1_e, x);
}

static VALUE rb_gsl_sf_exprel(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_exprel, x);
}

static VALUE rb_gsl_sf_exprel_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_exprel_e, x);
}

static VALUE rb_gsl_sf_exprel_2(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_exprel_2, x);
}

static VALUE rb_gsl_sf_exprel_2_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_exprel_2_e, x);
}

static VALUE rb_gsl_sf_exprel_n(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_exprel_n, n, x);
}

static VALUE rb_gsl_sf_exprel_n_e(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_exprel_n_e, n, x);
}

static VALUE rb_gsl_sf_exp_err_e(VALUE obj, VALUE x, VALUE dx)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_exp_err_e, x, dx);
}

static VALUE rb_gsl_sf_exp_err_e10_e(VALUE obj, VALUE x, VALUE dx)
{
  gsl_sf_result_e10 *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(x);   Need_Float(dx);
  v = Data_Make_Struct(cgsl_sf_result_e10, gsl_sf_result_e10, 0, free, rslt);
  /*status =*/ gsl_sf_exp_err_e10_e(NUM2DBL(x), NUM2DBL(dx), rslt);
  return v;
}

static VALUE rb_gsl_sf_exp_mult_err_e(VALUE obj, VALUE x, VALUE dx,
				      VALUE y, VALUE dy)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(x);   Need_Float(y);
  Need_Float(dx);   Need_Float(dy);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_exp_mult_err_e(NUM2DBL(x), NUM2DBL(dx), NUM2DBL(y), NUM2DBL(dy), rslt);
  return v;
}

static VALUE rb_gsl_sf_exp_mult_err_e10_e(VALUE obj, VALUE x, VALUE dx,
					  VALUE y, VALUE dy)
{
  gsl_sf_result_e10 *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(x);   Need_Float(y);
  Need_Float(dx);   Need_Float(dy);
  v = Data_Make_Struct(cgsl_sf_result_e10, gsl_sf_result_e10, 0, free, rslt);
  /*status =*/ gsl_sf_exp_mult_err_e10_e(NUM2DBL(x), NUM2DBL(dx), NUM2DBL(y), NUM2DBL(dy), rslt);
  return v;
}

void Init_gsl_sf_exp(VALUE module)
{
  rb_define_module_function(module, "exp",  rb_gsl_sf_exp, 1);
  rb_define_module_function(module, "exp_e",  rb_gsl_sf_exp_e, 1);
  rb_define_module_function(module, "exp_e10_e",  rb_gsl_sf_exp_e10_e, 1);
  rb_define_module_function(module, "exp_mult",  rb_gsl_sf_exp_mult, 2);
  rb_define_module_function(module, "exp_mult_e",  rb_gsl_sf_exp_mult_e, 2);
  rb_define_module_function(module, "exp_mult_e10_e",  rb_gsl_sf_exp_mult_e10_e, 2);
  rb_define_module_function(module, "expm1",  rb_gsl_sf_expm1, 1);
  rb_define_module_function(module, "expm1_e",  rb_gsl_sf_expm1_e, 1);
  rb_define_module_function(module, "exprel",  rb_gsl_sf_exprel, 1);
  rb_define_module_function(module, "exprel_e",  rb_gsl_sf_exprel_e, 1);
  rb_define_module_function(module, "exprel_2",  rb_gsl_sf_exprel_2, 1);
  rb_define_module_function(module, "exprel_2_e",  rb_gsl_sf_exprel_2_e, 1);
  rb_define_module_function(module, "exprel_n",  rb_gsl_sf_exprel_n, 2);
  rb_define_module_function(module, "exprel_n_e",  rb_gsl_sf_exprel_n_e, 2);
  rb_define_module_function(module, "exp_err_e",  rb_gsl_sf_exp_err_e, 2);
  rb_define_module_function(module, "exp_err_e10_e",  rb_gsl_sf_exp_err_e10_e, 2);
  rb_define_module_function(module, "exp_mult_err_e",  rb_gsl_sf_exp_mult_err_e, 4);
  rb_define_module_function(module, "exp_mult_err_e10_e",  rb_gsl_sf_exp_mult_err_e10_e, 4);
}
