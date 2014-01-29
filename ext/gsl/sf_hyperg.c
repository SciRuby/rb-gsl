/*
  sf_hyperg.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_sf.h"

/* Checked with Mathematica. */
static VALUE rb_gsl_sf_hyperg_0F1(VALUE obj, VALUE c, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_hyperg_0F1, c, x);
}

static VALUE rb_gsl_sf_hyperg_0F1_e(VALUE obj, VALUE c, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_hyperg_0F1_e, c, x);
}

static VALUE rb_gsl_sf_hyperg_1F1_int(VALUE obj, VALUE m, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_int_double(gsl_sf_hyperg_1F1_int, m, n, x);
}

static VALUE rb_gsl_sf_hyperg_1F1_int_e(VALUE obj, VALUE m, VALUE n, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  CHECK_FIXNUM(m); CHECK_FIXNUM(n); Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_hyperg_1F1_int_e(FIX2INT(m), FIX2INT(n), NUM2DBL(x), rslt);
  return v;
}

static VALUE rb_gsl_sf_hyperg_1F1(VALUE obj, VALUE a, VALUE b, VALUE x)
{
  return rb_gsl_sf_eval_double3(gsl_sf_hyperg_1F1, a, b, x);
}

static VALUE rb_gsl_sf_hyperg_1F1_e(VALUE obj, VALUE a, VALUE b, VALUE x)
{
  return rb_gsl_sf_eval_e_double3(gsl_sf_hyperg_1F1_e, a, b, x);
}

static VALUE rb_gsl_sf_hyperg_U_int(VALUE obj, VALUE m, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_int_double(gsl_sf_hyperg_U_int, m, n, x);
}

static VALUE rb_gsl_sf_hyperg_U_int_e(VALUE obj, VALUE m, VALUE n, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  CHECK_FIXNUM(m); CHECK_FIXNUM(n); Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_hyperg_U_int_e(FIX2INT(m), FIX2INT(n), NUM2DBL(x), rslt);
  return v;
}

static VALUE rb_gsl_sf_hyperg_U_int_e10_e(VALUE obj, VALUE m, VALUE n, VALUE x)
{
  gsl_sf_result_e10 *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  CHECK_FIXNUM(m); CHECK_FIXNUM(n); Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result_e10, gsl_sf_result_e10, 0, free, rslt);
  /*status =*/ gsl_sf_hyperg_U_int_e10_e(FIX2INT(m), FIX2INT(n), NUM2DBL(x), rslt);
  return v;
}

static VALUE rb_gsl_sf_hyperg_U(VALUE obj, VALUE a, VALUE b, VALUE x)
{
  return rb_gsl_sf_eval_double3(gsl_sf_hyperg_U, a, b, x);
}

static VALUE rb_gsl_sf_hyperg_U_e(VALUE obj, VALUE a, VALUE b, VALUE x)
{
  return rb_gsl_sf_eval_e_double3(gsl_sf_hyperg_U_e, a, b, x);
}

static VALUE rb_gsl_sf_hyperg_U_e10_e(VALUE obj, VALUE a, VALUE b, VALUE x)
{
  gsl_sf_result_e10 *rslt = NULL;
  VALUE v;
  int status;
  Need_Float(a); Need_Float(b); Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result_e10, gsl_sf_result_e10, 0, free, rslt);
  status = gsl_sf_hyperg_U_e10_e(NUM2DBL(a), NUM2DBL(b), NUM2DBL(x), rslt);
  return rb_ary_new3(2, v, INT2FIX(status));
}

static VALUE rb_gsl_sf_hyperg_2F1(VALUE obj, VALUE a, VALUE b, VALUE c, VALUE x)
{
  return rb_gsl_sf_eval_double4(gsl_sf_hyperg_2F1, a, b, c, x);
}

static VALUE rb_gsl_sf_hyperg_2F1_e(VALUE obj, VALUE a, VALUE b, VALUE c, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  int status;
  Need_Float(a); Need_Float(b); Need_Float(c); Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  status = gsl_sf_hyperg_2F1_e(NUM2DBL(a), NUM2DBL(b), NUM2DBL(c), NUM2DBL(x), rslt);
  return rb_ary_new3(2, v, INT2FIX(status));
}

static VALUE rb_gsl_sf_hyperg_2F1_conj(VALUE obj, VALUE aR, VALUE aI, VALUE c, VALUE x)
{
  return rb_gsl_sf_eval_double4(gsl_sf_hyperg_2F1_conj, aR, aI, c, x);
}

static VALUE rb_gsl_sf_hyperg_2F1_conj_e(VALUE obj, VALUE aR, VALUE aI, VALUE c, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(aR); Need_Float(aI);  Need_Float(c); Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_hyperg_2F1_conj_e(NUM2DBL(aR), NUM2DBL(aI), NUM2DBL(c), NUM2DBL(x), rslt);
  return v;
}

static VALUE rb_gsl_sf_hyperg_2F1_renorm(VALUE obj, VALUE a, VALUE b, VALUE c, VALUE x)
{
  return rb_gsl_sf_eval_double4(gsl_sf_hyperg_2F1_renorm, a, b, c, x);
}

static VALUE rb_gsl_sf_hyperg_2F1_renorm_e(VALUE obj, VALUE a, VALUE b, VALUE c, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(a); Need_Float(b);  Need_Float(c); Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_hyperg_2F1_renorm_e(NUM2DBL(a), NUM2DBL(b), NUM2DBL(c), NUM2DBL(x), rslt);
  return v;
}

static VALUE rb_gsl_sf_hyperg_2F1_conj_renorm(VALUE obj, VALUE aR, VALUE aI, VALUE c, VALUE x)
{
  return rb_gsl_sf_eval_double4(gsl_sf_hyperg_2F1_conj_renorm, aR, aI, c, x);
}

static VALUE rb_gsl_sf_hyperg_2F1_conj_renorm_e(VALUE obj, VALUE aR, VALUE aI, VALUE c, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(aR); Need_Float(aI);  Need_Float(c); Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_hyperg_2F1_conj_renorm_e(NUM2DBL(aR), NUM2DBL(aI), NUM2DBL(c), NUM2DBL(x), rslt);
  return v;
}

static VALUE rb_gsl_sf_hyperg_2F0(VALUE obj, VALUE a, VALUE b,  VALUE x)
{
  return rb_gsl_sf_eval_double3(gsl_sf_hyperg_2F0, a, b, x);
}

static VALUE rb_gsl_sf_hyperg_2F0_e(VALUE obj, VALUE a, VALUE b, VALUE x)
{
  return rb_gsl_sf_eval_e_double3(gsl_sf_hyperg_2F0_e, a, b, x);
}

void Init_gsl_sf_hyperg(VALUE module)
{
  rb_define_module_function(module, "hyperg_0F1",  rb_gsl_sf_hyperg_0F1, 2);
  rb_define_module_function(module, "hyperg_0F1_e",  rb_gsl_sf_hyperg_0F1_e, 2);
  rb_define_module_function(module, "hyperg_1F1_int",  rb_gsl_sf_hyperg_1F1_int, 3);
  rb_define_module_function(module, "hyperg_1F1_int_e",  rb_gsl_sf_hyperg_1F1_int_e, 3);
  rb_define_module_function(module, "hyperg_1F1",  rb_gsl_sf_hyperg_1F1, 3);
  rb_define_module_function(module, "hyperg_1F1_e",  rb_gsl_sf_hyperg_1F1_e, 3);
  rb_define_module_function(module, "hyperg_U_int",  rb_gsl_sf_hyperg_U_int, 3);
  rb_define_module_function(module, "hyperg_U_int_e",  rb_gsl_sf_hyperg_U_int_e, 3);
  rb_define_module_function(module, "hyperg_U_int_e10_e",  rb_gsl_sf_hyperg_U_int_e10_e, 3);
  rb_define_module_function(module, "hyperg_U",  rb_gsl_sf_hyperg_U, 3);
  rb_define_module_function(module, "hyperg_U_e",  rb_gsl_sf_hyperg_U_e, 3);
  rb_define_module_function(module, "hyperg_U_e10_e",  rb_gsl_sf_hyperg_U_e10_e, 3);
  rb_define_module_function(module, "hyperg_2F1",  rb_gsl_sf_hyperg_2F1, 4);
  rb_define_module_function(module, "hyperg_2F1_e",  rb_gsl_sf_hyperg_2F1_e, 4);
  rb_define_module_function(module, "hyperg_2F1_conj",  rb_gsl_sf_hyperg_2F1_conj, 4);
  rb_define_module_function(module, "hyperg_2F1_conj_e",  rb_gsl_sf_hyperg_2F1_conj_e, 4);
  rb_define_module_function(module, "hyperg_2F1_renorm",  rb_gsl_sf_hyperg_2F1_renorm, 4);
  rb_define_module_function(module, "hyperg_2F1_renorm_e",  rb_gsl_sf_hyperg_2F1_renorm_e, 4);
  rb_define_module_function(module, "hyperg_2F1_conj_renorm",  rb_gsl_sf_hyperg_2F1_conj_renorm, 4);
  rb_define_module_function(module, "hyperg_2F1_conj_renorm_e",  rb_gsl_sf_hyperg_2F1_conj_renorm_e, 4);
  rb_define_module_function(module, "hyperg_2F0",  rb_gsl_sf_hyperg_2F0, 3);
  rb_define_module_function(module, "hyperg_2F0_e",  rb_gsl_sf_hyperg_2F0_e, 3);
}
