/*
  sf_trigonometric.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_sf.h"

static VALUE rb_gsl_sf_sin(VALUE obj, VALUE x)
{
  if (COMPLEX_P(x) || VECTOR_COMPLEX_P(x) || MATRIX_COMPLEX_P(x)) 
    return rb_gsl_math_complex_eval(gsl_complex_sin, x);
  return rb_gsl_sf_eval1(gsl_sf_sin, x);
}

static VALUE rb_gsl_sf_sin_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_sin_e, x);
}

static VALUE rb_gsl_sf_cos(VALUE obj, VALUE x)
{
  if (COMPLEX_P(x) || VECTOR_COMPLEX_P(x) || MATRIX_COMPLEX_P(x)) 
    return rb_gsl_math_complex_eval(gsl_complex_cos, x);
  return rb_gsl_sf_eval1(gsl_sf_cos, x);
}

static VALUE rb_gsl_sf_cos_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_cos_e, x);
}

static VALUE rb_gsl_sf_hypot(VALUE obj, VALUE x, VALUE y)
{
  Need_Float(x); Need_Float(y);
  return rb_float_new(gsl_sf_hypot(NUM2DBL(x), NUM2DBL(y)));
}

static VALUE rb_gsl_sf_hypot_e(VALUE obj, VALUE x, VALUE y)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_hypot_e, x, y);
}

static VALUE rb_gsl_sf_sinc(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_sinc, x);
}

static VALUE rb_gsl_sf_sinc_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_sinc_e, x);
}

static VALUE rb_gsl_sf_complex_XXX_e(int argc, VALUE *argv, VALUE obj,
				     int (*f)(double, double, gsl_sf_result*, gsl_sf_result*))
{
  gsl_sf_result *r1, *r2;
  gsl_complex *z;
  double re, im;
  VALUE v1, v2;
  // local variable "status" declared and set, but never used
  //int status;
  switch (argc) {
  case 1:
    CHECK_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_complex, z);
    re = GSL_REAL(*z);
    im = GSL_IMAG(*z);
    break;
  case 2:
    Need_Float(argv[0]); Need_Float(argv[1]);
    re = NUM2DBL(argv[0]);
    im = NUM2DBL(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  v1 = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, r1);
  v2 = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, r2);
  /*status =*/ (*f)(re, im, r1, r2);
  return rb_ary_new3(2, v1, v2);
}

static VALUE rb_gsl_sf_complex_sin_e(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_sf_complex_XXX_e(argc, argv, obj, gsl_sf_complex_sin_e);
}

static VALUE rb_gsl_sf_complex_cos_e(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_sf_complex_XXX_e(argc, argv, obj, gsl_sf_complex_cos_e);
}

static VALUE rb_gsl_sf_complex_logsin_e(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_sf_complex_XXX_e(argc, argv, obj, gsl_sf_complex_logsin_e);
}

static VALUE rb_gsl_sf_lnsinh(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_lnsinh, x);
}

static VALUE rb_gsl_sf_lnsinh_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_lnsinh_e, x);
}

static VALUE rb_gsl_sf_lncosh(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_lncosh, x);
}

static VALUE rb_gsl_sf_lncosh_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_lncosh_e, x);
}

static VALUE rb_gsl_sf_polar_to_rect(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_sf_complex_XXX_e(argc, argv, obj, gsl_sf_polar_to_rect);
}

static VALUE rb_gsl_sf_rect_to_polar(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_sf_complex_XXX_e(argc, argv, obj, gsl_sf_rect_to_polar);
}

static VALUE rb_gsl_sf_angle_restrict_symm(VALUE obj, VALUE theta)
{
  return rb_gsl_sf_eval1(gsl_sf_angle_restrict_symm, theta);
}

static VALUE rb_gsl_sf_angle_restrict_pos(VALUE obj, VALUE theta)
{
  return rb_gsl_sf_eval1(gsl_sf_angle_restrict_pos, theta);
}

/*
static VALUE rb_gsl_sf_sin_err(VALUE obj, VALUE x, VALUE dx)
{
  return rb_float_new(gsl_sf_sin_err(NUM2DBL(x), NUM2DBL(dx)));
}

static VALUE rb_gsl_sf_cos_err(VALUE obj, VALUE x, VALUE dx)
{
  return rb_float_new(gsl_sf_cos_err(NUM2DBL(x), NUM2DBL(dx)));
}
*/

static VALUE rb_gsl_sf_sin_err_e(VALUE obj, VALUE x, VALUE dx)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(x);
  Need_Float(dx);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_sin_err_e(NUM2DBL(x), NUM2DBL(dx), rslt);
  return v;
}

static VALUE rb_gsl_sf_cos_err_e(VALUE obj, VALUE x, VALUE dx)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(x);
  Need_Float(dx);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_cos_err_e(NUM2DBL(x), NUM2DBL(dx), rslt);
  return v;
}

void Init_gsl_sf_trigonometric(VALUE module)
{
  rb_define_module_function(module, "sin",  rb_gsl_sf_sin, 1);
  rb_define_module_function(module, "sin_e",  rb_gsl_sf_sin_e, 1);
  rb_define_module_function(module, "cos",  rb_gsl_sf_cos, 1);
  rb_define_module_function(module, "cos_e",  rb_gsl_sf_cos_e, 1);
  rb_define_module_function(module, "hypot",  rb_gsl_sf_hypot, 2);
  rb_define_module_function(module, "hypot_e",  rb_gsl_sf_hypot_e, 2);
  rb_define_module_function(module, "sinc",  rb_gsl_sf_sinc, 1);
  rb_define_module_function(module, "sinc_e",  rb_gsl_sf_sinc_e, 1);
  rb_define_module_function(module, "complex_sin_e",  rb_gsl_sf_complex_sin_e, -1);
  rb_define_module_function(module, "complex_cos_e",  rb_gsl_sf_complex_cos_e, -1);
  rb_define_module_function(module, "complex_logsin_e",  rb_gsl_sf_complex_logsin_e, -1);
  rb_define_module_function(module, "lnsinh",  rb_gsl_sf_lnsinh, 1);
  rb_define_module_function(module, "lnsinh_e",  rb_gsl_sf_lnsinh_e, 1);
  rb_define_module_function(module, "lncosh",  rb_gsl_sf_lncosh, 1);
  rb_define_module_function(module, "lncosh_e",  rb_gsl_sf_lncosh_e, 1);
  rb_define_module_function(module, "polar_to_rect",  rb_gsl_sf_polar_to_rect, -1);
  rb_define_module_function(module, "rect_to_polar",  rb_gsl_sf_rect_to_polar, -1);
  rb_define_module_function(module, "angle_restrict_symm",  rb_gsl_sf_angle_restrict_symm, 1);
  rb_define_module_function(module, "angle_restrict_pos",  rb_gsl_sf_angle_restrict_pos, 1);

  /*  rb_define_module_function(module, "sin_err",  rb_gsl_sf_sin_err, 2);
      rb_define_module_function(module, "cos_err",  rb_gsl_sf_cos_err, 2);*/
  rb_define_module_function(module, "sin_err_e",  rb_gsl_sf_sin_err_e, 2);
  rb_define_module_function(module, "cos_err_e",  rb_gsl_sf_cos_err_e, 2);

}
