/*
  sf_airy.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "include/rb_gsl_sf.h"
#include "include/rb_gsl_array.h"

/* m: precision,  d(double), s(single), a(apporox) */
static VALUE rb_gsl_sf_airy_Ai(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) return eval_sf(gsl_sf_airy_Ai, argv[0]);
  else return rb_gsl_sf_eval_double_m(gsl_sf_airy_Ai, argv[0], argv[1]);
}

static VALUE rb_gsl_sf_airy_Ai_e(VALUE obj, VALUE x, VALUE m)
{
  return rb_gsl_sf_eval_e_m(gsl_sf_airy_Ai_e, x, m);
}

static VALUE rb_gsl_sf_airy_Bi(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) return eval_sf(gsl_sf_airy_Bi, argv[0]);
  else return rb_gsl_sf_eval_double_m(gsl_sf_airy_Bi, argv[0], argv[1]);
}

static VALUE rb_gsl_sf_airy_Bi_e(VALUE obj, VALUE x, VALUE m)
{
  return rb_gsl_sf_eval_e_m(gsl_sf_airy_Bi_e, x, m);
}

static VALUE rb_gsl_sf_airy_Ai_scaled(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) return eval_sf(gsl_sf_airy_Ai_scaled, argv[0]);
  else return rb_gsl_sf_eval_double_m(gsl_sf_airy_Ai_scaled, argv[0], argv[1]);
}

static VALUE rb_gsl_sf_airy_Ai_scaled_e(VALUE obj, VALUE x, VALUE m)
{
  return rb_gsl_sf_eval_e_m(gsl_sf_airy_Ai_scaled_e, x, m);
}

static VALUE rb_gsl_sf_airy_Bi_scaled(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) return eval_sf(gsl_sf_airy_Bi_scaled, argv[0]);
  else return rb_gsl_sf_eval_double_m(gsl_sf_airy_Bi_scaled, argv[0], argv[1]);
}

static VALUE rb_gsl_sf_airy_Bi_scaled_e(VALUE obj, VALUE x, VALUE m)
{
  return rb_gsl_sf_eval_e_m(gsl_sf_airy_Bi_scaled_e, x, m);
}

static VALUE rb_gsl_sf_airy_Ai_deriv(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) return eval_sf(gsl_sf_airy_Ai_deriv, argv[0]);
  else return rb_gsl_sf_eval_double_m(gsl_sf_airy_Ai_deriv, argv[0], argv[1]);
}

static VALUE rb_gsl_sf_airy_Ai_deriv_e(VALUE obj, VALUE x, VALUE m)
{
  return rb_gsl_sf_eval_e_m(gsl_sf_airy_Ai_deriv_e, x, m);
}

static VALUE rb_gsl_sf_airy_Bi_deriv(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) return eval_sf(gsl_sf_airy_Bi_deriv, argv[0]);
  else return rb_gsl_sf_eval_double_m(gsl_sf_airy_Bi_deriv, argv[0], argv[1]);
}

static VALUE rb_gsl_sf_airy_Bi_deriv_e(VALUE obj, VALUE x, VALUE m)
{
  return rb_gsl_sf_eval_e_m(gsl_sf_airy_Bi_deriv_e, x, m);
}

static VALUE rb_gsl_sf_airy_Ai_deriv_scaled(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) return eval_sf(gsl_sf_airy_Ai_deriv_scaled, argv[0]);
  else return rb_gsl_sf_eval_double_m(gsl_sf_airy_Ai_deriv_scaled, argv[0], argv[1]);
}

static VALUE rb_gsl_sf_airy_Ai_deriv_scaled_e(VALUE obj, VALUE x, VALUE m)
{
  return rb_gsl_sf_eval_e_m(gsl_sf_airy_Ai_deriv_scaled_e, x, m);
}

static VALUE rb_gsl_sf_airy_Bi_deriv_scaled(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) return eval_sf(gsl_sf_airy_Bi_deriv_scaled, argv[0]);
  else return rb_gsl_sf_eval_double_m(gsl_sf_airy_Bi_deriv_scaled, argv[0], argv[1]);
}

static VALUE rb_gsl_sf_airy_Bi_deriv_scaled_e(VALUE obj, VALUE x, VALUE m)
{
  return rb_gsl_sf_eval_e_m(gsl_sf_airy_Bi_deriv_scaled_e, x, m);
}

static VALUE rb_gsl_sf_airy_zero_Ai(VALUE obj, VALUE s)
{
  return rb_gsl_sf_eval1_uint(gsl_sf_airy_zero_Ai, s);
}

static VALUE rb_gsl_sf_airy_zero_Ai_e(VALUE obj, VALUE s)
{
  return rb_gsl_sf_eval_e_uint(gsl_sf_airy_zero_Ai_e, s);
}

static VALUE rb_gsl_sf_airy_zero_Bi(VALUE obj, VALUE s)
{
  return rb_gsl_sf_eval1_uint(gsl_sf_airy_zero_Bi, s);
}

static VALUE rb_gsl_sf_airy_zero_Bi_e(VALUE obj, VALUE s)
{
  return rb_gsl_sf_eval_e_uint(gsl_sf_airy_zero_Bi_e, s);
}

static VALUE rb_gsl_sf_airy_zero_Ai_deriv(VALUE obj, VALUE s)
{
  return rb_gsl_sf_eval1_uint(gsl_sf_airy_zero_Ai_deriv, s);
}

static VALUE rb_gsl_sf_airy_zero_Ai_deriv_e(VALUE obj, VALUE s)
{
  return rb_gsl_sf_eval_e_uint(gsl_sf_airy_zero_Ai_deriv_e, s);
}

static VALUE rb_gsl_sf_airy_zero_Bi_deriv(VALUE obj, VALUE s)
{
  return rb_gsl_sf_eval1_uint(gsl_sf_airy_zero_Bi_deriv, s);
}

static VALUE rb_gsl_sf_airy_zero_Bi_deriv_e(VALUE obj, VALUE s)
{
  return rb_gsl_sf_eval_e_uint(gsl_sf_airy_zero_Bi_deriv_e, s);
}

void Init_gsl_sf_airy(VALUE module)
{
  VALUE mgsl_sf_airy;

  rb_define_module_function(module, "airy_Ai",  rb_gsl_sf_airy_Ai, -1);
  rb_define_module_function(module, "airy_Ai_e",  rb_gsl_sf_airy_Ai_e, 2);
  rb_define_module_function(module, "airy_Bi",  rb_gsl_sf_airy_Bi, -1);
  rb_define_module_function(module, "airy_Bi_e",  rb_gsl_sf_airy_Bi_e, 2);
  rb_define_module_function(module, "airy_Ai_scaled",  rb_gsl_sf_airy_Ai_scaled, -1);
  rb_define_module_function(module, "airy_Ai_scaled_e",  rb_gsl_sf_airy_Ai_scaled_e, 2);
  rb_define_module_function(module, "airy_Bi_scaled",  rb_gsl_sf_airy_Bi_scaled, -1);
  rb_define_module_function(module, "airy_Bi_scaled_e",  rb_gsl_sf_airy_Bi_scaled_e, 2);
  rb_define_module_function(module, "airy_Ai_deriv",  rb_gsl_sf_airy_Ai_deriv, -1);
  rb_define_module_function(module, "airy_Ai_deriv_e",  rb_gsl_sf_airy_Ai_deriv_e, 2);
  rb_define_module_function(module, "airy_Bi_deriv",  rb_gsl_sf_airy_Bi_deriv, -1);
  rb_define_module_function(module, "airy_Bi_deriv_e",  rb_gsl_sf_airy_Bi_deriv_e, 2);
  rb_define_module_function(module, "airy_Ai_deriv_scaled",  rb_gsl_sf_airy_Ai_deriv_scaled, -1);
  rb_define_module_function(module, "airy_Ai_deriv_scaled_e",  rb_gsl_sf_airy_Ai_deriv_scaled_e, 2);
  rb_define_module_function(module, "airy_Bi_deriv_scaled",  rb_gsl_sf_airy_Bi_deriv_scaled, -1);
  rb_define_module_function(module, "airy_Bi_deriv_scaled_e",  rb_gsl_sf_airy_Bi_deriv_scaled_e, 2);
  rb_define_module_function(module, "airy_zero_Ai",  rb_gsl_sf_airy_zero_Ai, 1);
  rb_define_module_function(module, "airy_zero_Ai_e",  rb_gsl_sf_airy_zero_Ai_e, 1);
  rb_define_module_function(module, "airy_zero_Bi",  rb_gsl_sf_airy_zero_Bi, 1);
  rb_define_module_function(module, "airy_zero_Bi_e",  rb_gsl_sf_airy_zero_Bi_e, 1);
  rb_define_module_function(module, "airy_zero_Ai_deriv",  rb_gsl_sf_airy_zero_Ai_deriv, 1);
  rb_define_module_function(module, "airy_zero_Ai_deriv_e",  rb_gsl_sf_airy_zero_Ai_deriv_e, 1);
  rb_define_module_function(module, "airy_zero_Bi_deriv",  rb_gsl_sf_airy_zero_Bi_deriv, 1);
  rb_define_module_function(module, "airy_zero_Bi_deriv_e",  rb_gsl_sf_airy_zero_Bi_deriv_e, 1);

  mgsl_sf_airy = rb_define_module_under(module, "Airy");

  rb_define_module_function(mgsl_sf_airy, "Ai",  rb_gsl_sf_airy_Ai, -1);
  rb_define_module_function(mgsl_sf_airy, "Ai_e",  rb_gsl_sf_airy_Ai_e, 2);
  rb_define_module_function(mgsl_sf_airy, "Bi",  rb_gsl_sf_airy_Bi, -1);
  rb_define_module_function(mgsl_sf_airy, "Bi_e",  rb_gsl_sf_airy_Bi_e, 2);
  rb_define_module_function(mgsl_sf_airy, "Ai_scaled",  rb_gsl_sf_airy_Ai_scaled, -1);
  rb_define_module_function(mgsl_sf_airy, "Ai_scaled_e",  rb_gsl_sf_airy_Ai_scaled_e, 2);
  rb_define_module_function(mgsl_sf_airy, "Bi_scaled",  rb_gsl_sf_airy_Bi_scaled, -1);
  rb_define_module_function(mgsl_sf_airy, "Bi_scaled_e",  rb_gsl_sf_airy_Bi_scaled_e, 2);
  rb_define_module_function(mgsl_sf_airy, "Ai_deriv",  rb_gsl_sf_airy_Ai_deriv, -1);
  rb_define_module_function(mgsl_sf_airy, "Ai_deriv_e",  rb_gsl_sf_airy_Ai_deriv_e, 2);
  rb_define_module_function(mgsl_sf_airy, "Bi_deriv",  rb_gsl_sf_airy_Bi_deriv, -1);
  rb_define_module_function(mgsl_sf_airy, "Bi_deriv_e",  rb_gsl_sf_airy_Bi_deriv_e, 2);
  rb_define_module_function(mgsl_sf_airy, "Ai_deriv_scaled",  rb_gsl_sf_airy_Ai_deriv_scaled, -1);
  rb_define_module_function(mgsl_sf_airy, "Ai_deriv_scaled_e",  rb_gsl_sf_airy_Ai_deriv_scaled_e, 2);
  rb_define_module_function(mgsl_sf_airy, "Bi_deriv_scaled",  rb_gsl_sf_airy_Bi_deriv_scaled, -1);
  rb_define_module_function(mgsl_sf_airy, "Bi_deriv_scaled_e",  rb_gsl_sf_airy_Bi_deriv_scaled_e, 2);
  rb_define_module_function(mgsl_sf_airy, "zero_Ai",  rb_gsl_sf_airy_zero_Ai, 1);
  rb_define_module_function(mgsl_sf_airy, "zero_Ai_e",  rb_gsl_sf_airy_zero_Ai_e, 1);
  rb_define_module_function(mgsl_sf_airy, "zero_Bi",  rb_gsl_sf_airy_zero_Bi, 1);
  rb_define_module_function(mgsl_sf_airy, "zero_Bi_e",  rb_gsl_sf_airy_zero_Bi_e, 1);
  rb_define_module_function(mgsl_sf_airy, "zero_Ai_deriv",  rb_gsl_sf_airy_zero_Ai_deriv, 1);
  rb_define_module_function(mgsl_sf_airy, "zero_Ai_deriv_e",  rb_gsl_sf_airy_zero_Ai_deriv_e, 1);
  rb_define_module_function(mgsl_sf_airy, "zero_Bi_deriv",  rb_gsl_sf_airy_zero_Bi_deriv, 1);
  rb_define_module_function(mgsl_sf_airy, "zero_Bi_deriv_e",  rb_gsl_sf_airy_zero_Bi_deriv_e, 1);
}
