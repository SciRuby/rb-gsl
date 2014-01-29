/*
  sf_psi.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_sf.h"

static VALUE rb_gsl_sf_psi_int(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval1_int(gsl_sf_psi_int, n);
}

static VALUE rb_gsl_sf_psi_int_e(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval_e_int(gsl_sf_psi_int_e, n);
}

static VALUE rb_gsl_sf_psi(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_psi, x);
}

static VALUE rb_gsl_sf_psi_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_psi_e, x);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_sf_psi_1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_psi_1, x);
}
#endif

#ifdef GSL_1_4_9_LATER
static VALUE rb_gsl_sf_psi_1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_psi_1_e, x);
}
#endif

static VALUE rb_gsl_sf_psi_1piy(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_psi_1piy, x);
}

static VALUE rb_gsl_sf_psi_1piy_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_psi_1piy_e, x);
}

static VALUE rb_gsl_sf_psi_1_int(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval1_int(gsl_sf_psi_1_int, n);
}

static VALUE rb_gsl_sf_psi_1_int_e(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval_e_int(gsl_sf_psi_1_int_e, n);
}

static VALUE rb_gsl_sf_psi_n(VALUE obj, VALUE m, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_psi_n, m, x);
}

static VALUE rb_gsl_sf_psi_n_e(VALUE obj, VALUE m, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_psi_n_e, m, x);
}

void Init_gsl_sf_psi(VALUE module)
{
  rb_define_module_function(module, "psi_int",  rb_gsl_sf_psi_int, 1);
  rb_define_module_function(module, "psi_int_e",  rb_gsl_sf_psi_int_e, 1);
  rb_define_module_function(module, "psi_1piy",  rb_gsl_sf_psi_1piy, 1);
  rb_define_module_function(module, "psi_1piy_e",  rb_gsl_sf_psi_1piy_e, 1);
  rb_define_module_function(module, "psi_1_int",  rb_gsl_sf_psi_1_int, 1);
  rb_define_module_function(module, "psi_1_int_e",  rb_gsl_sf_psi_1_int_e, 1);
  rb_define_module_function(module, "psi_n",  rb_gsl_sf_psi_n, 2);
  rb_define_module_function(module, "psi_n_e",  rb_gsl_sf_psi_n_e, 2);

  rb_define_module_function(module, "psi",  rb_gsl_sf_psi, 1);
  rb_define_module_function(module, "psi_e",  rb_gsl_sf_psi_e, 1);

#ifdef GSL_1_6_LATER
    rb_define_module_function(module, "psi_1",  rb_gsl_sf_psi_1, 1);
#endif
#ifdef GSL_1_4_9_LATER
    rb_define_module_function(module, "psi_1_e",  rb_gsl_sf_psi_1_e, 1);
#endif
}
