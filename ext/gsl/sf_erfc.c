/*
  sf_erfc.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "include/rb_gsl_sf.h"

static VALUE rb_gsl_sf_erf(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_erf, x);
}

static VALUE rb_gsl_sf_erf_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_erf_e, x);
}

static VALUE rb_gsl_sf_erfc(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_erfc, x);
}

static VALUE rb_gsl_sf_erfc_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_erfc_e, x);
}

static VALUE rb_gsl_sf_log_erfc(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_log_erfc, x);
}

static VALUE rb_gsl_sf_log_erfc_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_log_erfc_e, x);
}

static VALUE rb_gsl_sf_erf_Z(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_erf_Z, x);
}

static VALUE rb_gsl_sf_erf_Z_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_erf_Z_e, x);
}

static VALUE rb_gsl_sf_erf_Q(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_erf_Q, x);
}

static VALUE rb_gsl_sf_erf_Q_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_erf_Q_e, x);
}

#ifdef GSL_1_4_LATER
static VALUE rb_gsl_sf_hazard(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_hazard, x);
}

static VALUE rb_gsl_sf_hazard_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_hazard_e, x);
}
#endif

void Init_gsl_sf_erfc(VALUE module)
{
  rb_define_module_function(module, "erf",  rb_gsl_sf_erf, 1);
  rb_define_module_function(module, "erf_e",  rb_gsl_sf_erf_e, 1);
  rb_define_module_function(module, "erfc",  rb_gsl_sf_erfc, 1);
  rb_define_module_function(module, "erfc_e",  rb_gsl_sf_erfc_e, 1);
  rb_define_module_function(module, "log_erfc",  rb_gsl_sf_log_erfc, 1);
  rb_define_module_function(module, "log_erfc_e",  rb_gsl_sf_log_erfc_e, 1);
  rb_define_module_function(module, "erf_Z",  rb_gsl_sf_erf_Z, 1);
  rb_define_module_function(module, "erf_Z_e",  rb_gsl_sf_erf_Z_e, 1);
  rb_define_module_function(module, "erf_Q",  rb_gsl_sf_erf_Q, 1);
  rb_define_module_function(module, "erf_Q_e",  rb_gsl_sf_erf_Q_e, 1);
#ifdef GSL_1_4_LATER
  rb_define_module_function(module, "hazard",  rb_gsl_sf_hazard, 1);
  rb_define_module_function(module, "hazard_e",  rb_gsl_sf_hazard_e, 1);
#endif
}
