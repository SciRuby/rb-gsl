/*
  sf_lambert.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_sf.h"

static VALUE rb_gsl_sf_lambert_W0(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_lambert_W0, x);
}

static VALUE rb_gsl_sf_lambert_W0_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_lambert_W0_e, x);
}

static VALUE rb_gsl_sf_lambert_Wm1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_lambert_Wm1, x);
}

static VALUE rb_gsl_sf_lambert_Wm1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_lambert_Wm1_e, x);
}

void Init_gsl_sf_lambert(VALUE module)
{
  VALUE mgsl_sf_lambert;
  rb_define_module_function(module, "lambert_W0",  rb_gsl_sf_lambert_W0, 1);
  rb_define_module_function(module, "lambert_W0_e",  rb_gsl_sf_lambert_W0_e, 1);
  rb_define_module_function(module, "lambert_Wm1",  rb_gsl_sf_lambert_Wm1, 1);
  rb_define_module_function(module, "lambert_Wm1_e",  rb_gsl_sf_lambert_Wm1_e, 1);

  mgsl_sf_lambert = rb_define_module_under(module, "Lambert");
  rb_define_module_function(mgsl_sf_lambert, "W0",  rb_gsl_sf_lambert_W0, 1);
  rb_define_module_function(mgsl_sf_lambert, "W0_e",  rb_gsl_sf_lambert_W0_e, 1);
  rb_define_module_function(mgsl_sf_lambert, "Wm1",  rb_gsl_sf_lambert_Wm1, 1);
  rb_define_module_function(mgsl_sf_lambert, "Wm1_e",  rb_gsl_sf_lambert_Wm1_e, 1);
}
