/*
  sf_clausen.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_sf.h"

static VALUE rb_gsl_sf_clausen(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_clausen, x);
}

static VALUE rb_gsl_sf_clausen_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_clausen_e, x);
}

void Init_gsl_sf_clausen(VALUE module)
{
  rb_define_module_function(module, "clausen",  rb_gsl_sf_clausen, 1);
  rb_define_module_function(module, "clausen_e",  rb_gsl_sf_clausen_e, 1);
}
