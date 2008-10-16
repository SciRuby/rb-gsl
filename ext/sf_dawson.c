/*
  sf_dawson.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "rb_gsl_sf.h"

static VALUE rb_gsl_sf_dawson(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_dawson, x);
}

static VALUE rb_gsl_sf_dawson_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_dawson_e, x);
}

void Init_gsl_sf_dawson(VALUE module)
{
  rb_define_module_function(module, "dawson",  rb_gsl_sf_dawson, 1);
  rb_define_module_function(module, "dawson_e",  rb_gsl_sf_dawson_e, 1);
}
