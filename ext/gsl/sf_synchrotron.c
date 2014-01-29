/*
  sf_synchrotron.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_sf.h"

static VALUE rb_gsl_sf_synchrotron_1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_synchrotron_1, x);
}

static VALUE rb_gsl_sf_synchrotron_1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_synchrotron_1_e, x);
}

static VALUE rb_gsl_sf_synchrotron_2(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_synchrotron_2, x);
}

static VALUE rb_gsl_sf_synchrotron_2_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_synchrotron_2_e, x);
}

void Init_gsl_sf_synchrotron(VALUE module)
{
  VALUE mgsl_sf_synch;

  rb_define_module_function(module, "synchrotron_1",  rb_gsl_sf_synchrotron_1, 1);
  rb_define_module_function(module, "synchrotron_1_e",  rb_gsl_sf_synchrotron_1_e, 1);
  rb_define_module_function(module, "synchrotron_2",  rb_gsl_sf_synchrotron_2, 1);
  rb_define_module_function(module, "synchrotron_2_e",  rb_gsl_sf_synchrotron_2_e, 1);

  mgsl_sf_synch = rb_define_module_under(module, "Synchrotron");
  rb_define_module_function(mgsl_sf_synch, "one",  rb_gsl_sf_synchrotron_1, 1);
  rb_define_module_function(mgsl_sf_synch, "one_e",  rb_gsl_sf_synchrotron_1_e, 1);
  rb_define_module_function(mgsl_sf_synch, "two",  rb_gsl_sf_synchrotron_2, 1);
  rb_define_module_function(mgsl_sf_synch, "two_e",  rb_gsl_sf_synchrotron_2_e, 1);
}
