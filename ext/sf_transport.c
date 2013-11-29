/*
  sf_transport.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_sf.h"

static VALUE rb_gsl_sf_transport_2(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_transport_2, x);
}

static VALUE rb_gsl_sf_transport_2_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_transport_2_e, x);
}

static VALUE rb_gsl_sf_transport_3(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_transport_3, x);
}

static VALUE rb_gsl_sf_transport_3_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_transport_3_e, x);
}

static VALUE rb_gsl_sf_transport_4(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_transport_4, x);
}

static VALUE rb_gsl_sf_transport_4_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_transport_4_e, x);
}

static VALUE rb_gsl_sf_transport_5(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_transport_5, x);
}

static VALUE rb_gsl_sf_transport_5_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_transport_5_e, x);
}

void Init_gsl_sf_transport(VALUE module)
{
  VALUE mgsl_sf_trans;

  rb_define_module_function(module, "transport_2",  rb_gsl_sf_transport_2, 1);
  rb_define_module_function(module, "transport_2_e",  rb_gsl_sf_transport_2_e, 1);
  rb_define_module_function(module, "transport_3",  rb_gsl_sf_transport_3, 1);
  rb_define_module_function(module, "transport_3_e",  rb_gsl_sf_transport_3_e, 1);
  rb_define_module_function(module, "transport_4",  rb_gsl_sf_transport_4, 1);
  rb_define_module_function(module, "transport_4_e",  rb_gsl_sf_transport_4_e, 1);
  rb_define_module_function(module, "transport_5",  rb_gsl_sf_transport_5, 1);
  rb_define_module_function(module, "transport_5_e",  rb_gsl_sf_transport_5_e, 1);

  mgsl_sf_trans = rb_define_module_under(module, "Transport");
  rb_define_module_function(mgsl_sf_trans, "two",  rb_gsl_sf_transport_2, 1);
  rb_define_module_function(mgsl_sf_trans, "two_e",  rb_gsl_sf_transport_2_e, 1);
  rb_define_module_function(mgsl_sf_trans, "three",  rb_gsl_sf_transport_3, 1);
  rb_define_module_function(mgsl_sf_trans, "three_e",  rb_gsl_sf_transport_3_e, 1);
  rb_define_module_function(mgsl_sf_trans, "four",  rb_gsl_sf_transport_4, 1);
  rb_define_module_function(mgsl_sf_trans, "four_e",  rb_gsl_sf_transport_4_e, 1);
  rb_define_module_function(mgsl_sf_trans, "five",  rb_gsl_sf_transport_5, 1);
  rb_define_module_function(mgsl_sf_trans, "fine_e",  rb_gsl_sf_transport_5_e, 1);
}
