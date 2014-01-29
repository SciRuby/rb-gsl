/*
  sf_debye.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "include/rb_gsl_sf.h"

static VALUE rb_gsl_sf_debye_1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_debye_1, x);
}

static VALUE rb_gsl_sf_debye_1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_debye_1_e, x);
}

static VALUE rb_gsl_sf_debye_2(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_debye_2, x);
}

static VALUE rb_gsl_sf_debye_2_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_debye_2_e, x);
}

static VALUE rb_gsl_sf_debye_3(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_debye_3, x);
}

static VALUE rb_gsl_sf_debye_3_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_debye_3_e, x);
}

static VALUE rb_gsl_sf_debye_4(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_debye_4, x);
}

static VALUE rb_gsl_sf_debye_4_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_debye_4_e, x);
}

static VALUE rb_gsl_sf_debye_n(int argc, VALUE *argv, VALUE obj)
{
  int n;
  VALUE x;
  switch (argc) {
  case 1:
    n = 1;
    x = argv[0];
    break;
  case 2:
    CHECK_FIXNUM(argv[0]);
    n = FIX2INT(argv[0]);
    x = argv[1];
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  switch (n) {
  case 1:
    return rb_gsl_sf_eval1(gsl_sf_debye_1, x);
    break;
  case 2:
    return rb_gsl_sf_eval1(gsl_sf_debye_2, x);
    break;
  case 3:
    return rb_gsl_sf_eval1(gsl_sf_debye_3, x);
    break;
  case 4:
    return rb_gsl_sf_eval1(gsl_sf_debye_4, x);
    break;
#ifdef GSL_1_8_LATER
  case 5:
    return rb_gsl_sf_eval1(gsl_sf_debye_5, x);
    break;
  case 6:
    return rb_gsl_sf_eval1(gsl_sf_debye_6, x);
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "n must be 1, 2, 3, or 4");
    break;
  }
}

#ifdef GSL_1_8_LATER
static VALUE rb_gsl_sf_debye_5(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_debye_5, x);
}

static VALUE rb_gsl_sf_debye_5_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_debye_5_e, x);
}
static VALUE rb_gsl_sf_debye_6(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_debye_6, x);
}

static VALUE rb_gsl_sf_debye_6_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_debye_6_e, x);
}

#endif

void Init_gsl_sf_debye(VALUE module)
{
  VALUE mgsl_sf_debye;
  rb_define_module_function(module, "debye_1",  rb_gsl_sf_debye_1, 1);
  rb_define_module_function(module, "debye_1_e",  rb_gsl_sf_debye_1_e, 1);
  rb_define_module_function(module, "debye_2",  rb_gsl_sf_debye_2, 1);
  rb_define_module_function(module, "debye_2_e",  rb_gsl_sf_debye_2_e, 1);
  rb_define_module_function(module, "debye_3",  rb_gsl_sf_debye_3, 1);
  rb_define_module_function(module, "debye_3_e",  rb_gsl_sf_debye_3_e, 1);
  rb_define_module_function(module, "debye_4",  rb_gsl_sf_debye_4, 1);
  rb_define_module_function(module, "debye_4_e",  rb_gsl_sf_debye_4_e, 1);
#ifdef GSL_1_8_LATER
  rb_define_module_function(module, "debye_5",  rb_gsl_sf_debye_5, 1);
  rb_define_module_function(module, "debye_5_e",  rb_gsl_sf_debye_5_e, 1);
  rb_define_module_function(module, "debye_6",  rb_gsl_sf_debye_6, 1);
  rb_define_module_function(module, "debye_6_e",  rb_gsl_sf_debye_6_e, 1);
#endif
  rb_define_module_function(module, "debye_n",  rb_gsl_sf_debye_n, -1);

  mgsl_sf_debye = rb_define_module_under(module, "Debye");
  rb_define_module_function(mgsl_sf_debye, "one",  rb_gsl_sf_debye_1, 1);
  rb_define_module_function(mgsl_sf_debye, "one_e",  rb_gsl_sf_debye_1_e, 1);
  rb_define_module_function(mgsl_sf_debye, "two",  rb_gsl_sf_debye_2, 1);
  rb_define_module_function(mgsl_sf_debye, "two_e",  rb_gsl_sf_debye_2_e, 1);
  rb_define_module_function(mgsl_sf_debye, "three",  rb_gsl_sf_debye_3, 1);
  rb_define_module_function(mgsl_sf_debye, "three_e",  rb_gsl_sf_debye_3_e, 1);
  rb_define_module_function(mgsl_sf_debye, "four",  rb_gsl_sf_debye_4, 1);
  rb_define_module_function(mgsl_sf_debye, "four_e",  rb_gsl_sf_debye_4_e, 1);
#ifdef GSL_1_8_LATER
  rb_define_module_function(mgsl_sf_debye, "five",  rb_gsl_sf_debye_5, 1);
  rb_define_module_function(mgsl_sf_debye, "five_e",  rb_gsl_sf_debye_5_e, 1);
  rb_define_module_function(mgsl_sf_debye, "six",  rb_gsl_sf_debye_6, 1);
  rb_define_module_function(mgsl_sf_debye, "six_e",  rb_gsl_sf_debye_6_e, 1);
#endif
  rb_define_module_function(mgsl_sf_debye, "n",  rb_gsl_sf_debye_n, -1);
}
