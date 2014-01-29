/*
  sf_expint.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "include/rb_gsl_sf.h"

static VALUE rb_gsl_sf_expint_E1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_expint_E1, x);
}

static VALUE rb_gsl_sf_expint_E1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_expint_E1_e, x);
}

static VALUE rb_gsl_sf_expint_E2(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_expint_E2, x);
}

static VALUE rb_gsl_sf_expint_E2_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_expint_E2_e, x);
}

static VALUE rb_gsl_sf_expint_Ei(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_expint_Ei, x);
}

static VALUE rb_gsl_sf_expint_Ei_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_expint_Ei_e, x);
}

static VALUE rb_gsl_sf_Shi(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_Shi, x);
}

static VALUE rb_gsl_sf_Shi_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_Shi_e, x);
}

static VALUE rb_gsl_sf_Chi(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_Chi, x);
}

static VALUE rb_gsl_sf_Chi_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_Chi_e, x);
}

static VALUE rb_gsl_sf_expint_3(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_expint_3, x);
}

static VALUE rb_gsl_sf_expint_3_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_expint_3_e, x);
}

static VALUE rb_gsl_sf_Si(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_Si, x);
}

static VALUE rb_gsl_sf_Si_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_Si_e, x);
}

static VALUE rb_gsl_sf_Ci(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_Ci, x);
}

static VALUE rb_gsl_sf_Ci_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_Ci_e, x);
}

static VALUE rb_gsl_sf_atanint(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_atanint, x);
}

static VALUE rb_gsl_sf_atanint_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_atanint_e, x);
}

#ifdef GSL_1_3_LATER
static VALUE rb_gsl_sf_expint_E1_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_expint_E1_scaled, x);
}

static VALUE rb_gsl_sf_expint_E1_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_expint_E1_scaled_e, x);
}

static VALUE rb_gsl_sf_expint_E2_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_expint_E2_scaled, x);
}

static VALUE rb_gsl_sf_expint_E2_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_expint_E2_scaled_e, x);
}

static VALUE rb_gsl_sf_expint_Ei_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_expint_Ei_scaled, x);
}

static VALUE rb_gsl_sf_expint_Ei_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_expint_Ei_scaled_e, x);
}
#endif

#ifdef GSL_1_10_LATER
static VALUE rb_gsl_sf_expint_En(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_expint_En, n, x);
}
static VALUE rb_gsl_sf_expint_En_e(VALUE obj, VALUE n, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE val;
  val = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  gsl_sf_expint_En_e(FIX2INT(n), NUM2DBL(x), rslt);
  return val;
}

#endif

void Init_gsl_sf_expint(VALUE module)
{
  VALUE mgsl_sf_expint;

  rb_define_module_function(module, "expint_E1",  rb_gsl_sf_expint_E1, 1);
  rb_define_module_function(module, "expint_E1_e",  rb_gsl_sf_expint_E1_e, 1);
  rb_define_module_function(module, "expint_E2",  rb_gsl_sf_expint_E2, 1);
  rb_define_module_function(module, "expint_E2_e",  rb_gsl_sf_expint_E2_e, 1);
  rb_define_module_function(module, "expint_Ei",  rb_gsl_sf_expint_Ei, 1);
  rb_define_module_function(module, "expint_Ei_e",  rb_gsl_sf_expint_Ei_e, 1);

  rb_define_module_function(module, "Shi",  rb_gsl_sf_Shi, 1);
  rb_define_module_function(module, "Shi_e",  rb_gsl_sf_Shi_e, 1);
  rb_define_module_function(module, "Chi",  rb_gsl_sf_Chi, 1);
  rb_define_module_function(module, "Chi_e",  rb_gsl_sf_Chi_e, 1);
  rb_define_module_function(module, "expint_3",  rb_gsl_sf_expint_3, 1);
  rb_define_module_function(module, "expint_3_e",  rb_gsl_sf_expint_3_e, 1);
  rb_define_module_function(module, "Si",  rb_gsl_sf_Si, 1);
  rb_define_module_function(module, "Si_e",  rb_gsl_sf_Si_e, 1);
  rb_define_module_function(module, "Ci",  rb_gsl_sf_Ci, 1);
  rb_define_module_function(module, "Ci_e",  rb_gsl_sf_Ci_e, 1);
  rb_define_module_function(module, "atanint",  rb_gsl_sf_atanint, 1);
  rb_define_module_function(module, "atanint_e",  rb_gsl_sf_atanint_e, 1);


  mgsl_sf_expint = rb_define_module_under(module, "Expint");
  rb_define_module_function(mgsl_sf_expint, "E1",  rb_gsl_sf_expint_E1, 1);
  rb_define_module_function(mgsl_sf_expint, "E1_e",  rb_gsl_sf_expint_E1_e, 1);
  rb_define_module_function(mgsl_sf_expint, "E2",  rb_gsl_sf_expint_E2, 1);
  rb_define_module_function(mgsl_sf_expint, "E2_e",  rb_gsl_sf_expint_E2_e, 1);
  rb_define_module_function(mgsl_sf_expint, "Ei",  rb_gsl_sf_expint_Ei, 1);
  rb_define_module_function(mgsl_sf_expint, "Ei_e",  rb_gsl_sf_expint_Ei_e, 1);
  rb_define_module_function(mgsl_sf_expint, "three",  rb_gsl_sf_expint_3, 1);
  rb_define_module_function(mgsl_sf_expint, "three_e",  rb_gsl_sf_expint_3_e, 1);

#ifdef GSL_1_3_LATER
  rb_define_module_function(module, "expint_E1_scaled",  rb_gsl_sf_expint_E1_scaled, 1);
  rb_define_module_function(module, "expint_E1_scaled_e",  rb_gsl_sf_expint_E1_scaled_e, 1);
  rb_define_module_function(module, "expint_E2_scaled",  rb_gsl_sf_expint_E2_scaled, 1);
  rb_define_module_function(module, "expint_E2_scaled_e",  rb_gsl_sf_expint_E2_scaled_e, 1);
  rb_define_module_function(module, "expint_Ei_scaled",  rb_gsl_sf_expint_Ei_scaled, 1);
  rb_define_module_function(module, "expint_Ei_scaled_e",  rb_gsl_sf_expint_Ei_scaled_e, 1);

  rb_define_module_function(mgsl_sf_expint, "E1_scaled",  rb_gsl_sf_expint_E1_scaled, 1);
  rb_define_module_function(mgsl_sf_expint, "E1_scaled_e",  rb_gsl_sf_expint_E1_scaled_e, 1);
  rb_define_module_function(mgsl_sf_expint, "E2_scaled",  rb_gsl_sf_expint_E2_scaled, 1);
  rb_define_module_function(mgsl_sf_expint, "E2_scaled_e",  rb_gsl_sf_expint_E2_scaled_e, 1);
  rb_define_module_function(mgsl_sf_expint, "Ei_scaled",  rb_gsl_sf_expint_Ei_scaled, 1);
  rb_define_module_function(mgsl_sf_expint, "Ei_scaled_e",  rb_gsl_sf_expint_Ei_scaled_e, 1);
#endif

#ifdef GSL_1_10_LATER
  rb_define_module_function(module, "expint_En",  rb_gsl_sf_expint_En, 2);
  rb_define_module_function(mgsl_sf_expint, "En",  rb_gsl_sf_expint_En, 2);
  rb_define_module_function(module, "expint_En_e",  rb_gsl_sf_expint_En_e, 2);
  rb_define_module_function(mgsl_sf_expint, "En_e",  rb_gsl_sf_expint_En_e, 2);  
#endif

}
