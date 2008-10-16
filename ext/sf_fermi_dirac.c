/*
  sf_fermi_dirac.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "rb_gsl_sf.h"

static VALUE rb_gsl_sf_fermi_dirac_m1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_fermi_dirac_m1, x);
}

static VALUE rb_gsl_sf_fermi_dirac_m1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_fermi_dirac_m1_e, x);
}

static VALUE rb_gsl_sf_fermi_dirac_0(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_fermi_dirac_0, x);
}

static VALUE rb_gsl_sf_fermi_dirac_0_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_fermi_dirac_0_e, x);
}

static VALUE rb_gsl_sf_fermi_dirac_1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_fermi_dirac_1, x);
}

static VALUE rb_gsl_sf_fermi_dirac_1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_fermi_dirac_1_e, x);
}

static VALUE rb_gsl_sf_fermi_dirac_2(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_fermi_dirac_2, x);
}

static VALUE rb_gsl_sf_fermi_dirac_2_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_fermi_dirac_2_e, x);
}

static VALUE rb_gsl_sf_fermi_dirac_int(VALUE obj, VALUE j, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_fermi_dirac_int, j, x);
}

static VALUE rb_gsl_sf_fermi_dirac_int_e(VALUE obj, VALUE j, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_fermi_dirac_int_e, j, x);
}

static VALUE rb_gsl_sf_fermi_dirac_mhalf(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_fermi_dirac_mhalf, x);
}

static VALUE rb_gsl_sf_fermi_dirac_mhalf_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_fermi_dirac_mhalf_e, x);
}

static VALUE rb_gsl_sf_fermi_dirac_half(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_fermi_dirac_half, x);
}

static VALUE rb_gsl_sf_fermi_dirac_half_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_fermi_dirac_half_e, x);
}

static VALUE rb_gsl_sf_fermi_dirac_3half(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_fermi_dirac_3half, x);
}

static VALUE rb_gsl_sf_fermi_dirac_3half_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_fermi_dirac_3half_e, x);
}

static VALUE rb_gsl_sf_fermi_dirac_inc_0(VALUE obj, VALUE x, VALUE b)
{
  return rb_float_new(gsl_sf_fermi_dirac_inc_0(NUM2DBL(x), NUM2DBL(b)));
}

static VALUE rb_gsl_sf_fermi_dirac_inc_0_e(VALUE obj, VALUE x, VALUE b)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_fermi_dirac_inc_0_e, x, b);
}

void Init_gsl_sf_fermi_dirac(VALUE module)
{
  VALUE mgsl_sf_fermi;

  rb_define_module_function(module, "fermi_dirac_m1",  rb_gsl_sf_fermi_dirac_m1, 1);
  rb_define_module_function(module, "fermi_dirac_m1_e",  rb_gsl_sf_fermi_dirac_m1_e, 1);
  rb_define_module_function(module, "fermi_dirac_0",  rb_gsl_sf_fermi_dirac_0, 1);
  rb_define_module_function(module, "fermi_dirac_0_e",  rb_gsl_sf_fermi_dirac_0_e, 1);
  rb_define_module_function(module, "fermi_dirac_1",  rb_gsl_sf_fermi_dirac_1, 1);
  rb_define_module_function(module, "fermi_dirac_1_e",  rb_gsl_sf_fermi_dirac_1_e, 1);
  rb_define_module_function(module, "fermi_dirac_2",  rb_gsl_sf_fermi_dirac_2, 1);
  rb_define_module_function(module, "fermi_dirac_2_e",  rb_gsl_sf_fermi_dirac_2_e, 1);
  rb_define_module_function(module, "fermi_dirac_int",  rb_gsl_sf_fermi_dirac_int, 2);
  rb_define_module_function(module, "fermi_dirac_int_e",  rb_gsl_sf_fermi_dirac_int_e, 2);
  rb_define_module_function(module, "fermi_dirac_mhalf",  rb_gsl_sf_fermi_dirac_mhalf, 1);
  rb_define_module_function(module, "fermi_dirac_mhalf_e",  rb_gsl_sf_fermi_dirac_mhalf_e, 1);
  rb_define_module_function(module, "fermi_dirac_half",  rb_gsl_sf_fermi_dirac_half, 1);
  rb_define_module_function(module, "fermi_dirac_half_e",  rb_gsl_sf_fermi_dirac_half_e, 1);
  rb_define_module_function(module, "fermi_dirac_3half",  rb_gsl_sf_fermi_dirac_3half, 1);
  rb_define_module_function(module, "fermi_dirac_3half_e",  rb_gsl_sf_fermi_dirac_3half_e, 1);
  rb_define_module_function(module, "fermi_dirac_inc_0",  rb_gsl_sf_fermi_dirac_inc_0, 2);
  rb_define_module_function(module, "fermi_dirac_inc_0_e",  rb_gsl_sf_fermi_dirac_inc_0_e, 2);

  mgsl_sf_fermi = rb_define_module_under(module, "Fermi_Dirac");
  rb_define_module_function(mgsl_sf_fermi, "m1",  rb_gsl_sf_fermi_dirac_m1, 1);
  rb_define_module_function(mgsl_sf_fermi, "m1_e",  rb_gsl_sf_fermi_dirac_m1_e, 1);
  rb_define_module_function(mgsl_sf_fermi, "zero",  rb_gsl_sf_fermi_dirac_0, 1);
  rb_define_module_function(mgsl_sf_fermi, "zero_e",  rb_gsl_sf_fermi_dirac_0_e, 1);
  rb_define_module_function(mgsl_sf_fermi, "one",  rb_gsl_sf_fermi_dirac_1, 1);
  rb_define_module_function(mgsl_sf_fermi, "one_e",  rb_gsl_sf_fermi_dirac_1_e, 1);
  rb_define_module_function(mgsl_sf_fermi, "two",  rb_gsl_sf_fermi_dirac_2, 1);
  rb_define_module_function(mgsl_sf_fermi, "two_e",  rb_gsl_sf_fermi_dirac_2_e, 1);
  rb_define_module_function(mgsl_sf_fermi, "int",  rb_gsl_sf_fermi_dirac_int, 2);
  rb_define_module_function(mgsl_sf_fermi, "int_e",  rb_gsl_sf_fermi_dirac_int_e, 2);
  rb_define_module_function(mgsl_sf_fermi, "mhalf",  rb_gsl_sf_fermi_dirac_mhalf, 1);
  rb_define_module_function(mgsl_sf_fermi, "mhalf_e",  rb_gsl_sf_fermi_dirac_mhalf_e, 1);
  rb_define_module_function(mgsl_sf_fermi, "half",  rb_gsl_sf_fermi_dirac_half, 1);
  rb_define_module_function(mgsl_sf_fermi, "half_e",  rb_gsl_sf_fermi_dirac_half_e, 1);
  rb_define_module_function(mgsl_sf_fermi, "threehalf",  rb_gsl_sf_fermi_dirac_3half, 1);
  rb_define_module_function(mgsl_sf_fermi, "threehalf_e",  rb_gsl_sf_fermi_dirac_3half_e, 1);
  rb_define_module_function(mgsl_sf_fermi, "inc_0",  rb_gsl_sf_fermi_dirac_inc_0, 2);
  rb_define_module_function(mgsl_sf_fermi, "inc_0_e",  rb_gsl_sf_fermi_dirac_inc_0_e, 2);

}
