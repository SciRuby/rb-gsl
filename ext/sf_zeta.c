/*
  sf_zeta.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_sf.h"

static VALUE rb_gsl_sf_zeta_int(VALUE obj, VALUE n)
{
  VALUE nn;
  if (TYPE(n) != T_FIXNUM) nn = INT2FIX(NUM2INT(n));
  else nn = n;
  return rb_gsl_sf_eval1_int(gsl_sf_zeta_int, nn);
}

static VALUE rb_gsl_sf_zeta_int_e(VALUE obj, VALUE n)
{
  VALUE nn;
  if (TYPE(n) != T_FIXNUM) nn = INT2FIX(NUM2INT(n));
  else nn = n;
  return rb_gsl_sf_eval_e_int(gsl_sf_zeta_int_e, nn);
}

static VALUE rb_gsl_sf_zeta(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_zeta, x);
}

static VALUE rb_gsl_sf_zeta_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_zeta_e, x);
}

static VALUE rb_gsl_sf_hzeta(VALUE obj, VALUE s, VALUE q)
{
  Need_Float(s); Need_Float(q);
  return rb_float_new(gsl_sf_hzeta(NUM2DBL(s), NUM2DBL(q)));
}

static VALUE rb_gsl_sf_hzeta_e(VALUE obj, VALUE s, VALUE q)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_hzeta_e, s, q);
}

static VALUE rb_gsl_sf_eta_int(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval1_int(gsl_sf_eta_int, n);
}

static VALUE rb_gsl_sf_eta_int_e(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval_e_int(gsl_sf_eta_int_e, n);
}

static VALUE rb_gsl_sf_eta(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_eta, x);
}

static VALUE rb_gsl_sf_eta_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_eta_e, x);
}

#ifdef GSL_1_4_9_LATER
static VALUE rb_gsl_sf_zetam1_int(VALUE obj, VALUE n)
{
  VALUE nn;
  if (TYPE(n) != T_FIXNUM) nn = INT2FIX(NUM2INT(n));
  else nn = n;
  return rb_gsl_sf_eval1_int(gsl_sf_zetam1_int, nn);
}

static VALUE rb_gsl_sf_zetam1_int_e(VALUE obj, VALUE n)
{
  VALUE nn;
  if (TYPE(n) != T_FIXNUM) nn = INT2FIX(NUM2INT(n));
  else nn = n;
  return rb_gsl_sf_eval_e_int(gsl_sf_zetam1_int_e, nn);
}

static VALUE rb_gsl_sf_zetam1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_zetam1, x);
}

static VALUE rb_gsl_sf_zetam1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_zetam1_e, x);
}
#endif

void Init_gsl_sf_zeta(VALUE module)
{
  rb_define_module_function(module, "zeta_int",  rb_gsl_sf_zeta_int, 1);
  rb_define_module_function(module, "zeta_int_e",  rb_gsl_sf_zeta_int_e, 1);
  rb_define_module_function(module, "zeta",  rb_gsl_sf_zeta, 1);
  rb_define_module_function(module, "zeta_e",  rb_gsl_sf_zeta_e, 1);

  rb_define_module_function(module, "hzeta",  rb_gsl_sf_hzeta, 2);
  rb_define_module_function(module, "hzeta_e",  rb_gsl_sf_hzeta_e, 2);
  rb_define_module_function(module, "eta_int",  rb_gsl_sf_eta_int, 1);
  rb_define_module_function(module, "eta_int_e",  rb_gsl_sf_eta_int_e, 1);
  rb_define_module_function(module, "eta",  rb_gsl_sf_eta, 1);
  rb_define_module_function(module, "eta_e",  rb_gsl_sf_eta_e, 1);

#ifdef GSL_1_4_9_LATER
  rb_define_module_function(module, "zetam1_int",  rb_gsl_sf_zetam1_int, 1);
  rb_define_module_function(module, "zetam1_int_e",  rb_gsl_sf_zetam1_int_e, 1);
  rb_define_module_function(module, "zetam1",  rb_gsl_sf_zetam1, 1);
  rb_define_module_function(module, "zetam1_e",  rb_gsl_sf_zetam1_e, 1);
#endif
}
