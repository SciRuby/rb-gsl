/*
  sf_gamma.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_sf.h"

static VALUE rb_gsl_sf_gamma(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_gamma, x);
}

static VALUE rb_gsl_sf_gamma_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_gamma_e, x);
}

static VALUE rb_gsl_sf_lngamma(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_lngamma, x);
}

static VALUE rb_gsl_sf_lngamma_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_lngamma_e, x);
}

static VALUE rb_gsl_sf_lngamma_sgn_e(VALUE obj, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  double sgn;
  Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_lngamma_sgn_e(NUM2DBL(x), rslt, &sgn);
  return rb_ary_new3(2, v, rb_float_new(sgn));
}

static VALUE rb_gsl_sf_gammastar(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_gammastar, x);
}

static VALUE rb_gsl_sf_gammastar_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_gammastar_e, x);
}

static VALUE rb_gsl_sf_gammainv(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_gammainv, x);
}

static VALUE rb_gsl_sf_gammainv_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_gammainv_e, x);
}

static VALUE rb_gsl_sf_lngamma_complex_e(int argc, VALUE *argv, VALUE obj)
{
  gsl_sf_result *lnr, *arg;
  gsl_complex *z;
  double re, im;
  VALUE vlnr, varg;
  int status;
  switch (argc) {
  case 1:
    CHECK_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_complex, z);
    re = GSL_REAL(*z);
    im = GSL_IMAG(*z);
    break;
  case 2:
    Need_Float(argv[0]); Need_Float(argv[1]);
    re = NUM2DBL(argv[0]);
    im = NUM2DBL(argv[1]);
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  }
  vlnr = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, lnr);
  varg = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, arg);
  status = gsl_sf_lngamma_complex_e(re, im, lnr, arg);
  return rb_ary_new3(3, vlnr, varg, INT2FIX(status));
}

static VALUE rb_gsl_sf_taylorcoeff(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_taylorcoeff, n, x);
}

static VALUE rb_gsl_sf_taylorcoeff_e(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_taylorcoeff_e, n, x);
}

static VALUE rb_gsl_sf_fact(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval1_uint(gsl_sf_fact, n);
}

static VALUE rb_gsl_sf_fact_e(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval_e_uint(gsl_sf_fact_e, n);
}

static VALUE rb_gsl_sf_doublefact(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval1_uint(gsl_sf_doublefact, n);
}

static VALUE rb_gsl_sf_doublefact_e(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval_e_uint(gsl_sf_doublefact_e, n);
}

static VALUE rb_gsl_sf_lnfact(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval1_uint(gsl_sf_lnfact, n);
}

static VALUE rb_gsl_sf_lnfact_e(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval_e_uint(gsl_sf_lnfact_e, n);
}

static VALUE rb_gsl_sf_lndoublefact(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval1_uint(gsl_sf_lndoublefact, n);
}

static VALUE rb_gsl_sf_lndoublefact_e(VALUE obj, VALUE n)
{
  return rb_gsl_sf_eval_e_uint(gsl_sf_lndoublefact_e, n);
}

static VALUE rb_gsl_sf_choose(VALUE obj, VALUE n, VALUE m)
{
  return rb_float_new(gsl_sf_choose(FIX2INT(n), FIX2INT(m)));
}

static VALUE rb_gsl_sf_choose_e(VALUE obj, VALUE n, VALUE m)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  CHECK_FIXNUM(n); CHECK_FIXNUM(m);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  gsl_sf_choose_e(FIX2INT(n), FIX2INT(m), rslt);
  return v;
}

static VALUE rb_gsl_sf_lnchoose(VALUE obj, VALUE n, VALUE m)
{
  CHECK_FIXNUM(n); CHECK_FIXNUM(m);
  return rb_float_new(gsl_sf_lnchoose(FIX2INT(n), FIX2INT(m)));
}

static VALUE rb_gsl_sf_lnchoose_e(VALUE obj, VALUE n, VALUE m)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  CHECK_FIXNUM(n); CHECK_FIXNUM(m);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_lnchoose_e(FIX2INT(n), FIX2INT(m), rslt);
  return v;
}

static VALUE rb_gsl_sf_poch(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_poch, a, x);
}

static VALUE rb_gsl_sf_poch_e(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_poch_e, a, x);
}

static VALUE rb_gsl_sf_lnpoch(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_lnpoch, a, x);
}

static VALUE rb_gsl_sf_lnpoch_e(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_lnpoch_e, a, x);
}

static VALUE rb_gsl_sf_lnpoch_sgn_e(VALUE obj, VALUE a, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  double sgn;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(a); Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_lnpoch_sgn_e(NUM2DBL(a), NUM2DBL(x), rslt, &sgn);
  return rb_ary_new3(2, v, rb_float_new(sgn));
}

static VALUE rb_gsl_sf_pochrel(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_pochrel, a, x);
}

static VALUE rb_gsl_sf_pochrel_e(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_pochrel_e, a, x);
}

static VALUE rb_gsl_sf_gamma_inc_Q(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_gamma_inc_Q, a, x);
}

static VALUE rb_gsl_sf_gamma_inc_Q_e(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_gamma_inc_Q_e, a, x);
}

static VALUE rb_gsl_sf_gamma_inc_P(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_gamma_inc_P, a, x);
}

static VALUE rb_gsl_sf_gamma_inc_P_e(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_gamma_inc_P_e, a, x);
}

static VALUE rb_gsl_sf_gamma_inc(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_gamma_inc_P, a, x);
}

#ifdef GSL_1_4_LATER
static VALUE rb_gsl_sf_gamma_inc_e(VALUE obj, VALUE a, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_gamma_inc_e, a, x);
}
#endif

static VALUE rb_gsl_sf_beta(VALUE obj, VALUE a, VALUE b)
{
  return rb_float_new(gsl_sf_beta(NUM2DBL(a), NUM2DBL(b)));
}

static VALUE rb_gsl_sf_beta_e(VALUE obj, VALUE a, VALUE b)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_beta_e, a, b);
}

static VALUE rb_gsl_sf_lnbeta(VALUE obj, VALUE a, VALUE b)
{
  Need_Float(a); Need_Float(b);
  return rb_float_new(gsl_sf_lnbeta(NUM2DBL(a), NUM2DBL(b)));
}

static VALUE rb_gsl_sf_lnbeta_e(VALUE obj, VALUE a, VALUE b)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_lnbeta_e, a, b);
}

static VALUE rb_gsl_sf_beta_inc(VALUE obj, VALUE a, VALUE b, VALUE x)
{
  Need_Float(a); Need_Float(b);
  return rb_gsl_sf_eval_double3(gsl_sf_beta_inc, a, b, x);
}

static VALUE rb_gsl_sf_beta_inc_e(VALUE obj, VALUE a, VALUE b, VALUE x)
{
  return rb_gsl_sf_eval_e_double3(gsl_sf_beta_inc_e, a, b, x);
}

double mygsl_binomial_coef(unsigned int n, unsigned int k);
double mygsl_binomial_coef(unsigned int n, unsigned int k)
{
  return floor(0.5 + exp(gsl_sf_lnfact(n) - gsl_sf_lnfact(k) - gsl_sf_lnfact(n-k))); 
}

static VALUE rb_gsl_sf_bincoef(VALUE obj, VALUE n, VALUE k)
{
  CHECK_FIXNUM(n); CHECK_FIXNUM(k);
  return rb_float_new(mygsl_binomial_coef(FIX2UINT(n), FIX2UINT(k)));
}

void Init_gsl_sf_gamma(VALUE module)
{
  rb_define_const(module, "GAMMA_XMAX", NUM2DBL(GSL_SF_GAMMA_XMAX));
  rb_define_module_function(module, "gamma",  rb_gsl_sf_gamma, 1);
  rb_define_module_function(module, "gamma_e",  rb_gsl_sf_gamma_e, 1);
  rb_define_module_function(module, "lngamma",  rb_gsl_sf_lngamma, 1);
  rb_define_module_function(module, "lngamma_e",  rb_gsl_sf_lngamma_e, 1);
  rb_define_module_function(module, "lngamma_sgn_e",  rb_gsl_sf_lngamma_sgn_e, 1);
  rb_define_module_function(module, "gammastar",  rb_gsl_sf_gammastar, 1);
  rb_define_module_function(module, "gammastar_e",  rb_gsl_sf_gammastar_e, 1);
  rb_define_module_function(module, "gammainv",  rb_gsl_sf_gammainv, 1);
  rb_define_module_function(module, "gammainv_e",  rb_gsl_sf_gammainv_e, 1);
  rb_define_module_function(module, "lngamma_complex_e",  rb_gsl_sf_lngamma_complex_e, -1);
  rb_define_module_function(module, "taylorcoeff",  rb_gsl_sf_taylorcoeff, 2);
  rb_define_module_function(module, "taylorcoeff_e",  rb_gsl_sf_taylorcoeff_e, 2);
  rb_define_module_function(module, "fact",  rb_gsl_sf_fact, 1);
  rb_define_module_function(module, "fact_e",  rb_gsl_sf_fact_e, 1);
  rb_define_module_function(module, "doublefact",  rb_gsl_sf_doublefact, 1);
  rb_define_module_function(module, "doublefact_e",  rb_gsl_sf_doublefact_e, 1);
  rb_define_module_function(module, "lnfact",  rb_gsl_sf_lnfact, 1);
  rb_define_module_function(module, "lnfact_e",  rb_gsl_sf_lnfact_e, 1);
  rb_define_module_function(module, "lndoublefact",  rb_gsl_sf_lndoublefact, 1);
  rb_define_module_function(module, "lndoublefact_e",  rb_gsl_sf_lndoublefact_e, 1);
  rb_define_module_function(module, "choose",  rb_gsl_sf_choose, 2);
  rb_define_module_function(module, "choose_e",  rb_gsl_sf_choose_e, 2);
  rb_define_module_function(module, "lnchoose",  rb_gsl_sf_lnchoose, 2);
  rb_define_module_function(module, "lnchoose_e",  rb_gsl_sf_lnchoose_e, 2);
  rb_define_module_function(module, "poch",  rb_gsl_sf_poch, 2);
  rb_define_module_function(module, "poch_e",  rb_gsl_sf_poch_e, 2);
  rb_define_module_function(module, "lnpoch",  rb_gsl_sf_lnpoch, 2);
  rb_define_module_function(module, "lnpoch_e",  rb_gsl_sf_lnpoch_e, 2);
  rb_define_module_function(module, "lnpoch_sgn_e",  rb_gsl_sf_lnpoch_sgn_e, 2);
  rb_define_module_function(module, "pochrel",  rb_gsl_sf_pochrel, 2);
  rb_define_module_function(module, "pochrel_e",  rb_gsl_sf_pochrel_e, 2);
  rb_define_module_function(module, "gamma_inc_P",  rb_gsl_sf_gamma_inc_P, 2);
  rb_define_module_function(module, "gamma_inc_P_e",  rb_gsl_sf_gamma_inc_P_e, 2);
  rb_define_module_function(module, "gamma_inc_Q",  rb_gsl_sf_gamma_inc_Q, 2);
  rb_define_module_function(module, "gamma_inc_Q_e",  rb_gsl_sf_gamma_inc_Q_e, 2);
  rb_define_module_function(module, "gamma_inc",  rb_gsl_sf_gamma_inc, 2);

#ifdef GSL_1_4_LATER
  rb_define_module_function(module, "gamma_inc_e",  rb_gsl_sf_gamma_inc_e, 2);
#endif
  rb_define_module_function(module, "beta",  rb_gsl_sf_beta, 2);
  rb_define_module_function(module, "beta_e",  rb_gsl_sf_beta_e, 2);
  rb_define_module_function(module, "lnbeta",  rb_gsl_sf_lnbeta, 2);
  rb_define_module_function(module, "lnbeta_e",  rb_gsl_sf_lnbeta_e, 2);
  rb_define_module_function(module, "beta_inc",  rb_gsl_sf_beta_inc, 3);
  rb_define_module_function(module, "beta_inc_e",  rb_gsl_sf_beta_inc_e, 3);

  rb_define_module_function(module, "bincoef",  rb_gsl_sf_bincoef, 2);
}
