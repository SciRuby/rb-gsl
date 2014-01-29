/*
  cdf.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#ifdef GSL_1_4_LATER
#include "include/rb_gsl_array.h"
#include "include/rb_gsl_common.h"
#include "include/rb_gsl_rng.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

VALUE rb_gsl_eval_pdf_cdf(VALUE xx, double (*f)(double));
VALUE rb_gsl_eval_pdf_cdf2(VALUE xx, VALUE aa, double (*f)(double, double));
VALUE rb_gsl_eval_pdf_cdf3(VALUE xx, VALUE aa, VALUE bb, 
			      double (*f)(double, double, double));
VALUE rb_gsl_eval_pdf_cdf2_uint(VALUE xx, VALUE aa, 
				       double (*f)(unsigned int, double));


static VALUE rb_gsl_cdf_gaussian_P(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) {
    return rb_gsl_eval_pdf_cdf(argv[0], gsl_cdf_ugaussian_P);
  } else if (argc == 2) {
    return rb_gsl_eval_pdf_cdf2(argv[0], argv[1], gsl_cdf_gaussian_P);
  } else {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  }
}
static VALUE rb_gsl_cdf_gaussian_Q(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) {
    return rb_gsl_eval_pdf_cdf(argv[0], gsl_cdf_ugaussian_Q);
  } else if (argc == 2) {
    return rb_gsl_eval_pdf_cdf2(argv[0], argv[1], gsl_cdf_gaussian_Q);
  } else {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  }
}

static VALUE rb_gsl_cdf_gaussian_Pinv(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) {
    return rb_gsl_eval_pdf_cdf(argv[0], gsl_cdf_ugaussian_Pinv);
  } else if (argc == 2) {
    return rb_gsl_eval_pdf_cdf2(argv[0], argv[1], gsl_cdf_gaussian_Pinv);
  } else {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  }
}

static VALUE rb_gsl_cdf_gaussian_Qinv(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) {
    return rb_gsl_eval_pdf_cdf(argv[0], gsl_cdf_ugaussian_Qinv);
  } else if (argc == 2) {
    return rb_gsl_eval_pdf_cdf2(argv[0], argv[1], gsl_cdf_gaussian_Qinv);
  } else {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  }
}

static VALUE rb_gsl_cdf_exponential_P(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_exponential_P);
}  

static VALUE rb_gsl_cdf_exponential_Q(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_exponential_Q);
}  

static VALUE rb_gsl_cdf_exponential_Pinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_exponential_Pinv);
}  

static VALUE rb_gsl_cdf_exponential_Qinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_exponential_Qinv);
}  

static VALUE rb_gsl_cdf_laplace_P(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_laplace_P);
}  

static VALUE rb_gsl_cdf_laplace_Q(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_laplace_Q);
}  

static VALUE rb_gsl_cdf_laplace_Pinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_laplace_Pinv);
}  

static VALUE rb_gsl_cdf_laplace_Qinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_laplace_Qinv);
}  

static VALUE rb_gsl_cdf_cauchy_P(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_cauchy_P);
}  

static VALUE rb_gsl_cdf_cauchy_Q(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_cauchy_Q);
}  

static VALUE rb_gsl_cdf_cauchy_Pinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_cauchy_Pinv);
}  

static VALUE rb_gsl_cdf_cauchy_Qinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_cauchy_Qinv);
}  

static VALUE rb_gsl_cdf_rayleigh_P(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_rayleigh_P);
}  

static VALUE rb_gsl_cdf_rayleigh_Q(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_rayleigh_Q);
}  

static VALUE rb_gsl_cdf_rayleigh_Pinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_rayleigh_Pinv);
}  

static VALUE rb_gsl_cdf_rayleigh_Qinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_rayleigh_Qinv);
}  

static VALUE rb_gsl_cdf_gamma_P(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_gamma_P);
}  

static VALUE rb_gsl_cdf_gamma_Q(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b,  gsl_cdf_gamma_Q);
}  

static VALUE rb_gsl_cdf_gamma_Pinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b,  gsl_cdf_gamma_Pinv);
}  

static VALUE rb_gsl_cdf_gamma_Qinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b,  gsl_cdf_gamma_Qinv);
}  

static VALUE rb_gsl_cdf_flat_P(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b,  gsl_cdf_flat_P);
}  

static VALUE rb_gsl_cdf_flat_Q(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b,  gsl_cdf_flat_Q);
}  

static VALUE rb_gsl_cdf_flat_Pinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b,  gsl_cdf_flat_Pinv);
}  

static VALUE rb_gsl_cdf_flat_Qinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b,  gsl_cdf_flat_Qinv);
}  

static VALUE rb_gsl_cdf_lognormal_P(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b,  gsl_cdf_lognormal_P);
}  

static VALUE rb_gsl_cdf_lognormal_Q(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b,  gsl_cdf_lognormal_Q);
}  

static VALUE rb_gsl_cdf_lognormal_Pinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b,  gsl_cdf_lognormal_Pinv);
}  

static VALUE rb_gsl_cdf_lognormal_Qinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b,  gsl_cdf_lognormal_Qinv);
}  

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_cdf_exppow_P(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_exppow_P);
}  

static VALUE rb_gsl_cdf_exppow_Q(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_exppow_Q);
}  
#endif

static VALUE rb_gsl_cdf_chisq_P(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_chisq_P);
}  

static VALUE rb_gsl_cdf_chisq_Q(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_chisq_Q);
}  

static VALUE rb_gsl_cdf_chisq_Pinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_chisq_Pinv);
}  

static VALUE rb_gsl_cdf_chisq_Qinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_chisq_Qinv);
}  

static VALUE rb_gsl_cdf_fdist_P(VALUE obj, VALUE x, VALUE nu1, VALUE nu2)
{
  return rb_gsl_eval_pdf_cdf3(x, nu1, nu2, gsl_cdf_fdist_P);
}  

static VALUE rb_gsl_cdf_fdist_Q(VALUE obj, VALUE x, VALUE nu1, VALUE nu2)
{
  return rb_gsl_eval_pdf_cdf3(x, nu1, nu2, gsl_cdf_fdist_Q);
}  

static VALUE rb_gsl_cdf_tdist_P(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_tdist_P);
}  

static VALUE rb_gsl_cdf_tdist_Q(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_tdist_Q);
}  

static VALUE rb_gsl_cdf_tdist_Pinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_tdist_Pinv);
}  

static VALUE rb_gsl_cdf_tdist_Qinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_tdist_Qinv);
}  

static VALUE rb_gsl_cdf_beta_P(VALUE obj, VALUE x, VALUE nu1, VALUE nu2)
{
  return rb_gsl_eval_pdf_cdf3(x, nu1, nu2, gsl_cdf_beta_P);
}  

static VALUE rb_gsl_cdf_beta_Q(VALUE obj, VALUE x, VALUE nu1, VALUE nu2)
{
  return rb_gsl_eval_pdf_cdf3(x, nu1, nu2, gsl_cdf_beta_Q);
}  

static VALUE rb_gsl_cdf_logistic_P(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_logistic_P);
}  

static VALUE rb_gsl_cdf_logistic_Q(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_logistic_Q);
}  

static VALUE rb_gsl_cdf_logistic_Pinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_logistic_Pinv);
}  

static VALUE rb_gsl_cdf_logistic_Qinv(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_cdf_logistic_Qinv);
}  

static VALUE rb_gsl_cdf_pareto_P(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_pareto_P);
}  

static VALUE rb_gsl_cdf_pareto_Q(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_pareto_Q);
}  

static VALUE rb_gsl_cdf_pareto_Pinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_pareto_Pinv);
}  

static VALUE rb_gsl_cdf_pareto_Qinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_pareto_Qinv);
}  

static VALUE rb_gsl_cdf_weibull_P(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_weibull_P);
}  

static VALUE rb_gsl_cdf_weibull_Q(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_weibull_Q);
}  

static VALUE rb_gsl_cdf_weibull_Pinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_weibull_Pinv);
}  

static VALUE rb_gsl_cdf_weibull_Qinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_weibull_Qinv);
}  

static VALUE rb_gsl_cdf_gumbel1_P(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_gumbel1_P);
}  

static VALUE rb_gsl_cdf_gumbel1_Q(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_gumbel1_Q);
}  

static VALUE rb_gsl_cdf_gumbel1_Pinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_gumbel1_Pinv);
}  

static VALUE rb_gsl_cdf_gumbel1_Qinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_gumbel1_Qinv);
}  

static VALUE rb_gsl_cdf_gumbel2_P(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_gumbel2_P);
}  

static VALUE rb_gsl_cdf_gumbel2_Q(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_gumbel2_Q);
}  

static VALUE rb_gsl_cdf_gumbel2_Pinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_gumbel2_Pinv);
}  

static VALUE rb_gsl_cdf_gumbel2_Qinv(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_cdf_gumbel2_Qinv);
}  

#ifdef GSL_1_8_LATER
static VALUE rb_gsl_cdf_binomial_P(VALUE obj, VALUE kk, VALUE pp, VALUE nn)
{
  unsigned int k, n;
  double p;
  k = NUM2UINT(kk);
  n = NUM2UINT(nn);
  p = NUM2DBL(pp);
  return rb_float_new(gsl_cdf_binomial_P(k, p, n));
}
static VALUE rb_gsl_cdf_binomial_Q(VALUE obj, VALUE kk, VALUE pp, VALUE nn)
{
  unsigned int k, n;
  double p;
  k = NUM2UINT(kk);
  n = NUM2UINT(nn);
  p = NUM2DBL(pp);
  return rb_float_new(gsl_cdf_binomial_Q(k, p, n));
}
static VALUE rb_gsl_cdf_poisson_P(VALUE obj, VALUE kk, VALUE mm)
{
  unsigned int k;
  double mu;
  k = NUM2UINT(kk);
  mu = NUM2DBL(mm);
  return rb_float_new(gsl_cdf_poisson_P(k, mu));
}
static VALUE rb_gsl_cdf_poisson_Q(VALUE obj, VALUE kk, VALUE mm)
{
  unsigned int k;
  double mu;
  k = NUM2UINT(kk);
  mu = NUM2DBL(mm);
  return rb_float_new(gsl_cdf_poisson_Q(k, mu));
}
static VALUE rb_gsl_cdf_geometric_P(VALUE obj, VALUE kk, VALUE mm)
{
  unsigned int k;
  double mu;
  k = NUM2UINT(kk);
  mu = NUM2DBL(mm);
  return rb_float_new(gsl_cdf_geometric_P(k, mu));
}
static VALUE rb_gsl_cdf_geometric_Q(VALUE obj, VALUE kk, VALUE mm)
{
  unsigned int k;
  double mu;
  k = NUM2UINT(kk);
  mu = NUM2DBL(mm);
  return rb_float_new(gsl_cdf_geometric_Q(k, mu));
}
static VALUE rb_gsl_cdf_negative_binomial_P(VALUE obj, VALUE kk, VALUE pp, VALUE nn)
{
  unsigned int k;
  double p, n;
  k = NUM2UINT(kk);
  p = NUM2DBL(pp);
  n = NUM2DBL(nn);
  return rb_float_new(gsl_cdf_negative_binomial_P(k, p, n));
}
static VALUE rb_gsl_cdf_negative_binomial_Q(VALUE obj, VALUE kk, VALUE pp, VALUE nn)
{
  unsigned int k;
  double p, n;
  k = NUM2UINT(kk);
  p = NUM2DBL(pp);
  n = NUM2DBL(nn);
  return rb_float_new(gsl_cdf_negative_binomial_Q(k, p, n));
}
static VALUE rb_gsl_cdf_pascal_P(VALUE obj, VALUE kk, VALUE pp, VALUE nn)
{
  unsigned int k, n;
  double p;
  k = NUM2UINT(kk);
  n = NUM2UINT(nn);
  p = NUM2DBL(pp);
  return rb_float_new(gsl_cdf_pascal_P(k, p, n));
}
static VALUE rb_gsl_cdf_pascal_Q(VALUE obj, VALUE kk, VALUE pp, VALUE nn)
{
  unsigned int k, n;
  double p;
  k = NUM2UINT(kk);
  n = NUM2UINT(nn);
  p = NUM2DBL(pp);
  return rb_float_new(gsl_cdf_pascal_Q(k, p, n));
}
static VALUE rb_gsl_cdf_hypergeometric_P(VALUE obj, VALUE kk, VALUE nn1, VALUE nn2, VALUE tt)
{
  unsigned int k, n1, n2, t;
  k = NUM2UINT(kk);
  n1 = NUM2UINT(nn1);
  n2 = NUM2UINT(nn2);
  t = NUM2UINT(tt);
  return rb_float_new(gsl_cdf_hypergeometric_P(k, n1, n2, t));
}
static VALUE rb_gsl_cdf_hypergeometric_Q(VALUE obj, VALUE kk, VALUE nn1, VALUE nn2, VALUE tt)
{
  unsigned int k, n1, n2, t;
  k = NUM2UINT(kk);
  n1 = NUM2UINT(nn1);
  n2 = NUM2UINT(nn2);
  t = NUM2UINT(tt);
  return rb_float_new(gsl_cdf_hypergeometric_Q(k, n1, n2, t));
}
static VALUE rb_gsl_cdf_beta_Pinv(VALUE obj, VALUE pp, VALUE aa, VALUE bb)
{
  double P, a, b;
  P = NUM2DBL(pp);
  a = NUM2DBL(aa);
  b = NUM2DBL(bb);
  return rb_float_new(gsl_cdf_beta_Pinv(P, a, b));
}
static VALUE rb_gsl_cdf_beta_Qinv(VALUE obj, VALUE qq, VALUE aa, VALUE bb)
{
  double Q, a, b;
  Q = NUM2DBL(qq);
  a = NUM2DBL(aa);
  b = NUM2DBL(bb);
  return rb_float_new(gsl_cdf_beta_Qinv(Q, a, b));
}
static VALUE rb_gsl_cdf_fdist_Pinv(VALUE obj, VALUE pp, VALUE aa, VALUE bb)
{
  double P, a, b;
  P = NUM2DBL(pp);
  a = NUM2DBL(aa);
  b = NUM2DBL(bb);
  return rb_float_new(gsl_cdf_fdist_Pinv(P, a, b));
}
static VALUE rb_gsl_cdf_fdist_Qinv(VALUE obj, VALUE qq, VALUE aa, VALUE bb)
{
  double Q, a, b;
  Q = NUM2DBL(qq);
  a = NUM2DBL(aa);
  b = NUM2DBL(bb);
  return rb_float_new(gsl_cdf_fdist_Qinv(Q, a, b));
}
#endif

void Init_gsl_cdf(VALUE module)
{
  VALUE mgsl_cdf;

  mgsl_cdf = rb_define_module_under(module, "Cdf");

  /***** Cumulative distribution functions *****/
  rb_define_module_function(module, "cdf_gaussian_P", rb_gsl_cdf_gaussian_P, -1);
  rb_define_module_function(module, "cdf_ugaussian_P", rb_gsl_cdf_gaussian_P, -1);
  rb_define_module_function(module, "cdf_gaussian_Q", rb_gsl_cdf_gaussian_Q, -1);
  rb_define_module_function(module, "cdf_ugaussian_Q", rb_gsl_cdf_gaussian_Q, -1);
  rb_define_module_function(module, "cdf_gaussian_Pinv", rb_gsl_cdf_gaussian_Pinv, -1);
  rb_define_module_function(module, "cdf_ugaussian_Pinv", rb_gsl_cdf_gaussian_Pinv, -1);
  rb_define_module_function(module, "cdf_gaussian_Qinv", rb_gsl_cdf_gaussian_Qinv, -1);
  rb_define_module_function(module, "cdf_ugaussian_Qinv", rb_gsl_cdf_gaussian_Qinv, -1);

  rb_define_module_function(module, "cdf_exponential_P", rb_gsl_cdf_exponential_P, 2);
 rb_define_module_function(module, "cdf_exponential_Q", rb_gsl_cdf_exponential_Q, 2);
  rb_define_module_function(module, "cdf_exponential_Pinv", rb_gsl_cdf_exponential_Pinv, 2);
  rb_define_module_function(module, "cdf_exponential_Qinv", rb_gsl_cdf_exponential_Qinv, 2);

  rb_define_module_function(module, "cdf_laplace_P", rb_gsl_cdf_laplace_P, 2);
  rb_define_module_function(module, "cdf_laplace_Q", rb_gsl_cdf_laplace_Q, 2);
  rb_define_module_function(module, "cdf_laplace_Pinv", rb_gsl_cdf_laplace_Pinv, 2);
  rb_define_module_function(module, "cdf_laplace_Qinv", rb_gsl_cdf_laplace_Qinv, 2);

  rb_define_module_function(module, "cdf_cauchy_P", rb_gsl_cdf_cauchy_P, 2);
  rb_define_module_function(module, "cdf_cauchy_Q", rb_gsl_cdf_cauchy_Q, 2);
  rb_define_module_function(module, "cdf_cauchy_Pinv", rb_gsl_cdf_cauchy_Pinv, 2);
  rb_define_module_function(module, "cdf_cauchy_Qinv", rb_gsl_cdf_cauchy_Qinv, 2);

  rb_define_module_function(module, "cdf_rayleigh_P", rb_gsl_cdf_rayleigh_P, 2);
  rb_define_module_function(module, "cdf_rayleigh_Q", rb_gsl_cdf_rayleigh_Q, 2);
  rb_define_module_function(module, "cdf_rayleigh_Pinv", rb_gsl_cdf_rayleigh_Pinv, 2);
  rb_define_module_function(module, "cdf_rayleigh_Qinv", rb_gsl_cdf_rayleigh_Qinv, 2);

  rb_define_module_function(module, "cdf_gamma_P", rb_gsl_cdf_gamma_P, 3);
  rb_define_module_function(module, "cdf_gamma_Q", rb_gsl_cdf_gamma_Q, 3);
  rb_define_module_function(module, "cdf_gamma_Pinv", rb_gsl_cdf_gamma_Pinv, 3);
  rb_define_module_function(module, "cdf_gamma_Qinv", rb_gsl_cdf_gamma_Qinv, 3);

  rb_define_module_function(module, "cdf_flat_P", rb_gsl_cdf_flat_P, 3);
  rb_define_module_function(module, "cdf_flat_Q", rb_gsl_cdf_flat_Q, 3);
  rb_define_module_function(module, "cdf_flat_Pinv", rb_gsl_cdf_flat_Pinv, 3);
  rb_define_module_function(module, "cdf_flat_Qinv", rb_gsl_cdf_flat_Qinv, 3);

  rb_define_module_function(module, "cdf_lognormal_P", rb_gsl_cdf_lognormal_P, 3);
  rb_define_module_function(module, "cdf_lognormal_Q", rb_gsl_cdf_lognormal_Q, 3);
  rb_define_module_function(module, "cdf_lognormal_Pinv", rb_gsl_cdf_lognormal_Pinv, 3);
  rb_define_module_function(module, "cdf_lognormal_Qinv", rb_gsl_cdf_lognormal_Qinv, 3);

  rb_define_module_function(module, "cdf_chisq_P", rb_gsl_cdf_chisq_P, 2);
  rb_define_module_function(module, "cdf_chisq_Q", rb_gsl_cdf_chisq_Q, 2);
  rb_define_module_function(module, "cdf_chisq_Pinv", rb_gsl_cdf_chisq_Pinv, 2);
  rb_define_module_function(module, "cdf_chisq_Qinv", rb_gsl_cdf_chisq_Qinv, 2);
  
  rb_define_module_function(module, "cdf_fdist_P", rb_gsl_cdf_fdist_P, 3);
  rb_define_module_function(module, "cdf_fdist_Q", rb_gsl_cdf_fdist_Q, 3);

  rb_define_module_function(module, "cdf_tdist_P", rb_gsl_cdf_tdist_P, 2);
  rb_define_module_function(module, "cdf_tdist_Q", rb_gsl_cdf_tdist_Q, 2);
  rb_define_module_function(module, "cdf_tdist_Pinv", rb_gsl_cdf_tdist_Pinv, 2);
  rb_define_module_function(module, "cdf_tdist_Qinv", rb_gsl_cdf_tdist_Qinv, 2);
  
  rb_define_module_function(module, "cdf_beta_P", rb_gsl_cdf_beta_P, 3);
  rb_define_module_function(module, "cdf_beta_Q", rb_gsl_cdf_beta_Q, 3);

  rb_define_module_function(module, "cdf_logistic_P", rb_gsl_cdf_logistic_P, 2);
  rb_define_module_function(module, "cdf_logistic_Q", rb_gsl_cdf_logistic_Q, 2);
  rb_define_module_function(module, "cdf_logistic_Pinv", rb_gsl_cdf_logistic_Pinv, 2);
  rb_define_module_function(module, "cdf_logistic_Qinv", rb_gsl_cdf_logistic_Qinv, 2);
  
  rb_define_module_function(module, "cdf_pareto_P", rb_gsl_cdf_pareto_P, 3);
  rb_define_module_function(module, "cdf_pareto_Q", rb_gsl_cdf_pareto_Q, 3);
  rb_define_module_function(module, "cdf_pareto_Pinv", rb_gsl_cdf_pareto_Pinv, 3);
  rb_define_module_function(module, "cdf_pareto_Qinv", rb_gsl_cdf_pareto_Qinv, 3);

  rb_define_module_function(module, "cdf_weibull_P", rb_gsl_cdf_weibull_P, 3);
  rb_define_module_function(module, "cdf_weibull_Q", rb_gsl_cdf_weibull_Q, 3);
  rb_define_module_function(module, "cdf_weibull_Pinv", rb_gsl_cdf_weibull_Pinv, 3);
  rb_define_module_function(module, "cdf_weibull_Qinv", rb_gsl_cdf_weibull_Qinv, 3);

  rb_define_module_function(module, "cdf_gumbel1_P", rb_gsl_cdf_gumbel1_P, 3);
  rb_define_module_function(module, "cdf_gumbel1_Q", rb_gsl_cdf_gumbel1_Q, 3);
  rb_define_module_function(module, "cdf_gumbel1_Pinv", rb_gsl_cdf_gumbel1_Pinv, 3);
  rb_define_module_function(module, "cdf_gumbel1_Qinv", rb_gsl_cdf_gumbel1_Qinv, 3);

  rb_define_module_function(module, "cdf_gumbel2_P", rb_gsl_cdf_gumbel2_P, 3);
  rb_define_module_function(module, "cdf_gumbel2_Q", rb_gsl_cdf_gumbel2_Q, 3);
  rb_define_module_function(module, "cdf_gumbel2_Pinv", rb_gsl_cdf_gumbel2_Pinv, 3);
  rb_define_module_function(module, "cdf_gumbel2_Qinv", rb_gsl_cdf_gumbel2_Qinv, 3);

  /*****/

  rb_define_module_function(mgsl_cdf, "gaussian_P", rb_gsl_cdf_gaussian_P, -1);
  rb_define_module_function(mgsl_cdf, "gaussian_Q", rb_gsl_cdf_gaussian_Q, -1);
  rb_define_module_function(mgsl_cdf, "gaussian_Pinv", rb_gsl_cdf_gaussian_Pinv, -1);
  rb_define_module_function(mgsl_cdf, "gaussian_Qinv", rb_gsl_cdf_gaussian_Qinv, -1);
  rb_define_module_function(mgsl_cdf, "ugaussian_P", rb_gsl_cdf_gaussian_P, -1);
  rb_define_module_function(mgsl_cdf, "ugaussian_Q", rb_gsl_cdf_gaussian_Q, -1);
  rb_define_module_function(mgsl_cdf, "ugaussian_Pinv", rb_gsl_cdf_gaussian_Pinv, -1);
  rb_define_module_function(mgsl_cdf, "ugaussian_Qinv", rb_gsl_cdf_gaussian_Qinv, -1);

  rb_define_module_function(mgsl_cdf, "exponential_P", rb_gsl_cdf_exponential_P, 2);
 rb_define_module_function(mgsl_cdf, "exponential_Q", rb_gsl_cdf_exponential_Q, 2);
  rb_define_module_function(mgsl_cdf, "exponential_Pinv", rb_gsl_cdf_exponential_Pinv, 2);
  rb_define_module_function(mgsl_cdf, "exponential_Qinv", rb_gsl_cdf_exponential_Qinv, 2);

 rb_define_module_function(mgsl_cdf, "laplace_P", rb_gsl_cdf_laplace_P, 2);
 rb_define_module_function(mgsl_cdf, "laplace_Q", rb_gsl_cdf_laplace_Q, 2);
  rb_define_module_function(mgsl_cdf, "laplace_Pinv", rb_gsl_cdf_laplace_Pinv, 2);
  rb_define_module_function(mgsl_cdf, "laplace_Qinv", rb_gsl_cdf_laplace_Qinv, 2);

  rb_define_module_function(mgsl_cdf, "cauchy_P", rb_gsl_cdf_cauchy_P, 2);
  rb_define_module_function(mgsl_cdf, "cauchy_Q", rb_gsl_cdf_cauchy_Q, 2);
  rb_define_module_function(mgsl_cdf, "cauchy_Pinv", rb_gsl_cdf_cauchy_Pinv, 2);
  rb_define_module_function(mgsl_cdf, "cauchy_Qinv", rb_gsl_cdf_cauchy_Qinv, 2);

  rb_define_module_function(mgsl_cdf, "rayleigh_P", rb_gsl_cdf_rayleigh_P, 2);
  rb_define_module_function(mgsl_cdf, "rayleigh_Q", rb_gsl_cdf_rayleigh_Q, 2);
  rb_define_module_function(mgsl_cdf, "rayleigh_Pinv", rb_gsl_cdf_rayleigh_Pinv, 2);
  rb_define_module_function(mgsl_cdf, "rayleigh_Qinv", rb_gsl_cdf_rayleigh_Qinv, 2);

  rb_define_module_function(mgsl_cdf, "gamma_P", rb_gsl_cdf_gamma_P, 3);
  rb_define_module_function(mgsl_cdf, "gamma_Q", rb_gsl_cdf_gamma_Q, 3);
  rb_define_module_function(mgsl_cdf, "gamma_Pinv", rb_gsl_cdf_gamma_Pinv, 3);
  rb_define_module_function(mgsl_cdf, "gamma_Qinv", rb_gsl_cdf_gamma_Qinv, 3);

  rb_define_module_function(mgsl_cdf, "flat_P", rb_gsl_cdf_flat_P, 3);
  rb_define_module_function(mgsl_cdf, "flat_Q", rb_gsl_cdf_flat_Q, 3);
  rb_define_module_function(mgsl_cdf, "flat_Pinv", rb_gsl_cdf_flat_Pinv, 3);
  rb_define_module_function(mgsl_cdf, "flat_Qinv", rb_gsl_cdf_flat_Qinv, 3);

  rb_define_module_function(mgsl_cdf, "lognormal_P", rb_gsl_cdf_lognormal_P, 3);
  rb_define_module_function(mgsl_cdf, "lognormal_Q", rb_gsl_cdf_lognormal_Q, 3);
  rb_define_module_function(mgsl_cdf, "lognormal_Pinv", rb_gsl_cdf_lognormal_Pinv, 3);
  rb_define_module_function(mgsl_cdf, "lognormal_Qinv", rb_gsl_cdf_lognormal_Qinv, 3);

  rb_define_module_function(mgsl_cdf, "chisq_P", rb_gsl_cdf_chisq_P, 2);
  rb_define_module_function(mgsl_cdf, "chisq_Q", rb_gsl_cdf_chisq_Q, 2);
  rb_define_module_function(mgsl_cdf, "chisq_Pinv", rb_gsl_cdf_chisq_Pinv, 2);
  rb_define_module_function(mgsl_cdf, "chisq_Qinv", rb_gsl_cdf_chisq_Qinv, 2);
  
  rb_define_module_function(mgsl_cdf, "fdist_P", rb_gsl_cdf_fdist_P, 3);
  rb_define_module_function(mgsl_cdf, "fdist_Q", rb_gsl_cdf_fdist_Q, 3);

  rb_define_module_function(mgsl_cdf, "tdist_P", rb_gsl_cdf_tdist_P, 2);
  rb_define_module_function(mgsl_cdf, "tdist_Q", rb_gsl_cdf_tdist_Q, 2);
  rb_define_module_function(mgsl_cdf, "tdist_Pinv", rb_gsl_cdf_tdist_Pinv, 2);
  rb_define_module_function(mgsl_cdf, "tdist_Qinv", rb_gsl_cdf_tdist_Qinv, 2);
  
  rb_define_module_function(mgsl_cdf, "beta_P", rb_gsl_cdf_beta_P, 3);
  rb_define_module_function(mgsl_cdf, "beta_Q", rb_gsl_cdf_beta_Q, 3);

  rb_define_module_function(mgsl_cdf, "logistic_P", rb_gsl_cdf_logistic_P, 2);
  rb_define_module_function(mgsl_cdf, "logistic_Q", rb_gsl_cdf_logistic_Q, 2);
  rb_define_module_function(mgsl_cdf, "logistic_Pinv", rb_gsl_cdf_logistic_Pinv, 2);
  rb_define_module_function(mgsl_cdf, "logistic_Qinv", rb_gsl_cdf_logistic_Qinv, 2);
  
  rb_define_module_function(mgsl_cdf, "pareto_P", rb_gsl_cdf_pareto_P, 3);
  rb_define_module_function(mgsl_cdf, "pareto_Q", rb_gsl_cdf_pareto_Q, 3);
  rb_define_module_function(mgsl_cdf, "pareto_Pinv", rb_gsl_cdf_pareto_Pinv, 3);
  rb_define_module_function(mgsl_cdf, "pareto_Qinv", rb_gsl_cdf_pareto_Qinv, 3);

  rb_define_module_function(mgsl_cdf, "weibull_P", rb_gsl_cdf_weibull_P, 3);
  rb_define_module_function(mgsl_cdf, "weibull_Q", rb_gsl_cdf_weibull_Q, 3);
  rb_define_module_function(mgsl_cdf, "weibull_Pinv", rb_gsl_cdf_weibull_Pinv, 3);
  rb_define_module_function(mgsl_cdf, "weibull_Qinv", rb_gsl_cdf_weibull_Qinv, 3);

  rb_define_module_function(mgsl_cdf, "gumbel1_P", rb_gsl_cdf_gumbel1_P, 3);
  rb_define_module_function(mgsl_cdf, "gumbel1_Q", rb_gsl_cdf_gumbel1_Q, 3);
  rb_define_module_function(mgsl_cdf, "gumbel1_Pinv", rb_gsl_cdf_gumbel1_Pinv, 3);
  rb_define_module_function(mgsl_cdf, "gumbel1_Qinv", rb_gsl_cdf_gumbel1_Qinv, 3);

  rb_define_module_function(mgsl_cdf, "gumbel2_P", rb_gsl_cdf_gumbel2_P, 3);
  rb_define_module_function(mgsl_cdf, "gumbel2_Q", rb_gsl_cdf_gumbel2_Q, 3);
  rb_define_module_function(mgsl_cdf, "gumbel2_Pinv", rb_gsl_cdf_gumbel2_Pinv, 3);
  rb_define_module_function(mgsl_cdf, "gumbel2_Qinv", rb_gsl_cdf_gumbel2_Qinv, 3);


#ifdef GSL_1_6_LATER
  rb_define_module_function(module, "cdf_exppow_P", rb_gsl_cdf_exppow_P, 3);
  rb_define_module_function(module, "cdf_exppow_Q", rb_gsl_cdf_exppow_Q, 3);

  rb_define_module_function(mgsl_cdf, "exppow_P", rb_gsl_cdf_exppow_P, 3);
  rb_define_module_function(mgsl_cdf, "exppow_Q", rb_gsl_cdf_exppow_Q, 3);
#endif

#ifdef GSL_1_8_LATER
  rb_define_module_function(module, "cdf_binomial_P", rb_gsl_cdf_binomial_P, 3);
  rb_define_module_function(module, "cdf_binomial_Q", rb_gsl_cdf_binomial_Q, 3);
  rb_define_module_function(mgsl_cdf, "binomial_P", rb_gsl_cdf_binomial_P, 3);
  rb_define_module_function(mgsl_cdf, "binomial_Q", rb_gsl_cdf_binomial_Q, 3);

  rb_define_module_function(module, "cdf_poisson_P", rb_gsl_cdf_poisson_P, 2);
  rb_define_module_function(module, "cdf_poisson_Q", rb_gsl_cdf_poisson_Q, 2);
  rb_define_module_function(mgsl_cdf, "poisson_P", rb_gsl_cdf_poisson_P, 2);
  rb_define_module_function(mgsl_cdf, "poisson_Q", rb_gsl_cdf_poisson_Q, 2);

  rb_define_module_function(module, "cdf_geometric_P", rb_gsl_cdf_geometric_P, 2);
  rb_define_module_function(module, "cdf_geometric_Q", rb_gsl_cdf_geometric_Q, 2);
  rb_define_module_function(mgsl_cdf, "geometric_P", rb_gsl_cdf_geometric_P, 2);
  rb_define_module_function(mgsl_cdf, "geometric_Q", rb_gsl_cdf_geometric_Q, 2);

  rb_define_module_function(module, "cdf_negative_binomial_P", rb_gsl_cdf_negative_binomial_P, 3);
  rb_define_module_function(module, "cdf_negative_binomial_Q", rb_gsl_cdf_negative_binomial_Q, 3);
  rb_define_module_function(mgsl_cdf, "negative_binomial_P", rb_gsl_cdf_negative_binomial_P, 3);
  rb_define_module_function(mgsl_cdf, "negative_binomial_Q", rb_gsl_cdf_negative_binomial_Q, 3);

  rb_define_module_function(module, "cdf_pascal_P", rb_gsl_cdf_pascal_P, 3);
  rb_define_module_function(module, "cdf_pascal_Q", rb_gsl_cdf_pascal_Q, 3);
  rb_define_module_function(mgsl_cdf, "pascal_P", rb_gsl_cdf_pascal_P, 3);
  rb_define_module_function(mgsl_cdf, "pascal_Q", rb_gsl_cdf_pascal_Q, 3);

  rb_define_module_function(module, "cdf_hypergeometric_P", rb_gsl_cdf_hypergeometric_P, 4);
  rb_define_module_function(module, "cdf_hypergeometric_Q", rb_gsl_cdf_hypergeometric_Q, 4);
  rb_define_module_function(mgsl_cdf, "hypergeometric_P", rb_gsl_cdf_hypergeometric_P, 4);
  rb_define_module_function(mgsl_cdf, "hypergeometric_Q", rb_gsl_cdf_hypergeometric_Q, 4);

  rb_define_module_function(module, "cdf_beta_Pinv", rb_gsl_cdf_beta_Pinv, 3);
  rb_define_module_function(mgsl_cdf, "beta_Pinv", rb_gsl_cdf_beta_Pinv, 3);
  rb_define_module_function(module, "cdf_beta_Qinv", rb_gsl_cdf_beta_Qinv, 3);
  rb_define_module_function(mgsl_cdf, "beta_Qinv", rb_gsl_cdf_beta_Qinv, 3);

  rb_define_module_function(module, "cdf_fdist_Pinv", rb_gsl_cdf_fdist_Pinv, 3);
  rb_define_module_function(mgsl_cdf, "fdist_Pinv", rb_gsl_cdf_fdist_Pinv, 3);
  rb_define_module_function(module, "cdf_fdist_Qinv", rb_gsl_cdf_fdist_Qinv, 3);
  rb_define_module_function(mgsl_cdf, "fdist_Qinv", rb_gsl_cdf_fdist_Qinv, 3);
#endif

}
#endif
