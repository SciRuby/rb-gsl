/*
  randist.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_array.h"
#include "include/rb_gsl_common.h"
#include "include/rb_gsl_rng.h"
#include <gsl/gsl_randist.h>

VALUE rb_gsl_eval_pdf_cdf(VALUE xx, double (*f)(double));
VALUE rb_gsl_eval_pdf_cdf2(VALUE xx, VALUE aa, double (*f)(double, double));
VALUE rb_gsl_eval_pdf_cdf3(VALUE xx, VALUE aa, VALUE bb, 
			   double (*f)(double, double, double));
VALUE rb_gsl_eval_pdf_cdf2_uint(VALUE xx, VALUE aa, 
				double (*f)(unsigned int, double));

static VALUE rb_gsl_ran_eval0(int argc, VALUE *argv, VALUE obj,
			      double (*f)(const gsl_rng*))
{
  gsl_rng *r;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 1:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 0:
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0)", argc);
      break;
    }
    Data_Get_Struct(obj, gsl_rng, r);  
    break;
  }
  return rb_float_new((*f)(r));
}

static VALUE rb_gsl_ran_eval1(int argc, VALUE *argv, VALUE obj,
			      double (*f)(const gsl_rng*, double))
{
  gsl_rng *r;
  gsl_vector *v;
  size_t n, i;
  double a;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 3:
      n = NUM2INT(argv[2]);
      a = NUM2DBL(argv[1]);
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);    
      v = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) 
	gsl_vector_set(v, i, (*f)(r, a));
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    case 2:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      a = NUM2DBL(argv[1]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 2:
      n = NUM2INT(argv[1]);
      a = NUM2DBL(argv[0]);
      Data_Get_Struct(obj, gsl_rng, r);    
      v = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) 
	gsl_vector_set(v, i, (*f)(r, a));
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    case 1:
      a = NUM2DBL(argv[0]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
      break;
    }
    Data_Get_Struct(obj, gsl_rng, r);  
    break;
  }
  return rb_float_new((*f)(r, a));
}

static VALUE rb_gsl_ran_eval1_uint(int argc, VALUE *argv, VALUE obj,
			      unsigned int (*f)(const gsl_rng*, double))
{
  gsl_rng *r;
  gsl_vector_int *v;
  size_t n, i;
  double a;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 3:
      n = NUM2INT(argv[2]);
      a = NUM2DBL(argv[1]);
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);    
      v = gsl_vector_int_alloc(n);
      for (i = 0; i < n; i++) 
	gsl_vector_int_set(v, i, (int) (*f)(r, a));
      return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, v);
      break;
    case 2:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      a = NUM2DBL(argv[1]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 2:
      n = NUM2INT(argv[1]);
      a = NUM2DBL(argv[0]);
      Data_Get_Struct(obj, gsl_rng, r);    
      v = gsl_vector_int_alloc(n);
      for (i = 0; i < n; i++) 
	gsl_vector_int_set(v, i, (*f)(r, a));
      return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, v);
      break;
    case 1:
      a = NUM2DBL(argv[0]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
      break;
    }
    Data_Get_Struct(obj, gsl_rng, r);  
    break;
  }
  return INT2FIX((*f)(r, a));
}

static VALUE rb_gsl_ran_eval2(int argc, VALUE *argv, VALUE obj,
			      double (*f)(const gsl_rng*, double, double))
{
  gsl_rng *r;
  gsl_vector *v;
  size_t n, i;
  double a, b;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 4:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      a = NUM2DBL(argv[1]); b = NUM2DBL(argv[2]);
      n = NUM2INT(argv[3]);
      v = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) gsl_vector_set(v, i, (*f)(r, a, b));
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    case 3:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      a = NUM2DBL(argv[1]); b = NUM2DBL(argv[2]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 3)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 3:
      Data_Get_Struct(obj, gsl_rng, r);  
      a = NUM2DBL(argv[0]); b = NUM2DBL(argv[1]);
      n = NUM2INT(argv[2]);
      v = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) gsl_vector_set(v, i, (*f)(r, a, b));
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    case 2:
      a = NUM2DBL(argv[0]); b = NUM2DBL(argv[1]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)", argc);
      break;
    }
    Data_Get_Struct(obj, gsl_rng, r);  
    break;
  }
  return rb_float_new((*f)(r, a, b));
}

static VALUE rb_gsl_ran_eval3(int argc, VALUE *argv, VALUE obj,
			      double (*f)(const gsl_rng*, double, double, double))
{
  gsl_rng *r;
  gsl_vector *v;
  size_t n, i;
  double a, b, c;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 5:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      a = NUM2DBL(argv[1]); b = NUM2DBL(argv[2]); c = NUM2DBL(argv[3]);
      n = NUM2INT(argv[4]);
      v = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) gsl_vector_set(v, i, (*f)(r, a, b, c));
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    case 4:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      a = NUM2DBL(argv[1]); b = NUM2DBL(argv[2]); c = NUM2DBL(argv[3]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 4)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 4:
      Data_Get_Struct(obj, gsl_rng, r);  
      a = NUM2DBL(argv[0]); b = NUM2DBL(argv[1]); c = NUM2DBL(argv[2]);
      n = NUM2INT(argv[3]);
      v = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) gsl_vector_set(v, i, (*f)(r, a, b, c));
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    case 3:
      a = NUM2DBL(argv[0]); b = NUM2DBL(argv[1]); c = NUM2DBL(argv[2]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 3)", argc);
      break;
    }
    Data_Get_Struct(obj, gsl_rng, r);  
    break;
  }
  return rb_float_new((*f)(r, a, b, c));
}

static VALUE rb_gsl_ran_gaussian(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  gsl_vector *v;
  size_t n, i;
  double sigma = 1.0;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 3:
      n = NUM2INT(argv[2]);
      sigma = NUM2DBL(argv[1]);
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);    
      v = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) 
	gsl_vector_set(v, i, gsl_ran_gaussian(r, sigma));
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    case 2:
      sigma = NUM2DBL(argv[1]);
      /* no break */
    case 1:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);    
      return rb_float_new(gsl_ran_gaussian(r, sigma));
      break;
    default:
      rb_raise(rb_eArgError, "too many arguments (%d for 1 or 2)", argc);
      return Qnil;
      break;
    }
    break;
  default:
    Data_Get_Struct(obj, gsl_rng, r);    
    switch (argc) {
    case 2:
      n = NUM2INT(argv[1]);
      sigma = NUM2DBL(argv[0]);
      v = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) 
	gsl_vector_set(v, i, gsl_ran_gaussian(r, sigma));
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    case 1:
      sigma = NUM2DBL(argv[0]);
      /* no break */
    case 0:
      return rb_float_new(gsl_ran_gaussian(r, sigma));
      break;
    default:
      rb_raise(rb_eArgError, "too many arguments (%d for 0 or 1)", argc);
      return Qnil;
      break;
    }
    break;
  }
  return Qnil;
}

static VALUE rb_gsl_ran_gaussian_ratio_method(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  double sigma = 1.0;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 2:
      sigma = NUM2DBL(argv[1]);
      /* no break */
    case 1:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);    
      return rb_float_new(gsl_ran_gaussian_ratio_method(r, sigma));
      break;
    default:
      rb_raise(rb_eArgError, "too many arguments (%d for 1 or 2)", argc);
      break;
    }
    break;
  default:
    Data_Get_Struct(obj, gsl_rng, r);    
    switch (argc) {
    case 1:
      sigma = NUM2DBL(argv[0]);
      /* no break */
    case 0:
      return rb_float_new(gsl_ran_gaussian_ratio_method(r, sigma));
      break;
    default:
      rb_raise(rb_eArgError, "too many arguments (%d for 0 or 1)", argc);
      break;
    }
    break;
  }
  return Qnil;
}

static VALUE rb_gsl_ran_gaussian_pdf(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) {
    return rb_gsl_eval_pdf_cdf(argv[0], gsl_ran_ugaussian_pdf);
  } else if (argc == 2) {
    return rb_gsl_eval_pdf_cdf2(argv[0], argv[1], gsl_ran_gaussian_pdf);
  } else {
    rb_raise(rb_eArgError, "wrong number of arguments (1 or 2)");
  }
}

static VALUE rb_gsl_ran_gaussian_tail(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  gsl_vector *v;
  size_t n, i;
  double a, sigma = 1.0;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 4:
      n = NUM2INT(argv[3]);
      sigma = NUM2DBL(argv[2]);
      a = NUM2DBL(argv[1]);
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);    
      v = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) 
	gsl_vector_set(v, i, gsl_ran_gaussian_tail(r, a, sigma));
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    case 3:
      sigma = NUM2DBL(argv[2]);
      /* no break */
    case 2:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);    
      a = NUM2DBL(argv[1]);
      return rb_float_new(gsl_ran_gaussian_tail(r, a, sigma));
      break;
    default:
      rb_raise(rb_eArgError, "too many arguments (%d for 2 or 3)", argc);
      break;
    }
    break;
  default:
    Data_Get_Struct(obj, gsl_rng, r);    
    switch (argc) {
    case 3:
      n = NUM2INT(argv[2]);
      sigma = NUM2DBL(argv[1]);
      a = NUM2DBL(argv[0]);
      v = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) 
	gsl_vector_set(v, i, gsl_ran_gaussian_tail(r, a, sigma));
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    case 2:
      sigma = NUM2DBL(argv[1]);
      /* no break */
    case 1:
      a = NUM2DBL(argv[0]);
      return rb_float_new(gsl_ran_gaussian_tail(r, a, sigma));
      break;
    default:
      rb_raise(rb_eArgError, "too many arguments (%d for 1 or 2)", argc);
      break;
    }
    break;
  }
}

static VALUE rb_gsl_ran_gaussian_tail_pdf(int argc, VALUE *argv, VALUE obj)
{
  switch (argc) {
  case 2:
    return rb_gsl_eval_pdf_cdf2(argv[0], argv[1], gsl_ran_ugaussian_tail_pdf);
    break;
  case 3:
    return rb_gsl_eval_pdf_cdf3(argv[0], argv[1], argv[2], gsl_ran_gaussian_tail_pdf);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
}

static VALUE rb_gsl_ran_bivariate_gaussian(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  double sigmax, sigmay, x, y, rho;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 4:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r); 
      sigmax = NUM2DBL(argv[1]); sigmay = NUM2DBL(argv[2]);
      rho = NUM2DBL(argv[3]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 4)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 3:
      Data_Get_Struct(obj, gsl_rng, r);  
      sigmax = NUM2DBL(argv[0]); sigmay = NUM2DBL(argv[1]);
      rho = NUM2DBL(argv[2]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 3)", argc);
      break;
    }
    break;
  }
  gsl_ran_bivariate_gaussian(r, sigmax, sigmay, rho, &x, &y);
  return rb_ary_new3(2, rb_float_new(x), rb_float_new(y));
}

static VALUE rb_gsl_ran_bivariate_gaussian_pdf(VALUE obj, VALUE x, VALUE y, VALUE sx, VALUE sy, VALUE rr)
{
  Need_Float(x); Need_Float(y);
  Need_Float(sx); Need_Float(sy); Need_Float(rr);
  return rb_float_new(gsl_ran_bivariate_gaussian_pdf(NUM2DBL(x), NUM2DBL(y), NUM2DBL(sx), NUM2DBL(sy), NUM2DBL(rr)));
}

static VALUE rb_gsl_ran_exponential(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval1(argc, argv, obj, gsl_ran_exponential);
}

static VALUE rb_gsl_ran_exponential_pdf(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_ran_exponential_pdf);
}

static VALUE rb_gsl_ran_laplace(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval1(argc, argv, obj, gsl_ran_laplace);
}

static VALUE rb_gsl_ran_laplace_pdf(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_ran_laplace_pdf);
}

static VALUE rb_gsl_ran_exppow(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_exppow);
}

static VALUE rb_gsl_ran_exppow_pdf(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_ran_exppow_pdf);
}

static VALUE rb_gsl_ran_cauchy(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval1(argc, argv, obj, gsl_ran_cauchy);
}

static VALUE rb_gsl_ran_cauchy_pdf(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_ran_cauchy_pdf);
}

static VALUE rb_gsl_ran_rayleigh(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval1(argc, argv, obj, gsl_ran_rayleigh);
}

static VALUE rb_gsl_ran_rayleigh_pdf(VALUE obj, VALUE x, VALUE a)
{
  return rb_gsl_eval_pdf_cdf2(x, a, gsl_ran_rayleigh_pdf);
}

static VALUE rb_gsl_ran_rayleigh_tail(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_rayleigh_tail);
}

static VALUE rb_gsl_ran_rayleigh_tail_pdf(VALUE obj, VALUE x, VALUE a, VALUE s)
{
  return rb_gsl_eval_pdf_cdf3(x, a, s, gsl_ran_rayleigh_tail_pdf);
}

static VALUE rb_gsl_ran_landau(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval0(argc, argv, obj, gsl_ran_landau);
}

static VALUE rb_gsl_ran_landau_pdf(VALUE obj, VALUE x)
{
  return rb_gsl_eval_pdf_cdf(x, gsl_ran_landau_pdf);
}

static VALUE rb_gsl_ran_levy(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_levy);
}

static VALUE rb_gsl_ran_levy_skew(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval3(argc, argv, obj, gsl_ran_levy_skew);
}

static VALUE rb_gsl_ran_gamma(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_gamma);
}

static VALUE rb_gsl_ran_gamma_pdf(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_ran_gamma_pdf);
}

static VALUE rb_gsl_ran_flat(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_flat);
}

static VALUE rb_gsl_ran_flat_pdf(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_ran_flat_pdf);
}

static VALUE rb_gsl_ran_lognormal(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_lognormal);
}

static VALUE rb_gsl_ran_lognormal_pdf(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_ran_lognormal_pdf);
}

static VALUE rb_gsl_ran_chisq(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval1(argc, argv, obj, gsl_ran_chisq);
}

static VALUE rb_gsl_ran_chisq_pdf(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_ran_chisq_pdf);
}

static VALUE rb_gsl_ran_fdist(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_fdist);
}

static VALUE rb_gsl_ran_fdist_pdf(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_ran_fdist_pdf);
}

static VALUE rb_gsl_ran_tdist(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval1(argc, argv, obj, gsl_ran_tdist);
}

static VALUE rb_gsl_ran_tdist_pdf(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_ran_tdist_pdf);
}

static VALUE rb_gsl_ran_beta(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_beta);
}

static VALUE rb_gsl_ran_beta_pdf(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_ran_beta_pdf);
}

static VALUE rb_gsl_ran_logistic(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval1(argc, argv, obj, gsl_ran_logistic);
}

static VALUE rb_gsl_ran_logistic_pdf(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2(x, mu, gsl_ran_logistic_pdf);
}

static VALUE rb_gsl_ran_pareto(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_pareto);
}

static VALUE rb_gsl_ran_pareto_pdf(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_ran_pareto_pdf);
}

static VALUE rb_gsl_ran_weibull(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_weibull);
}

static VALUE rb_gsl_ran_weibull_pdf(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_ran_weibull_pdf);
}

static VALUE rb_gsl_ran_gumbel1(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_gumbel1);
}

static VALUE rb_gsl_ran_gumbel1_pdf(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_ran_gumbel1_pdf);
}

static VALUE rb_gsl_ran_gumbel2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_gumbel2);
}

static VALUE rb_gsl_ran_gumbel2_pdf(VALUE obj, VALUE x, VALUE a, VALUE b)
{
  return rb_gsl_eval_pdf_cdf3(x, a, b, gsl_ran_gumbel2_pdf);
}

static VALUE rb_gsl_ran_poisson(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval1_uint(argc, argv, obj, gsl_ran_poisson);
}

static VALUE rb_gsl_ran_poisson_pdf(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2_uint(x, mu, gsl_ran_poisson_pdf);
}

static VALUE rb_gsl_ran_bernoulli(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval1_uint(argc, argv, obj, gsl_ran_bernoulli);
}

static VALUE rb_gsl_ran_bernoulli_pdf(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2_uint(x, mu, gsl_ran_bernoulli_pdf);
}

static VALUE rb_gsl_ran_binomial(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r;
  double p;
  unsigned int n;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 3:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      p = NUM2DBL(argv[1]);
      n = FIX2UINT(argv[2]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 3)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 2:
      p = NUM2DBL(argv[0]);
      n = FIX2UINT(argv[1]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)", argc);
      break;
    }
    Data_Get_Struct(obj, gsl_rng, r);  
    break;
  }
  return UINT2NUM(gsl_ran_binomial(r, p, n));
}

#ifdef GSL_1_4_LATER
static VALUE rb_gsl_ran_binomial_tpe(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  double p;
  unsigned int n;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 3:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      p = NUM2DBL(argv[1]);
      n = FIX2UINT(argv[2]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 3)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 2:
      p = NUM2DBL(argv[0]);
      n = FIX2UINT(argv[1]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)", argc);
      break;
    }
    Data_Get_Struct(obj, gsl_rng, r);  
    break;
  }
  return UINT2NUM(gsl_ran_binomial_tpe(r, p, n));
}
#endif

static VALUE rb_gsl_ran_binomial_pdf(VALUE obj, VALUE x, VALUE p, VALUE n)
{
  return rb_float_new(gsl_ran_binomial_pdf(NUM2UINT(x), NUM2DBL(p), NUM2UINT(n)));
}

static VALUE rb_gsl_ran_negative_binomial(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r;
  double p;
  unsigned int n;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 3:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      p = NUM2DBL(argv[1]);
      n = FIX2UINT(argv[2]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 3)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 2:
      p = NUM2DBL(argv[0]);
      n = FIX2UINT(argv[1]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)", argc);
      break;
    }
    Data_Get_Struct(obj, gsl_rng, r);  
    break;
  }
  return UINT2NUM(gsl_ran_negative_binomial(r, p, n));
}

static VALUE rb_gsl_ran_negative_binomial_pdf(VALUE obj, VALUE x, VALUE p, VALUE n)
{
  return rb_float_new(gsl_ran_negative_binomial_pdf(NUM2UINT(x), NUM2DBL(p), NUM2DBL(n)));
}

static VALUE rb_gsl_ran_pascal(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  double p;
  unsigned int n;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 3:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      p = NUM2DBL(argv[1]);
      n = FIX2UINT(argv[2]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 3)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 2:
      p = NUM2DBL(argv[0]);
      n = FIX2UINT(argv[1]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)", argc);
      break;
    }
    Data_Get_Struct(obj, gsl_rng, r);  
    break;
  }
  return UINT2NUM(gsl_ran_pascal(r, p, n));
}

static VALUE rb_gsl_ran_pascal_pdf(VALUE obj, VALUE x, VALUE p, VALUE n)
{
  return rb_float_new(gsl_ran_pascal_pdf(NUM2UINT(x), NUM2DBL(p), NUM2UINT(n)));
}

static VALUE rb_gsl_ran_geometric(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval1_uint(argc, argv, obj, gsl_ran_geometric);
}

static VALUE rb_gsl_ran_geometric_pdf(VALUE obj, VALUE x, VALUE p)
{
  return rb_gsl_eval_pdf_cdf2_uint(x, p, gsl_ran_geometric_pdf);
}

static VALUE rb_gsl_ran_hypergeometric(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r;
  unsigned int n1, n2, t;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 4:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      n1 = FIX2UINT(argv[1]);
      n2 = FIX2UINT(argv[2]);
      t = FIX2UINT(argv[3]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 4)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 3:
      n1 = FIX2UINT(argv[0]);
      n2 = FIX2UINT(argv[1]);
      t = FIX2UINT(argv[2]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 3)", argc);
      break;
    }
    Data_Get_Struct(obj, gsl_rng, r);  
    break;
  }
  return UINT2NUM(gsl_ran_hypergeometric(r, n1, n2, t));
}

static VALUE rb_gsl_ran_hypergeometric_pdf(VALUE obj, VALUE k, VALUE n1, VALUE n2, VALUE t)
{
  return rb_float_new(gsl_ran_hypergeometric_pdf(NUM2UINT(k), NUM2UINT(n1),  NUM2UINT(n2),  NUM2UINT(t)));
}

static VALUE rb_gsl_ran_logarithmic(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval1_uint(argc, argv, obj, gsl_ran_logarithmic);
}

static VALUE rb_gsl_ran_logarithmic_pdf(VALUE obj, VALUE x, VALUE mu)
{
  return rb_gsl_eval_pdf_cdf2_uint(x, mu, gsl_ran_logarithmic_pdf);
}

static VALUE rb_gsl_ran_dir_2d(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  double x, y;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 1:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 0:
      Data_Get_Struct(obj, gsl_rng, r);  
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0)", argc);
      break;
    }
  }
  gsl_ran_dir_2d(r, &x, &y);
  return rb_ary_new3(2, rb_float_new(x), rb_float_new(y));
}

static VALUE rb_gsl_ran_dir_2d_trig_method(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  double x, y;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 1:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 0:
      Data_Get_Struct(obj, gsl_rng, r);  
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0)", argc);
      break;
    }
  }
  gsl_ran_dir_2d_trig_method(r, &x, &y);
  return rb_ary_new3(2, rb_float_new(x), rb_float_new(y));
}

static VALUE rb_gsl_ran_dir_3d(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  double x, y, z;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 1:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 0:
      Data_Get_Struct(obj, gsl_rng, r);  
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0)", argc);
      break;
    }
  }
  gsl_ran_dir_3d(r, &x, &y, &z);
  return rb_ary_new3(3, rb_float_new(x), rb_float_new(y), rb_float_new(z));
}

static VALUE rb_gsl_ran_dir_nd(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  size_t n;
  gsl_vector *v;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 2:
      CHECK_RNG(argv[0]); CHECK_FIXNUM(argv[1]);
      Data_Get_Struct(argv[0], gsl_rng, r);  
      n = FIX2INT(argv[1]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 1:
      CHECK_FIXNUM(argv[0]);
      n = FIX2INT(argv[0]);
      Data_Get_Struct(obj, gsl_rng, r);  
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0)", argc);
      break;
    }
  }
  v = gsl_vector_alloc(n);
  gsl_ran_dir_nd(r, n, v->data);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

static VALUE rb_gsl_ran_shuffle(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  gsl_vector *v = NULL;
  gsl_permutation *p = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    switch (argc) {
    case 2:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);
      if (VECTOR_P(argv[1])) {
	Data_Get_Struct(argv[1], gsl_vector, v);
	gsl_ran_shuffle(r, v->data, v->size, sizeof(double));
      } else if (PERMUTATION_P(argv[1])) {
	Data_Get_Struct(argv[1], gsl_permutation, p);
	gsl_ran_shuffle(r, p->data, p->size, sizeof(size_t));
      } else {
	rb_raise(rb_eTypeError, "wrong argument type %s (Vector or Permutation expected)", rb_class2name(CLASS_OF(argv[1])));
      }
      break;
    case 3:
      CHECK_RNG(argv[0]);
      CHECK_FIXNUM(argv[2]);
      Data_Get_Struct(argv[0], gsl_rng, r);
      if (VECTOR_P(argv[1])) {
	Data_Get_Struct(argv[1], gsl_vector, v);
	gsl_ran_shuffle(r, v->data, FIX2INT(argv[2]), sizeof(double));
      } else if (PERMUTATION_P(argv[1])) {
	Data_Get_Struct(argv[1], gsl_permutation, p);
	gsl_ran_shuffle(r, p->data, FIX2INT(argv[2]), sizeof(size_t));
      } else {
	rb_raise(rb_eTypeError, 
		 "wrong argument type %s (Vector or Permutation expected)", 
		 rb_class2name(CLASS_OF(argv[1])));
      }
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
      break;
    }
    break;
  default:
    Data_Get_Struct(obj, gsl_rng, r);  
    switch (argc) {
    case 1:
      if (VECTOR_P(argv[0])) {
	Data_Get_Struct(argv[0], gsl_vector, v);
	gsl_ran_shuffle(r, v->data, v->size, sizeof(double));
      } else if (PERMUTATION_P(argv[0])) {
	Data_Get_Struct(argv[0], gsl_permutation, p);
	gsl_ran_shuffle(r, p->data, p->size, sizeof(size_t));
      } else {
	rb_raise(rb_eTypeError, "wrong argument type %s (Vector or Permutation expected)", rb_class2name(CLASS_OF(argv[0])));
      }
      break;
    case 2:
      CHECK_FIXNUM(argv[1]);
      if (VECTOR_P(argv[0])) {
	Data_Get_Struct(argv[0], gsl_vector, v);
	gsl_ran_shuffle(r, v->data, FIX2INT(argv[1]), sizeof(double));
      } else if (PERMUTATION_P(argv[0])) {
	Data_Get_Struct(argv[0], gsl_permutation, p);
	gsl_ran_shuffle(r, p->data, FIX2INT(argv[1]), sizeof(size_t));
      } else {
	rb_raise(rb_eTypeError, 
		 "wrong argument type %s (Vector or Permutation expected)", 
		 rb_class2name(CLASS_OF(argv[0])));
      }
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
      break;
    }
    break;
  }
  return Qnil;
}

static VALUE rb_gsl_ran_choose(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  gsl_vector *v = NULL, *v2 = NULL;
  size_t k, n;
  Data_Get_Struct(obj, gsl_rng, r);  
  switch (argc) {
  case 2:
    CHECK_VECTOR(argv[0]);
    CHECK_FIXNUM(argv[1]);
    Data_Get_Struct(argv[0], gsl_vector, v);  
    n = v->size;
    k = FIX2INT(argv[1]);
    if (k > n) rb_raise(rb_eArgError, "the argument 1 must be less than or equal to the size of the vector.");
    v2 = gsl_vector_alloc(k);
    gsl_ran_choose(r, v2->data, k, v->data, n, sizeof(double));
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v2);;
    break;
  case 1:
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, v);  
    n = v->size;
    k = v->size;
    if (k > n) rb_raise(rb_eArgError, "the argument 1 must be less than or equal to the size of the vector.");
    v2 = gsl_vector_alloc(k);
    gsl_ran_choose(r, v2->data, k, v->data, n, sizeof(double));
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v2);;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
}

static VALUE rb_gsl_ran_choose_singleton(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  gsl_vector *v = NULL, *v2 = NULL;
  size_t k, n;
  switch (argc) {
  case 3:
    CHECK_RNG(argv[0]);
    CHECK_VECTOR(argv[1]);
    CHECK_FIXNUM(argv[2]);
    Data_Get_Struct(argv[0], gsl_rng, r);  
    Data_Get_Struct(argv[1], gsl_vector, v);  
    n = v->size;
    k = FIX2INT(argv[2]);
    if (k > n) rb_raise(rb_eArgError, "the argument 1 must be less than or equal to the size of the vector.");
    v2 = gsl_vector_alloc(k);
    gsl_ran_choose(r, v2->data, k, v->data, n, sizeof(double));
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v2);;
    break;
  case 2:
    CHECK_RNG(argv[0]);
    CHECK_VECTOR(argv[1]);
    Data_Get_Struct(argv[0], gsl_rng, r);  
    Data_Get_Struct(argv[1], gsl_vector, v);  
    n = v->size;
    k = v->size;
    if (k > n) rb_raise(rb_eArgError, "the argument 1 must be less than or equal to the size of the vector.");
    v2 = gsl_vector_alloc(k);
    gsl_ran_choose(r, v2->data, k, v->data, n, sizeof(double));
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v2);;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
}

static VALUE rb_gsl_ran_sample(VALUE obj, VALUE vv, VALUE kk)
{
  gsl_rng *r = NULL;	
  gsl_vector *v = NULL, *v2 = NULL;
  size_t k, n;
  Data_Get_Struct(obj, gsl_rng, r);  
  Data_Get_Struct(vv, gsl_vector, v);  
  n = v->size;
  k = FIX2INT(kk);
  v2 = gsl_vector_alloc(k);
  gsl_ran_sample(r, v2->data, k, v->data, n, sizeof(double));
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v2);;
}

#ifdef GSL_1_3_LATER
static VALUE rb_gsl_ran_dirichlet(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  gsl_vector *v = NULL, *v2 = NULL;
  Data_Get_Struct(obj, gsl_rng, r);    
  if (argc == 1) {
    Data_Get_Struct(argv[0], gsl_vector, v);    
    v2 = gsl_vector_alloc(v->size);
    gsl_ran_dirichlet(r, v->size, v->data, v2->data);
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v2);
  } else if (argc == 2) {
    Data_Get_Struct(argv[0], gsl_vector, v);    
    Data_Get_Struct(argv[1], gsl_vector, v2);    
    gsl_ran_dirichlet(r, v->size, v->data, v2->data);
    return obj;
  } else {
    rb_raise(rb_eArgError, "wrong number of arguments (1 or 2 GSL_Vectors)");
  }
}

static VALUE rb_gsl_ran_dirichlet_pdf(VALUE obj, VALUE a, VALUE t)
{
  gsl_vector *alpha = NULL, *theta = NULL;
  Data_Get_Struct(a, gsl_vector, alpha);    
  Data_Get_Struct(t, gsl_vector, theta);    
  return rb_float_new(gsl_ran_dirichlet_pdf(alpha->size, alpha->data, theta->data));
}
static VALUE rb_gsl_ran_dirichlet_lnpdf(VALUE obj, VALUE a, VALUE t)
{
  gsl_vector *alpha = NULL, *theta = NULL;
  Data_Get_Struct(a, gsl_vector, alpha);    
  Data_Get_Struct(t, gsl_vector, theta);    
  return rb_float_new(gsl_ran_dirichlet_lnpdf(alpha->size, alpha->data, theta->data));
}
#endif

static VALUE rb_gsl_ran_discrete_new(VALUE klass, VALUE vv)
{
  gsl_vector *v = NULL;
  gsl_ran_discrete_t *g = NULL;
  Data_Get_Struct(vv, gsl_vector, v);
  g = gsl_ran_discrete_preproc(v->size, v->data);
  return Data_Wrap_Struct(klass, 0, gsl_ran_discrete_free, g);
}

static VALUE rb_gsl_ran_discrete(VALUE obj, VALUE gg)
{
  gsl_rng *r = NULL;
  gsl_ran_discrete_t *g = NULL;
  Data_Get_Struct(obj, gsl_rng, r);
  Data_Get_Struct(gg, gsl_ran_discrete_t, g);
  return INT2FIX(gsl_ran_discrete(r, g));
}

static VALUE rb_gsl_ran_discrete_pdf(VALUE obj, VALUE k, VALUE gg)
{
  gsl_ran_discrete_t *g = NULL;
  Data_Get_Struct(gg, gsl_ran_discrete_t, g);
  return rb_float_new(gsl_ran_discrete_pdf(FIX2INT(k), g));
}

#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif

/*****/
VALUE rb_gsl_eval_pdf_cdf(VALUE xx, double (*f)(double))
{
  VALUE x, ary;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  size_t i, j, n;  
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2;
#endif
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch(TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new((*f)(NUM2DBL(xx)));
    break;
  case T_ARRAY:
    //    n = RARRAY(xx)->len;
    n = RARRAY_LEN(xx);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(x);
      rb_ary_store(ary, i, rb_float_new((*f)(NUM2DBL(x))));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      xx = na_change_type(xx, NA_DFLOAT);
      GetNArray(xx, na);
      ptr1 = (double *) na->ptr;
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*f)(ptr1[i]);
      return ary;
    }
#endif
    if (VECTOR_P(xx)) {
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	gsl_vector_set(vnew, i, (*f)(gsl_vector_get(v, i)));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MATRIX_P(xx)) {
      Data_Get_Struct(xx, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  gsl_matrix_set(mnew, i, j, (*f)(gsl_matrix_get(m, i, j)));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  }
  /* never reach here */
  return Qnil;
}

VALUE rb_gsl_eval_pdf_cdf2(VALUE xx, VALUE aa, 
			      double (*f)(double, double))
{
  VALUE x, ary;
  double a;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  size_t i, j, n;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2;
#endif
  Need_Float(aa);
  a = NUM2DBL(aa);
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch(TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new((*f)(NUM2DBL(xx), a));
    break;
  case T_ARRAY:
    //    n = RARRAY(xx)->len;
    n = RARRAY_LEN(xx);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(x);
      rb_ary_store(ary, i, rb_float_new((*f)(NUM2DBL(x), a)));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      xx = na_change_type(xx, NA_DFLOAT);
      GetNArray(xx, na);
      ptr1 = (double *) na->ptr;
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*f)(ptr1[i], a);
      return ary;
    }
#endif
    if (VECTOR_P(xx)) {
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	gsl_vector_set(vnew, i, (*f)(gsl_vector_get(v, i), a));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MATRIX_P(xx)) {
      Data_Get_Struct(xx, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  gsl_matrix_set(mnew, i, j, (*f)(gsl_matrix_get(m, i, j), a));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  }
  /* never reach here */
  return Qnil;
}

VALUE rb_gsl_eval_pdf_cdf3(VALUE xx, VALUE aa, VALUE bb, 
			      double (*f)(double, double, double))
{
  VALUE x, ary;
  double a, b;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  size_t i, j, n;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2;
#endif
  Need_Float(aa);  Need_Float(bb);
  a = NUM2DBL(aa); b = NUM2DBL(bb);
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch(TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new((*f)(NUM2DBL(xx), a, b));
    break;
  case T_ARRAY:
    //    n = RARRAY(xx)->len;
    n = RARRAY_LEN(xx);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(x);
      rb_ary_store(ary, i, rb_float_new((*f)(NUM2DBL(x), a, b)));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      xx = na_change_type(xx, NA_DFLOAT);
      GetNArray(xx, na);
      ptr1 = (double *) na->ptr;
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*f)(ptr1[i], a, b);
      return ary;
    }
#endif
    if (VECTOR_P(xx)) {
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	gsl_vector_set(vnew, i, (*f)(gsl_vector_get(v, i), a, b));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MATRIX_P(xx)) {
      Data_Get_Struct(xx, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  gsl_matrix_set(mnew, i, j, (*f)(gsl_matrix_get(m, i, j), a, b));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  }
  /* never reach here */
  return Qnil;
}

VALUE rb_gsl_eval_pdf_cdf2_uint(VALUE xx, VALUE aa, 
				       double (*f)(unsigned int, double))
{
  VALUE x, ary;
  double a;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_vector_int *vi = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  gsl_matrix_int *mi = NULL;
  size_t i, j, n;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  char *ptr1, *ptr2;
#endif
  Need_Float(aa);
  a = NUM2DBL(aa);
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch(TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new((*f)(NUM2UINT(xx), a));
    break;
  case T_ARRAY:
    //    n = RARRAY(xx)->len;
    n = RARRAY_LEN(xx);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      x = rb_ary_entry(xx, i);
      rb_ary_store(ary, i, rb_float_new((*f)(NUM2UINT(x), a)));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      GetNArray(xx, na);
      n = na->total;
      ary = na_make_object(na->type, na->rank, na->shape, CLASS_OF(xx));
      ptr1 = (char *)na->ptr;
      ptr2 = (char *)NA_STRUCT(ary)->ptr;
      switch (na->type) {
      case NA_LINT:
	for (i = 0; i < n; i++) 
	  ((int*)ptr2)[i] = (*f)((unsigned int) ((int*)ptr1)[i], a);
	break;
      default:
	for (i = 0; i < n; i++) 
	  ((double*)ptr2)[i] = (*f)((unsigned int) ((double*)ptr1)[i], a);
	break;
      }
      return ary;
    }
#endif
    if (VECTOR_P(xx)) {
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	gsl_vector_set(vnew, i, (*f)((unsigned int) gsl_vector_get(v, i), a));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (VECTOR_INT_P(xx)) {
      Data_Get_Struct(xx, gsl_vector_int, vi);
      vnew = gsl_vector_alloc(vi->size);
      for (i = 0; i < vi->size; i++) {
	gsl_vector_set(vnew, i, (*f)((unsigned int) gsl_vector_int_get(vi, i), a));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MATRIX_P(xx)) {
      Data_Get_Struct(xx, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  gsl_matrix_set(mnew, i, j, (*f)((unsigned int) gsl_matrix_get(m, i, j), a));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else if (MATRIX_INT_P(xx)) {
      Data_Get_Struct(xx, gsl_matrix_int, mi);
      mnew = gsl_matrix_alloc(mi->size1, mi->size2);
      for (i = 0; i < mi->size1; i++) {
	for (j = 0; j < mi->size2; j++) {
	  gsl_matrix_set(mnew, i, j, (*f)((unsigned int) gsl_matrix_int_get(mi, i, j), a));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  }
  /* never reach here */
  return Qnil;
}

/*
static VALUE rb_gsl_ran_erlang(VALUE obj, VALUE a, VALUE n)
{
  gsl_rng *r = NULL;
  Data_Get_Struct(obj, gsl_rng, r);
  return rb_float_new(gsl_ran_erlang(r, NUM2DBL(a), NUM2DBL(n)));
}

static VALUE rb_gsl_ran_erlang_pdf(VALUE obj, VALUE x, VALUE a, VALUE n)
{
  return rb_float_new(gsl_ran_erlang_pdf(NUM2DBL(x), NUM2DBL(a), NUM2DBL(n)));
}
*/

#ifdef GSL_1_8_LATER

static VALUE rb_gsl_ran_gaussian_ziggurat(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;	
  gsl_vector *v;
  size_t n, i;
  double sigma = 1.0;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 3:
      n = NUM2INT(argv[2]);
      sigma = NUM2DBL(argv[1]);
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);    
      v = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) 
	gsl_vector_set(v, i, gsl_ran_gaussian_ziggurat(r, sigma));
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    case 2:
      sigma = NUM2DBL(argv[1]);
      /* no break */
    case 1:
      CHECK_RNG(argv[0]);
      Data_Get_Struct(argv[0], gsl_rng, r);    
      return rb_float_new(gsl_ran_gaussian_ziggurat(r, sigma));
      break;
    default:
      rb_raise(rb_eArgError, "too many arguments (%d for 1 or 2)", argc);
      return Qnil;
      break;
    }
    break;
  default:
    Data_Get_Struct(obj, gsl_rng, r);    
    switch (argc) {
    case 2:
      n = NUM2INT(argv[1]);
      sigma = NUM2DBL(argv[0]);
      v = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) 
	gsl_vector_set(v, i, gsl_ran_gaussian_ziggurat(r, sigma));
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    case 1:
      sigma = NUM2DBL(argv[0]);
      /* no break */
    case 0:
      return rb_float_new(gsl_ran_gaussian_ziggurat(r, sigma));
      break;
    default:
      rb_raise(rb_eArgError, "too many arguments (%d for 0 or 1)", argc);
      return Qnil;
      break;
    }
    break;
  }
  return Qnil;
}

static VALUE rb_gsl_ran_gamma_mt(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_ran_eval2(argc, argv, obj, gsl_ran_gamma_mt);
}
#endif

void Init_gsl_ran(VALUE module)
{
  VALUE mgsl_ran;
  VALUE cgsl_ran_discrete;

  mgsl_ran = rb_define_module_under(module, "Ran");
  
  rb_define_module_function(mgsl_ran, "gaussian", rb_gsl_ran_gaussian, -1);
  rb_define_module_function(mgsl_ran, "ugaussian", rb_gsl_ran_gaussian, -1);
  rb_define_method(cgsl_rng, "gaussian", rb_gsl_ran_gaussian, -1);
  rb_define_alias(cgsl_rng, "ugaussian", "gaussian");

  rb_define_module_function(mgsl_ran, "gaussian_ratio_method", 
			    rb_gsl_ran_gaussian_ratio_method, -1);
  rb_define_module_function(mgsl_ran, "ugaussian_ratio_method",
			    rb_gsl_ran_gaussian_ratio_method, -1);
  rb_define_method(cgsl_rng, "gaussian_ratio_method", 
		   rb_gsl_ran_gaussian_ratio_method, -1);
  rb_define_alias(cgsl_rng, "ugaussian_ratio_method", "gaussian_ratio_method");
  rb_define_module_function(mgsl_ran,  "gaussian_pdf", rb_gsl_ran_gaussian_pdf, -1);
  rb_define_module_function(mgsl_ran,  "ugaussian_pdf", rb_gsl_ran_gaussian_pdf, -1);

  rb_define_module_function(mgsl_ran, "gaussian_tail", rb_gsl_ran_gaussian_tail, -1);
  rb_define_module_function(mgsl_ran, "ugaussian_tail", rb_gsl_ran_gaussian_tail, -1);
  rb_define_method(cgsl_rng, "gaussian_tail", rb_gsl_ran_gaussian_tail, -1);
  rb_define_alias(cgsl_rng, "ugaussian_tail", "gaussian_tail");
  rb_define_module_function(mgsl_ran,  "gaussian_tail_pdf", 
			    rb_gsl_ran_gaussian_tail_pdf, -1);
  rb_define_module_function(mgsl_ran,  "ugaussian_tail_pdf", 
			    rb_gsl_ran_gaussian_tail_pdf, -1);

  rb_define_module_function(mgsl_ran, "bivariate_gaussian", 
		   rb_gsl_ran_bivariate_gaussian, -1);
  rb_define_method(cgsl_rng, "bivariate_gaussian", 
		   rb_gsl_ran_bivariate_gaussian, -1);
  rb_define_module_function(mgsl_ran,  "bivariate_gaussian_pdf", 
			    rb_gsl_ran_bivariate_gaussian_pdf, 5);

  rb_define_module_function(mgsl_ran, "exponential", rb_gsl_ran_exponential, -1);
  rb_define_method(cgsl_rng, "exponential", rb_gsl_ran_exponential, -1);
  rb_define_module_function(mgsl_ran,  "exponential_pdf", 
			    rb_gsl_ran_exponential_pdf, 2);

  rb_define_module_function(mgsl_ran, "laplace", rb_gsl_ran_laplace, -1);
  rb_define_method(cgsl_rng, "laplace", rb_gsl_ran_laplace, -1);
  rb_define_module_function(mgsl_ran,  "laplace_pdf", rb_gsl_ran_laplace_pdf, 2);

  rb_define_module_function(mgsl_ran, "exppow", rb_gsl_ran_exppow, -1);
  rb_define_method(cgsl_rng, "exppow", rb_gsl_ran_exppow, -1);
  rb_define_module_function(mgsl_ran,  "exppow_pdf", rb_gsl_ran_exppow_pdf, 3);

  rb_define_module_function(mgsl_ran, "cauchy", rb_gsl_ran_cauchy, -1);
  rb_define_method(cgsl_rng, "cauchy", rb_gsl_ran_cauchy, -1);
  rb_define_module_function(mgsl_ran,  "cauchy_pdf", rb_gsl_ran_cauchy_pdf, 2);

  rb_define_module_function(mgsl_ran, "rayleigh", rb_gsl_ran_rayleigh, -1);
  rb_define_method(cgsl_rng, "rayleigh", rb_gsl_ran_rayleigh, -1);
  rb_define_module_function(mgsl_ran,  "rayleigh_pdf", rb_gsl_ran_rayleigh_pdf, 2);

  rb_define_module_function(mgsl_ran, "rayleigh_tail", rb_gsl_ran_rayleigh_tail, -1);
  rb_define_method(cgsl_rng, "rayleigh_tail", rb_gsl_ran_rayleigh_tail, -1);
  rb_define_module_function(mgsl_ran,  "rayleigh_tail_pdf", 
			    rb_gsl_ran_rayleigh_tail_pdf, 3);

  rb_define_module_function(mgsl_ran, "landau", rb_gsl_ran_landau, -1);
  rb_define_method(cgsl_rng, "landau", rb_gsl_ran_landau, -1);
  rb_define_module_function(mgsl_ran,  "landau_pdf", rb_gsl_ran_landau_pdf, 1);

  rb_define_method(cgsl_rng, "levy", rb_gsl_ran_levy, -1);
  rb_define_method(cgsl_rng, "levy_skew", rb_gsl_ran_levy_skew, -1);
  rb_define_module_function(mgsl_ran, "levy", rb_gsl_ran_levy, -1);
  rb_define_module_function(mgsl_ran, "levy_skew", rb_gsl_ran_levy_skew, -1);

  rb_define_module_function(mgsl_ran, "gamma", rb_gsl_ran_gamma, -1);
  rb_define_method(cgsl_rng, "gamma", rb_gsl_ran_gamma, -1);
  rb_define_module_function(mgsl_ran,  "gamma_pdf", rb_gsl_ran_gamma_pdf, 3);

  rb_define_module_function(mgsl_ran, "flat", rb_gsl_ran_flat, -1);
  rb_define_method(cgsl_rng, "flat", rb_gsl_ran_flat, -1);
  rb_define_module_function(mgsl_ran,  "flat_pdf", rb_gsl_ran_flat_pdf, 3);

  rb_define_module_function(mgsl_ran, "lognormal", rb_gsl_ran_lognormal, -1);
  rb_define_method(cgsl_rng, "lognormal", rb_gsl_ran_lognormal, -1);
  rb_define_module_function(mgsl_ran,  "lognormal_pdf", rb_gsl_ran_lognormal_pdf, 3);

  rb_define_module_function(mgsl_ran, "chisq", rb_gsl_ran_chisq, -1);
  rb_define_method(cgsl_rng, "chisq", rb_gsl_ran_chisq, -1);
  rb_define_module_function(mgsl_ran,  "chisq_pdf", rb_gsl_ran_chisq_pdf, 2);

  rb_define_module_function(mgsl_ran, "fdist", rb_gsl_ran_fdist, -1);
  rb_define_method(cgsl_rng, "fdist", rb_gsl_ran_fdist, -1);
  rb_define_module_function(mgsl_ran,  "fdist_pdf", rb_gsl_ran_fdist_pdf, 3);

  rb_define_module_function(mgsl_ran, "tdist", rb_gsl_ran_tdist, -1);
  rb_define_method(cgsl_rng, "tdist", rb_gsl_ran_tdist, -1);
  rb_define_module_function(mgsl_ran,  "tdist_pdf", rb_gsl_ran_tdist_pdf, 2);

  rb_define_module_function(mgsl_ran, "beta", rb_gsl_ran_beta, -1);
  rb_define_method(cgsl_rng, "beta", rb_gsl_ran_beta, -1);
  rb_define_module_function(mgsl_ran,  "beta_pdf", rb_gsl_ran_beta_pdf, 3);

  rb_define_module_function(mgsl_ran, "logistic", rb_gsl_ran_logistic, -1);
  rb_define_method(cgsl_rng, "logistic", rb_gsl_ran_logistic, -1);
  rb_define_module_function(mgsl_ran,  "logistic_pdf", rb_gsl_ran_logistic_pdf, 2);

  rb_define_module_function(mgsl_ran, "pareto", rb_gsl_ran_pareto, -1);
  rb_define_method(cgsl_rng, "pareto", rb_gsl_ran_pareto, -1);
  rb_define_module_function(mgsl_ran,  "pareto_pdf", rb_gsl_ran_pareto_pdf, 3);

  rb_define_module_function(mgsl_ran, "weibull", rb_gsl_ran_weibull, -1);
  rb_define_method(cgsl_rng, "weibull", rb_gsl_ran_weibull, -1);
  rb_define_module_function(mgsl_ran,  "weibull_pdf", rb_gsl_ran_weibull_pdf, 3);

  rb_define_module_function(mgsl_ran, "gumbel1", rb_gsl_ran_gumbel1, -1);
  rb_define_method(cgsl_rng, "gumbel1", rb_gsl_ran_gumbel1, -1);
  rb_define_module_function(mgsl_ran,  "gumbel1_pdf", rb_gsl_ran_gumbel1_pdf, 3);

  rb_define_module_function(mgsl_ran, "gumbel2", rb_gsl_ran_gumbel2, -1);
  rb_define_method(cgsl_rng, "gumbel2", rb_gsl_ran_gumbel2, -1);
  rb_define_module_function(mgsl_ran,  "gumbel2_pdf", rb_gsl_ran_gumbel2_pdf, 3);

  rb_define_module_function(mgsl_ran, "poisson", rb_gsl_ran_poisson, -1);
  rb_define_method(cgsl_rng, "poisson", rb_gsl_ran_poisson, -1);
  rb_define_module_function(mgsl_ran,  "poisson_pdf", rb_gsl_ran_poisson_pdf, 2);

  rb_define_module_function(mgsl_ran, "bernoulli", rb_gsl_ran_bernoulli, -1);
  rb_define_method(cgsl_rng, "bernoulli", rb_gsl_ran_bernoulli, -1);
  rb_define_module_function(mgsl_ran,  "bernoulli_pdf", rb_gsl_ran_bernoulli_pdf, 2);

  rb_define_module_function(mgsl_ran, "binomial", rb_gsl_ran_binomial, -1);
  rb_define_method(cgsl_rng, "binomial", rb_gsl_ran_binomial, -1);
#ifdef GSL_1_4_LATER
  rb_define_module_function(mgsl_ran, "binomial_tpe", rb_gsl_ran_binomial_tpe, -1);
  rb_define_method(cgsl_rng, "binomial_tpe", rb_gsl_ran_binomial_tpe, -1);
#endif
  rb_define_module_function(mgsl_ran,  "binomial_pdf", rb_gsl_ran_binomial_pdf, 3);

  rb_define_module_function(mgsl_ran, "negative_binomial", 
			    rb_gsl_ran_negative_binomial, -1);
  rb_define_method(cgsl_rng, "negative_binomial", rb_gsl_ran_negative_binomial, -1);
  rb_define_module_function(mgsl_ran,  "negative_binomial_pdf", rb_gsl_ran_negative_binomial_pdf, 3);

  rb_define_module_function(mgsl_ran, "pascal", rb_gsl_ran_pascal, -1);
  rb_define_method(cgsl_rng, "pascal", rb_gsl_ran_pascal, -1);
  rb_define_module_function(mgsl_ran,  "pascal_pdf", rb_gsl_ran_pascal_pdf, 3);

  rb_define_module_function(mgsl_ran, "geometric", rb_gsl_ran_geometric, -1);
  rb_define_method(cgsl_rng, "geometric", rb_gsl_ran_geometric, -1);
  rb_define_module_function(mgsl_ran,  "geometric_pdf", rb_gsl_ran_geometric_pdf, 2);

  rb_define_module_function(mgsl_ran, "hypergeometric", rb_gsl_ran_hypergeometric, -1);
  rb_define_method(cgsl_rng, "hypergeometric", rb_gsl_ran_hypergeometric, -1);
  rb_define_module_function(mgsl_ran,  "hypergeometric_pdf", 
			    rb_gsl_ran_hypergeometric_pdf, 4);

  rb_define_module_function(mgsl_ran, "logarithmic", rb_gsl_ran_logarithmic, -1);
  rb_define_method(cgsl_rng, "logarithmic", rb_gsl_ran_logarithmic, -1);
  rb_define_module_function(mgsl_ran,  "logarithmic_pdf",
			    rb_gsl_ran_logarithmic_pdf, 2);

  rb_define_module_function(mgsl_ran, "dir_2d", rb_gsl_ran_dir_2d, -1);
  rb_define_module_function(mgsl_ran, "dir_2d_trig_method", 
			    rb_gsl_ran_dir_2d_trig_method, -1);
  rb_define_module_function(mgsl_ran, "dir_3d", rb_gsl_ran_dir_3d, -1);
  rb_define_module_function(mgsl_ran, "dir_nd", rb_gsl_ran_dir_nd, -1);

  rb_define_method(cgsl_rng, "dir_2d", rb_gsl_ran_dir_2d, -1);
  rb_define_method(cgsl_rng, "dir_2d_trig_method", rb_gsl_ran_dir_2d_trig_method, -1);
  rb_define_method(cgsl_rng, "dir_3d", rb_gsl_ran_dir_3d, -1);
  rb_define_method(cgsl_rng, "dir_nd", rb_gsl_ran_dir_nd, -1);

  rb_define_method(cgsl_rng, "shuffle", rb_gsl_ran_shuffle, -1);
  rb_define_module_function(mgsl_ran, "shuffle", rb_gsl_ran_shuffle, -1);
  rb_define_method(cgsl_rng, "choose", rb_gsl_ran_choose, -1);
  rb_define_singleton_method(mgsl_ran, "choose", rb_gsl_ran_choose_singleton, -1);
  rb_define_method(cgsl_rng, "sample", rb_gsl_ran_sample, 2);

  /*****/

  cgsl_ran_discrete = rb_define_class_under(mgsl_ran, "Discrete", cGSL_Object);
  rb_define_singleton_method(cgsl_ran_discrete, "alloc", rb_gsl_ran_discrete_new, 1);
  rb_define_singleton_method(cgsl_ran_discrete, "preproc", rb_gsl_ran_discrete_new, 1);
  rb_define_method(cgsl_rng, "discrete", rb_gsl_ran_discrete, 1);
  rb_define_module_function(mgsl_ran,  "discrete_pdf", rb_gsl_ran_discrete_pdf, 2);

#ifdef GSL_1_3_LATER
  rb_define_method(cgsl_rng, "dirichlet", rb_gsl_ran_dirichlet, -1);
  rb_define_module_function(mgsl_ran,  "dirichlet_pdf", rb_gsl_ran_dirichlet_pdf, 2);
  rb_define_module_function(mgsl_ran,  "dirichlet_lnpdf", rb_gsl_ran_dirichlet_lnpdf, 2);
#endif

  /*  rb_define_method(cgsl_rng, "erlang", rb_gsl_ran_erlang, 2);
  rb_define_method(module, "ran_erlang_pdf", rb_gsl_ran_erlang_pdf, 3);
  rb_define_method(mgsl_ran, "erlang_pdf", rb_gsl_ran_erlang_pdf, 3);*/

#ifdef GSL_1_8_LATER
  rb_define_module_function(mgsl_ran, "gaussian_ziggurat", rb_gsl_ran_gaussian_ziggurat, -1);
  rb_define_method(cgsl_rng, "gaussian_ziggurat", rb_gsl_ran_gaussian_ziggurat, -1);
  rb_define_module_function(mgsl_ran, "gamma_mt", rb_gsl_ran_gamma_mt, -1);
  rb_define_method(cgsl_rng, "gamma_mt", rb_gsl_ran_gamma_mt, -1);
#endif

}
