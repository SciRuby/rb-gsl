/*
  rb_gsl_sf.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY
*/
#ifndef ___RB_GSL_SF_H___
#define ___RB_GSL_SF_H___

#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include "rb_gsl.h"

EXTERN VALUE cgsl_sf_result, cgsl_sf_result_e10;

VALUE rb_gsl_sf_result_new(VALUE klass);

VALUE rb_gsl_sf_eval1(double (*func)(double), VALUE argv);
VALUE rb_gsl_sf_eval_int_double(double (*func)(int, double), VALUE jj, VALUE argv);
VALUE rb_gsl_sf_eval_double_double(double (*func)(double, double), VALUE ff, VALUE argv);
VALUE rb_gsl_sf_eval1_uint(double (*func)(unsigned int), VALUE argv);
VALUE rb_gsl_sf_eval_double_m(double (*func)(double, gsl_mode_t), VALUE argv, VALUE m);;
VALUE rb_gsl_sf_eval_int_double_double(double (*func)(int, double, double), VALUE jj, 
				       VALUE ff, VALUE argv);
VALUE rb_gsl_sf_eval_int_int_double(double (*func)(int, int, double), VALUE jj,
				    VALUE jj2, VALUE argv);
VALUE rb_gsl_sf_eval_double3(double (*func)(double, double, double), 
					   VALUE ff, VALUE ff2, VALUE argv);
VALUE rb_gsl_sf_eval_double4(double (*func)(double, double, double, double), 
			     VALUE ff, VALUE ff2, VALUE ff3, VALUE argv);
VALUE rb_gsl_sf_eval_double_int(double (*func)(double, int), VALUE argv, VALUE jj);
VALUE rb_gsl_sf_eval1_int(double (*func)(int), VALUE argv);
VALUE rb_gsl_sf_eval_e(int (*func)(double, gsl_sf_result*), VALUE x);
VALUE rb_gsl_sf_eval_e_int_double(int (*func)(int, double, gsl_sf_result*), 
				  VALUE n, VALUE x);
VALUE rb_gsl_sf_eval_e_int_int_double(int (*func)(int, int, double, gsl_sf_result*), 
				      VALUE n1, VALUE n2, VALUE x);
VALUE rb_gsl_sf_eval_e_double2(int (*func)(double, double, gsl_sf_result*), 
			       VALUE x1, VALUE x2);
VALUE rb_gsl_sf_eval_e_uint(int (*func)(unsigned int, gsl_sf_result*), VALUE x);
VALUE rb_gsl_sf_eval_e_int_uint(int (*func)(int, unsigned int, gsl_sf_result*), 
				VALUE n, VALUE x);
VALUE rb_gsl_sf_eval_e_double_uint(int (*func)(double, unsigned int, gsl_sf_result*), 
				   VALUE y, VALUE x);
VALUE rb_gsl_sf_eval_e_m(int (*func)(double, gsl_mode_t, gsl_sf_result*), 
			 VALUE x, VALUE m);
VALUE rb_gsl_sf_eval_e_double2_m(int (*func)(double, double, gsl_mode_t, gsl_sf_result*), 
				 VALUE x1, VALUE x2, VALUE m);
VALUE rb_gsl_sf_eval_e_double3_m(int (*func)(double, double, double, gsl_mode_t, gsl_sf_result*), 
				 VALUE x1, VALUE x2, VALUE x3, VALUE m);
VALUE rb_gsl_sf_eval_double2_m(double (*func)(double, double, gsl_mode_t), 
			       VALUE argv, VALUE x2, VALUE m);
VALUE rb_gsl_sf_eval_double3_m(double (*func)(double, double, double, gsl_mode_t), 
			       VALUE argv, VALUE x2, VALUE x3, VALUE m);

VALUE rb_gsl_sf_eval_e_double4_m(int (*func)(double, double, double, double, gsl_mode_t, gsl_sf_result*), 
				 VALUE x1, VALUE x2, VALUE x3, VALUE x4, VALUE m);

VALUE rb_gsl_sf_eval_e_int(int (*func)(int, gsl_sf_result*), VALUE x);
VALUE rb_gsl_sf_eval_e_double3(int (*func)(double, double, double, gsl_sf_result*), 
			      VALUE x1, VALUE x2, VALUE x3);
VALUE rb_gsl_sf_eval_e_int_double2(int (*func)(int, double, double, gsl_sf_result*), 
				   VALUE n, VALUE x1, VALUE x2);
VALUE rb_gsl_sf_eval_e_double2(int (*func)(double, double, gsl_sf_result*), 
			       VALUE x1, VALUE x2);

VALUE eval_sf(double (*func)(double, gsl_mode_t), VALUE argv);

VALUE rb_gsl_sf_eval_double4_m(double (*func)(double, double, double, double,
					      gsl_mode_t), 
			       VALUE argv, VALUE x2, VALUE x3, VALUE x4, VALUE m);

VALUE rb_gsl_sf_eval_complex(double (*f)(double), VALUE obj);

void Init_gsl_sf_airy(VALUE module);
void Init_gsl_sf_bessel(VALUE module);
void Init_gsl_sf_clausen(VALUE module);
void Init_gsl_sf_coulomb(VALUE module);
void Init_gsl_sf_coupling(VALUE module);
void Init_gsl_sf_dawson(VALUE module);
void Init_gsl_sf_debye(VALUE module);
void Init_gsl_sf_dilog(VALUE module);
void Init_gsl_sf_elementary(VALUE module);
void Init_gsl_sf_ellint(VALUE module);
void Init_gsl_sf_elljac(VALUE module);
void Init_gsl_sf_erfc(VALUE module);
void Init_gsl_sf_exp(VALUE module);
void Init_gsl_sf_expint(VALUE module);
void Init_gsl_sf_fermi_dirac(VALUE module);
void Init_gsl_sf_gamma(VALUE module);
void Init_gsl_sf_gegenbauer(VALUE module);
void Init_gsl_sf_hyperg(VALUE module);
void Init_gsl_sf_laguerre(VALUE module);
void Init_gsl_sf_lambert(VALUE module);
void Init_gsl_sf_legendre(VALUE module);
void Init_gsl_sf_log(VALUE module);
void Init_gsl_sf_power(VALUE module);
void Init_gsl_sf_psi(VALUE module);
void Init_gsl_sf_synchrotron(VALUE module);
void Init_gsl_sf_transport(VALUE module);
void Init_gsl_sf_trigonometric(VALUE module);
void Init_gsl_sf_zeta(VALUE module);


#ifdef GSL_1_9_LATER
#include <gsl/gsl_sf_mathieu.h>
void Init_sf_mathieu(VALUE module);
#endif

#endif
