/*
  rb_gsl.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY
*/

#ifndef ___RB_GSL_H___
#define ___RB_GSL_H___

#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_mode.h>
#include "rb_gsl_poly.h"
#include "rb_gsl_math.h"
#include "rb_gsl_sf.h"
#include "rb_gsl_linalg.h"
#include "rb_gsl_eigen.h"
#include "rb_gsl_fft.h"
#include "rb_gsl_integration.h"
#include "rb_gsl_rng.h"
#include "rb_gsl_statistics.h"
#include "rb_gsl_histogram.h"
#include "rb_gsl_odeiv.h"
#include "rb_gsl_interp.h"
#include "rb_gsl_cheb.h"
#include "rb_gsl_root.h"
#include "rb_gsl_fit.h"
#include "rb_gsl_const.h"

void Init_gsl_error(VALUE module);
void Init_gsl_math(VALUE module);
void Init_gsl_complex(VALUE module);
void Init_gsl_array(VALUE module);
void Init_gsl_blas(VALUE module);
void Init_gsl_sort(VALUE module);
void Init_gsl_poly(VALUE module);
void Init_gsl_poly_int(VALUE module);
void Init_gsl_poly2(VALUE module);
void Init_gsl_rational(VALUE module);
void Init_gsl_sf(VALUE module);
void Init_gsl_linalg(VALUE module);
void Init_gsl_eigen(VALUE module);
void Init_gsl_fft(VALUE module);
void Init_gsl_signal(VALUE module);
void Init_gsl_function(VALUE module);
void Init_gsl_integration(VALUE module);

void Init_gsl_rng(VALUE module);
void Init_gsl_qrng(VALUE module);
void Init_gsl_ran(VALUE module);
void Init_gsl_cdf(VALUE module);
void Init_gsl_stats(VALUE module);

void Init_gsl_histogram(VALUE module);
void Init_gsl_histogram2d(VALUE module);
void Init_gsl_histogram3d(VALUE module);
void Init_gsl_ntuple(VALUE module);
void Init_gsl_monte(VALUE module);
void Init_gsl_siman(VALUE module);

void Init_gsl_odeiv(VALUE module);
void Init_gsl_interp(VALUE module);
void Init_gsl_spline(VALUE module);
void Init_gsl_diff(VALUE module);
#ifdef GSL_1_4_9_LATER
void Init_gsl_deriv(VALUE module);
#endif

void Init_gsl_cheb(VALUE module);
void Init_gsl_sum(VALUE module);
void Init_gsl_dht(VALUE module);

void Init_gsl_root(VALUE module);
void Init_gsl_multiroot(VALUE module);
void Init_gsl_min(VALUE module);
void Init_gsl_multimin(VALUE module);
void Init_gsl_fit(VALUE module);
void Init_gsl_multifit(VALUE module);

void Init_gsl_const(VALUE module);

void Init_gsl_ieee(VALUE module);

#ifdef HAVE_NARRAY_H
void Init_gsl_narray(VALUE module);
#endif

void Init_wavelet(VALUE module);

void Init_gsl_graph(VALUE module);

#ifdef HAVE_TENSOR_TENSOR_H
void Init_tensor_init(VALUE module);
void Init_tensor_int_init(VALUE module);
#endif

void Init_gsl_dirac(VALUE module);

EXTERN VALUE cGSL_Object;

void Init_tamu_anova(VALUE module);

#ifdef HAVE_TAMU_ANOVA_TAMU_ANOVA_H
#include "tamu_anova/tamu_anova.h"
#endif

#ifdef HAVE_OOL_OOL_VERSION_H
void Init_ool(VALUE module);
#endif

#ifdef HAVE_JACOBI_H
void Init_jacobi(VALUE module);
#endif

#ifdef HAVE_GSL_GSL_CQP_H
void Init_cqp(VALUE module);
#endif

void Init_fresnel(VALUE module);

#ifdef GSL_1_9_LATER
void Init_bspline(VALUE module);
#endif

#ifdef HAVE_ALF_ALF_H
#include "alf/alf.h"
#endif
void Init_alf(VALUE module);
void Init_geometry(VALUE module);

#ifdef GSL_1_14_LATER
#include <gsl/gsl_multiset.h>
extern VALUE cMultiset;
void Init_multiset(VALUE module);

#endif

#endif
