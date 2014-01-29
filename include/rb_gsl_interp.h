/*
  rb_gsl_interp.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or INTERPOLATIONNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef ___RB_GSL_INTERP_H___
#define ___RB_GSL_INTERP_H___

#include "rb_gsl.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif


typedef struct {
  gsl_interp *p;
  gsl_interp_accel *a;
} rb_gsl_interp;

typedef struct {
  gsl_spline *s;
  gsl_interp_accel *a;
} rb_gsl_spline;

enum {
  GSL_INTERP_LINEAR,
  GSL_INTERP_POLYNOMIAL,
  GSL_INTERP_CSPLINE,
  GSL_INTERP_CSPLINE_PERIODIC,
  GSL_INTERP_AKIMA,
  GSL_INTERP_AKIMA_PERIODIC,
};

const gsl_interp_type* get_interp_type(VALUE t);

#endif
