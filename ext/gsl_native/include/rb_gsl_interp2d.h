/*
  rb_gsl_interp2d.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada
    (C) Copyright 2015 - by Ruby Science Foundation

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or INTERPOLATIONNESS FOR A PARTICULAR PURPOSE.
*/
#ifdef GSL_2_0_LATER
#ifndef ___RB_GSL_INTERP2D_H___
#define ___RB_GSL_INTERP2D_H___

#include "rb_gsl.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

typedef struct {
  gsl_interp2d *p;
  gsl_interp_accel *xacc, *yacc;
} rb_gsl_interp2d;

typedef struct {
  gsl_spline2d *s;
  gsl_interp_accel *xacc, *yacc;
} rb_gsl_spline2d;

enum {
  GSL_INTERP2D_BICUBIC,
  GSL_INTERP2D_BILINEAR
};

const gsl_interp2d_type* get_interp2d_type(VALUE);
static void rb_gsl_interp2d_free(rb_gsl_interp2d*);

#endif
#endif