/*
  rb_gsl_histogram.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef ___RB_GSL_HISTOGRAM_H___
#define ___RB_GSL_HISTOGRAM_H___

#include <gsl/gsl_vector.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include "rb_gsl.h"

EXTERN VALUE cgsl_histogram;
EXTERN VALUE cgsl_histogram_range;
EXTERN VALUE cgsl_histogram_bin;
EXTERN VALUE cgsl_histogram2d;
EXTERN VALUE cgsl_histogram2d_view;

typedef struct {
  gsl_histogram h;
} mygsl_histogram2d_view;

#ifndef HISTOGRAM2D_P
#define HISTOGRAM2D_P(x) (rb_obj_is_kind_of(x,cgsl_histogram2d)?1:0)
#endif

#ifndef CHECK_HISTOGRAM2D
#define CHECK_HISTOGRAM2D(x) if(!rb_obj_is_kind_of(x,cgsl_histogram2d))\
    rb_raise(rb_eTypeError, "wrong type (Histogram2d expected)");
#endif


#ifndef HISTOGRAM3D_P
#define HISTOGRAM3D_P(x) (rb_obj_is_kind_of(x,cgsl_histogram3d)?1:0)
#endif

#ifndef CHECK_HISTOGRAM3D
#define CHECK_HISTOGRAM3D(x) if(!rb_obj_is_kind_of(x,cgsl_histogram3d))\
    rb_raise(rb_eTypeError, "wrong type (Histogram3d expected)");
#endif

#include "rb_gsl_histogram3d.h"

int
mygsl_histogram_equal_bins_p (const gsl_histogram * h1, const gsl_histogram * h2);
int 
mygsl_histogram_add (gsl_histogram * h1, const gsl_histogram * h2);
int 
mygsl_histogram_sub (gsl_histogram * h1, const gsl_histogram * h2);
int 
mygsl_histogram_mul (gsl_histogram * h1, const gsl_histogram * h2);
int 
mygsl_histogram_div (gsl_histogram * h1, const gsl_histogram * h2);

#endif
