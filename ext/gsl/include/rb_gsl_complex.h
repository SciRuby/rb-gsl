/*
  rb_gsl_complex.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef ___RB_GSL_COMPLEX_H___
#define ___RB_GSL_COMPLEX_H___

#include <math.h>
#include <ruby.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

EXTERN VALUE cgsl_complex;
VALUE rb_gsl_complex_pow(int argc, VALUE *argv, VALUE obj);
VALUE rb_gsl_complex_pow_real(int argc, VALUE *argv, VALUE obj);

#endif
