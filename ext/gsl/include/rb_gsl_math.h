/*
  rb_gsl_math.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef ___RB_GSL_MATH_H___
#define ___RB_GSL_MATH_H___

#include "rb_gsl.h"

VALUE rb_gsl_math_complex_eval(gsl_complex (*func)(gsl_complex), VALUE obj);

#endif
