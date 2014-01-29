/*
  rb_gsl_linalg.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY
*/

#ifndef ___RB_GSL_LINALG_H___
#define ___RB_GSL_LINALG_H___

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include "rb_gsl_with_narray.h"

VALUE rb_gsl_linalg_complex_LU_decomp(int argc, VALUE *argv, VALUE obj);
VALUE rb_gsl_linalg_complex_LU_decomp2(int argc, VALUE *argv, VALUE obj);

#endif
