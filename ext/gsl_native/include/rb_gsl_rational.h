/*
  rb_gsl_rational.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or POLYNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef ___RB_GSL_RATIONAL_H___
#define ___RB_GSL_RATIONAL_H___

#include "rb_gsl_poly.h"

typedef struct ___gsl_rational
{
  VALUE num, den;
  gsl_poly *pnum;
  gsl_poly *pden;
} gsl_rational;

gsl_rational* gsl_rational_alloc();
gsl_rational* gsl_rational_new(const gsl_poly *num, const gsl_poly *den);
gsl_rational* gsl_rational_new2(const gsl_poly *num, const gsl_poly *den);
void gsl_rational_free(gsl_rational *r);

#endif
