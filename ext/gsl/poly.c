/*
  poly.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_poly.h"
#include "include/rb_gsl_array.h"
#include "include/rb_gsl_common.h"
#ifdef HAVE_NARARY_H
#include "narray.h"
#endif

void Init_gsl_poly_init(VALUE module);
void Init_gsl_poly_int_init(VALUE module);
void Init_gsl_poly2(VALUE module);

#define BASE_DOUBLE
#include "include/templates_on.h"
#include "poly_source.c"
#include "include/templates_off.h"
void Init_gsl_poly(VALUE module)
{
  Init_gsl_poly_init(module);
}

#undef  BASE_DOUBLE

#define BASE_INT
#include "include/templates_on.h"
#include "poly_source.c"
#include "include/templates_off.h"
void Init_gsl_poly_int(VALUE module)
{
  Init_gsl_poly_int_init(module);
}
#undef  BASE_INT
