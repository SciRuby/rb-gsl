/*
  sf_elljac.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "include/rb_gsl_sf.h"

static VALUE rb_gsl_sf_elljac_e(VALUE obj, VALUE n, VALUE m)
{
  double sn, cn, dn;
  // local variable "status" declared and set, but never used
  //int status;
  Need_Float(n); Need_Float(m);
  /*status =*/ gsl_sf_elljac_e(NUM2DBL(n), NUM2DBL(m), &sn, &cn, &dn);
  return rb_ary_new3(3, rb_float_new(sn),
         rb_float_new(cn), rb_float_new(dn));
}

void Init_gsl_sf_elljac(VALUE module)
{
  rb_define_module_function(module, "elljac_e",  rb_gsl_sf_elljac_e, 2);
  rb_define_module_function(module, "elljac",  rb_gsl_sf_elljac_e, 2);
}
