/*
  blas.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include <gsl/gsl_blas.h>
#include "include/rb_gsl_common.h"
#include "include/rb_gsl_array.h"

void Init_gsl_blas1(VALUE module);
void Init_gsl_blas2(VALUE module);
void Init_gsl_blas3(VALUE module);

void Init_gsl_blas(VALUE module)
{
  VALUE mgsl_blas;
  mgsl_blas = rb_define_module_under(module, "Blas");

  Init_gsl_blas1(mgsl_blas);
  Init_gsl_blas2(mgsl_blas);
  Init_gsl_blas3(mgsl_blas);
}
