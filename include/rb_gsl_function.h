/*
  rb_gsl_function.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef ___RB_GSL_FUNCTION_H___
#define ___RB_GSL_FUNCTION_H___

#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include "rb_gsl.h"

EXTERN VALUE cgsl_function;
EXTERN VALUE cgsl_function_fdf;
extern ID RBGSL_ID_call, RBGSL_ID_arity;
void gsl_function_mark(gsl_function *f);
void gsl_function_free(gsl_function *f);
#endif
