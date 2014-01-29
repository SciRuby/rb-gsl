/*
  rb_gsl_poly.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or POLYNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef ___RB_GSL_POLY_H___
#define ___RB_GSL_POLY_H___

#include <gsl/gsl_poly.h>
#include "rb_gsl_complex.h"
#include "rb_gsl_array.h"

EXTERN VALUE cgsl_poly;
EXTERN VALUE cgsl_poly_int;
EXTERN VALUE cgsl_poly_dd;
EXTERN VALUE cgsl_poly_taylor;
EXTERN VALUE cgsl_poly_workspace;
EXTERN VALUE cgsl_rational;

typedef gsl_vector gsl_poly;
typedef gsl_vector_int gsl_poly_int;
/*
typedef struct ___gsl_rational 
{
  VALUE num, den;
  gsl_poly *pnum;
  gsl_poly *pden;
} gsl_rational;
*/

int gsl_poly_conv(const double *a, size_t na, const double *b, size_t nb,
		  double *c, size_t *nc);

gsl_vector* gsl_poly_deconv_vector(const gsl_vector *c, const gsl_vector *a, gsl_vector **r);
gsl_vector* gsl_poly_deriv(const gsl_vector *v);
gsl_vector* gsl_poly_integ(const gsl_vector *v);
gsl_vector* gsl_poly_reduce(const gsl_vector *v);
VALUE rb_gsl_poly_complex_solve(int argc, VALUE *argv, VALUE obj);
gsl_vector* gsl_poly_conv_vector(const gsl_vector *v1, const gsl_vector *v2);
gsl_poly* gsl_poly_add(const gsl_poly *a, const gsl_poly *b);
VALUE rb_gsl_poly_deconv(VALUE obj, VALUE bb);

VALUE rb_gsl_poly_complex_solve2(int argc, VALUE *argv, VALUE obj);


int gsl_poly_int_conv(const int *a, size_t na, const int *b, size_t nb,
		  int *c, size_t *nc);

gsl_vector_int* gsl_poly_int_deconv_vector(const gsl_vector_int *c, const gsl_vector_int *a, gsl_vector_int **r);
gsl_vector_int* gsl_poly_int_deriv(const gsl_vector_int *v);
gsl_vector_int* gsl_poly_int_integ(const gsl_vector_int *v);
gsl_vector_int* gsl_poly_int_reduce(const gsl_vector_int *v);
VALUE rb_gsl_poly_int_complex_solve(int argc, VALUE *argv, VALUE obj);
gsl_vector_int* gsl_poly_int_conv_vector(const gsl_vector_int *v1, const gsl_vector_int *v2);
gsl_poly_int* gsl_poly_int_add(const gsl_poly_int *a, const gsl_poly_int *b);
VALUE rb_gsl_poly_int_deconv(VALUE obj, VALUE bb);

gsl_poly* get_poly_get(VALUE obj, int *flag);
gsl_poly_int* get_poly_int_get(VALUE obj, int *flag);
#endif
