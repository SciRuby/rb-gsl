/*
  poly2.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_poly.h"
#include "rb_gsl_math.h"
#include "rb_gsl_array.h"
#include <gsl/gsl_sf_gamma.h>

static gsl_poly_int* mygsl_poly_hermite(int n1)
{
  size_t n;
  gsl_vector_int *p1, *p2, *p0;
  int coef1[2] = {0, 2};
  int coef2[3] = {-2, 0, 4};
  if (n1 < 0) rb_raise(rb_eArgError, "order must be >= 0");
  p0 = gsl_vector_int_calloc(n1 + 1);
  switch (n1) {
  case 0:
    gsl_vector_int_set(p0, 0, 1);
    break;
  case 1:
    memcpy(p0->data, coef1, 2*sizeof(int));
    break;
  case 2:
    memcpy(p0->data, coef2, 3*sizeof(int));
    break;
  default:
    p1 = gsl_vector_int_calloc(n1 + 1);
    p2 = gsl_vector_int_calloc(n1 + 1);
    memcpy(p1->data, coef2, 3*sizeof(int));
    memcpy(p2->data, coef1, 2*sizeof(int));
    for (n = 2; (int) n < n1; n++) {
      gsl_vector_int_memcpy(p0, p1);
      mygsl_vector_int_shift_scale2(p0, n);
      gsl_vector_int_scale(p2, 2*n);

      gsl_vector_int_sub(p0, p2);
      /* save for the next iteration */
      gsl_vector_int_memcpy(p2, p1);
      gsl_vector_int_memcpy(p1, p0);
    }
    gsl_vector_int_free(p2);
    gsl_vector_int_free(p1);    
    break;
  }
  return p0;
}

static gsl_poly_int* mygsl_poly_cheb(int n1)
{
  size_t n;
  gsl_vector_int *p1, *p2, *p0;
  int coef1[2] = {0, 1};
  int coef2[3] = {-1, 0, 2};
  if (n1 < 0) rb_raise(rb_eArgError, "order must be >= 0");
  p0 = gsl_vector_int_calloc(n1 + 1);
  switch (n1) {
  case 0:
    gsl_vector_int_set(p0, 0, 1);
    break;
  case 1:
    memcpy(p0->data, coef1, 2*sizeof(int));
    break;
  case 2:
    memcpy(p0->data, coef2, 3*sizeof(int));
    break;
  default:
    p1 = gsl_vector_int_calloc(n1 + 1);
    p2 = gsl_vector_int_calloc(n1 + 1);
    memcpy(p1->data, coef2, 3*sizeof(int));
    memcpy(p2->data, coef1, 2*sizeof(int));
    for (n = 2; (int) n < n1; n++) {
      gsl_vector_int_memcpy(p0, p1);
      mygsl_vector_int_shift_scale2(p0, n);

      gsl_vector_int_sub(p0, p2);
      /* save for the next iteration */
      gsl_vector_int_memcpy(p2, p1);
      gsl_vector_int_memcpy(p1, p0);
    }
    gsl_vector_int_free(p2);
    gsl_vector_int_free(p1);    
    break;
  }
  return p0;
}

static gsl_poly_int* mygsl_poly_chebII(int n1)
{
  size_t n;
  gsl_vector_int *p1, *p2, *p0;
  int coef1[2] = {0, 2};
  int coef2[3] = {-1, 0, 4};
  if (n1 < 0) rb_raise(rb_eArgError, "order must be >= 0");
  p0 = gsl_vector_int_calloc(n1 + 1);
  switch (n1) {
  case 0:
    gsl_vector_int_set(p0, 0, 1);
    break;
  case 1:
    memcpy(p0->data, coef1, 2*sizeof(int));
    break;
  case 2:
    memcpy(p0->data, coef2, 3*sizeof(int));
    break;
  default:
    p1 = gsl_vector_int_calloc(n1 + 1);
    p2 = gsl_vector_int_calloc(n1 + 1);
    memcpy(p1->data, coef2, 3*sizeof(int));
    memcpy(p2->data, coef1, 2*sizeof(int));
    for (n = 2; (int) n < n1; n++) {
      gsl_vector_int_memcpy(p0, p1);
      mygsl_vector_int_shift_scale2(p0, n);
      gsl_vector_int_sub(p0, p2);
      /* save for the next iteration */
      gsl_vector_int_memcpy(p2, p1);
      gsl_vector_int_memcpy(p1, p0);
    }
    gsl_vector_int_free(p2);
    gsl_vector_int_free(p1);    
    break;
  }
  return p0;
}

static gsl_poly_int* mygsl_poly_laguerre(int n)
{
  size_t m, k;
  int val;
  gsl_vector_int *p0;
  if (n < 0) rb_raise(rb_eArgError, "order must be >= 0");
  p0 = gsl_vector_int_calloc(n + 1);
  switch (n) {
  case 0:
    gsl_vector_int_set(p0, 0, 1);
    break;
  case 1:
    gsl_vector_int_set(p0, 0, 1);
    gsl_vector_int_set(p0, 1, -1);
    break;
  default:
    k = gsl_sf_fact(n);
    for (m = 0; (int) m <= n; m++) {
      val = k*k/gsl_sf_fact(n-m)/gsl_pow_2(gsl_sf_fact(m));
      if (m%2 == 1) val *= -1;
      gsl_vector_int_set(p0, m, val);
    }
    break;
  }
  return p0;
}

static gsl_poly_int* mygsl_poly_bessel(int n)
{
  size_t k;
  gsl_vector_int *p0;
  if (n < 0) rb_raise(rb_eArgError, "order must be >= 0");
  p0 = gsl_vector_int_calloc(n + 1);
  for (k = 0; (int) k <= n; k++) {
    gsl_vector_int_set(p0, k, gsl_sf_fact(n+k)/gsl_sf_fact(n-k)/gsl_sf_fact(k)/((int) pow(2, k)));
  }
  return p0;
}

static gsl_poly_int* mygsl_poly_bell(int n1)
{
  size_t n, j;
  gsl_vector_int *p1, *p0;
  int coef1[2] = {0, 1};
  int coef2[3] = {0, 1, 1};
  if (n1 < 0) rb_raise(rb_eArgError, "order must be >= 0");
  p0 = gsl_vector_int_calloc(n1 + 1);
  switch (n1) {
  case 0:
    gsl_vector_int_set(p0, 0, 1);
    break;
  case 1:
    memcpy(p0->data, coef1, 2*sizeof(int));
    break;
  case 2:
    memcpy(p0->data, coef2, 3*sizeof(int));
    break;
  default:
    p1 = gsl_vector_int_calloc(n1 + 1);
    memcpy(p1->data, coef2, 3*sizeof(int));
    for (n = 2; (int) n < n1; n++) {
      gsl_vector_int_memcpy(p0, p1);
      mygsl_vector_int_shift(p0, n);
      for (j = 0; j < n; j++) {
	gsl_vector_int_set(p1, j, gsl_vector_int_get(p1, j+1)*(j+1));
      }
      gsl_vector_int_set(p1, n, 0);
      mygsl_vector_int_shift(p1, n);
      gsl_vector_int_add(p0, p1);
      /* save for the next iteration */
      gsl_vector_int_memcpy(p1, p0);
    }
    gsl_vector_int_free(p1);    
    break;
  }
  return p0;
}

static VALUE rb_gsl_poly_define_poly(VALUE klass, VALUE order,
				     gsl_poly_int* (*f)(int n1)) {
  int n1;
  gsl_poly_int *pnew = NULL;
  CHECK_FIXNUM(order);
  n1 = FIX2INT(order);
  if (n1 < 0) rb_raise(rb_eArgError, "order must be >= 0");
  pnew = (*f)(n1);
  return Data_Wrap_Struct(cgsl_poly_int, 0, gsl_vector_int_free, pnew);
}

static VALUE rb_gsl_poly_hermite(VALUE klass, VALUE order)
{
  return rb_gsl_poly_define_poly(klass, order, mygsl_poly_hermite);
}

static VALUE rb_gsl_poly_cheb(VALUE klass, VALUE order)
{
  return rb_gsl_poly_define_poly(klass, order, mygsl_poly_cheb);
}

static VALUE rb_gsl_poly_chebII(VALUE klass, VALUE order)
{
  return rb_gsl_poly_define_poly(klass, order, mygsl_poly_chebII);
}

static VALUE rb_gsl_poly_laguerre(VALUE klass, VALUE order)
{
  return rb_gsl_poly_define_poly(klass, order, mygsl_poly_laguerre);
}

static VALUE rb_gsl_poly_bessel(VALUE klass, VALUE order)
{
  return rb_gsl_poly_define_poly(klass, order, mygsl_poly_bessel);
}

static VALUE rb_gsl_poly_bell(VALUE klass, VALUE order)
{
  return rb_gsl_poly_define_poly(klass, order, mygsl_poly_bell);
}

void Init_gsl_poly2(VALUE module)
{
  rb_define_singleton_method(cgsl_poly, "hermite", rb_gsl_poly_hermite, 1);
  rb_define_singleton_method(cgsl_poly, "cheb", rb_gsl_poly_cheb, 1);
  rb_define_singleton_method(cgsl_poly, "chebyshev", rb_gsl_poly_cheb, 1);
  rb_define_singleton_method(cgsl_poly, "cheb_I", rb_gsl_poly_cheb, 1);
  rb_define_singleton_method(cgsl_poly, "chebyshev_I", rb_gsl_poly_cheb, 1);
  rb_define_singleton_method(cgsl_poly, "cheb_II", rb_gsl_poly_chebII, 1);
  rb_define_singleton_method(cgsl_poly, "chebyshev_II", rb_gsl_poly_chebII, 1);
  rb_define_singleton_method(cgsl_poly, "bessel", rb_gsl_poly_bessel, 1);
  rb_define_singleton_method(cgsl_poly, "bell", rb_gsl_poly_bell, 1);
  rb_define_singleton_method(cgsl_poly, "laguerre", rb_gsl_poly_laguerre, 1);
}
