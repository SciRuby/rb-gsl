/*
  rational.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_rational.h"
#include "rb_gsl_array.h"

VALUE cgsl_rational;
static gsl_rational* gsl_rational_div_poly(const gsl_rational *r1, const gsl_poly *p);

static void gsl_rational_mark(gsl_rational *r);
gsl_rational* gsl_rational_alloc()
{
  gsl_rational *r = NULL;
  r = ALLOC(gsl_rational);
  r->num = (VALUE) NULL;
  r->den = (VALUE) NULL;
  return r;
}

gsl_rational* gsl_rational_new(const gsl_poly *num, const gsl_poly *den)
{
  gsl_rational *r = NULL;
  r = gsl_rational_alloc();
  r->pnum = make_vector_clone(num);
  r->pden = make_vector_clone(den);
  r->num = Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, r->pnum);
  r->den = Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, r->pden);
  return r;
}

gsl_rational* gsl_rational_new2(const gsl_poly *num, const gsl_poly *den)
{
  gsl_rational *r = NULL;
  r = gsl_rational_alloc();
  r->pnum = (gsl_poly *) num;
  r->pden = (gsl_poly *) den;
  r->num = Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, r->pnum);
  r->den = Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, r->pden);
  return r;
}

void gsl_rational_free(gsl_rational *r)
{
  free((gsl_rational *) r);
}


static gsl_rational* gsl_rational_add(const gsl_rational *r, const gsl_rational *r2)
{
  gsl_rational *rnew = NULL;
  gsl_poly *dennew = NULL, *numnew = NULL;
  gsl_poly *a = NULL, *b = NULL;
  int flag = 0;
  if (rbgsl_vector_equal(r->pden, r2->pden, 1e-10)) {
    dennew = r->pden;
    numnew = gsl_poly_add(r->pnum, r2->pnum);
    flag = 0;
  } else {
    dennew = gsl_poly_conv_vector(r->pden, r2->pden);
    a = gsl_poly_conv_vector(r->pden, r2->pnum);
    b = gsl_poly_conv_vector(r2->pden, r->pnum);
    numnew = gsl_poly_add(a, b);
    gsl_vector_free(a);
    gsl_vector_free(b);
    flag = 1;
  }
  rnew = gsl_rational_new(numnew, dennew);
  gsl_vector_free(numnew);
  if (flag == 1) gsl_vector_free(dennew);
  return rnew;
}

static gsl_rational* gsl_rational_add_poly(const gsl_rational *r, const gsl_poly *p)
{
  gsl_rational *rnew = NULL;
  gsl_poly *numnew = NULL;
  gsl_poly *a = NULL;
  a = gsl_poly_conv_vector(r->pden, p);
  numnew = gsl_poly_add(a, r->pnum);
  rnew = gsl_rational_new(numnew, r->pden);
  gsl_vector_free(a);
  gsl_vector_free(numnew);
  return rnew;
}

static gsl_rational* gsl_rational_mul(const gsl_rational *r1, const gsl_rational *r2)
{
  gsl_rational *rnew = NULL;
  gsl_poly *num = NULL, *den = NULL;
  num = gsl_poly_conv_vector(r1->pnum, r2->pnum);
  den = gsl_poly_conv_vector(r1->pden, r2->pden);
  rnew = gsl_rational_new2(num, den);
  return rnew;
}

static gsl_rational* gsl_rational_mul_poly(const gsl_rational *r1, const gsl_poly *p)
{
  gsl_rational *rnew = NULL;
  gsl_poly *num = NULL;
  num = gsl_poly_conv_vector(r1->pnum, p);
  rnew = gsl_rational_new(num, r1->pden);
  gsl_vector_free(num);
  return rnew;
}

static gsl_rational* gsl_rational_div(const gsl_rational *r1, const gsl_rational *r2)
{
  gsl_rational *rnew = NULL;
  gsl_poly *num = NULL, *den = NULL;
  num = gsl_poly_conv_vector(r1->pnum, r2->pden);
  den = gsl_poly_conv_vector(r1->pden, r2->pnum);
  rnew = gsl_rational_new2(num, den);
  return rnew;
}

static gsl_rational* gsl_rational_div_poly(const gsl_rational *r1, const gsl_poly *p)
{
  gsl_rational *rnew = NULL;
  gsl_poly *den = NULL;
  den = gsl_poly_conv_vector(r1->pden, p);
  rnew = gsl_rational_new(r1->pnum, den);
  gsl_vector_free(den);
  return rnew;
}

static gsl_rational* gsl_rational_inverse(const gsl_rational *r)
{
  return gsl_rational_new(r->pden, r->pnum);
}

/*
static gsl_poly* gsl_rational_partial_fraction(const gsl_rational *r)
{
  gsl_poly *num = NULL, *r = NULL;
  num = gsl_poly_deconv_vector(r->pnum, r->pden, &r);
  if (gsl_vector_isnull(r)) {
    gsl_vector_free(num);
    gsl_vector_free(r);
    return NULL;
  }

  gsl_vector_free(num);
  return r;
}
*/

/*****/
gsl_poly* get_poly_get(VALUE obj, int *flag);


static void gsl_rational_mark(gsl_rational *r)
{
  rb_gc_mark(r->num);
  rb_gc_mark(r->den);
}

static VALUE rb_gsl_rational_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_poly *den = NULL, *num = NULL;
  gsl_rational *r = NULL;
  int flag1 = 0, flag2 = 0;
  switch (argc) {
  case 0:
    r = gsl_rational_alloc();
    break;
  case 2:
    den = get_poly_get(argv[0], &flag1);
    num = get_poly_get(argv[1], &flag2);
    r = gsl_rational_new(den, num);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 2)", argc);
    break;
  }
  if (flag1 == 1) gsl_vector_free(den);
  if (flag2 == 1) gsl_vector_free(num);
  return Data_Wrap_Struct(klass, gsl_rational_mark, gsl_rational_free, r);
}
/* singleton */
static VALUE rb_gsl_poly_make_rational(VALUE obj, VALUE other)
{
  gsl_rational *rnew = NULL;
  gsl_poly *p, *p2;
  size_t i;
  Data_Get_Struct(obj, gsl_poly, p);
  if (VECTOR_P(other)) {
    Data_Get_Struct(other, gsl_vector, p2);
    rnew = gsl_rational_new(p, p2);
  } else {
    switch (TYPE(other)) {
    case T_ARRAY:
      p2 = gsl_vector_alloc(RARRAY_LEN(other));
      for (i = 0; i < p2->size; i++)
	gsl_vector_set(p2, i, NUM2DBL(rb_ary_entry(other, i)));
      rnew = gsl_rational_new(p, p2);
      gsl_vector_free(p2);
      break;
    case T_FLOAT:
    case T_FIXNUM:
      p2 = make_vector_clone(p);
      gsl_vector_scale(p2, 1.0/NUM2DBL(other));
      return Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, p2);
      break;
    default:
      rb_raise(rb_eTypeError, "wrong argument type %s", 
	       rb_class2name(CLASS_OF(other)));
      break;
    }
  }
  return Data_Wrap_Struct(cgsl_rational, gsl_rational_mark, gsl_rational_free, rnew);
}

static VALUE rb_gsl_rational_print(VALUE obj)
{
  gsl_rational *r = NULL;
  Data_Get_Struct(obj, gsl_rational, r);
  gsl_vector_print(r->pnum, cgsl_vector);
  gsl_vector_print(r->pden, cgsl_vector);
  return Qnil;
}

static VALUE rb_gsl_rational_to_s(VALUE obj)
{
  gsl_rational *r = NULL;
  VALUE str;
  Data_Get_Struct(obj, gsl_rational, r);
  str = rb_gsl_vector_to_s(r->num);
  rb_str_concat(str, rb_str_new2("\n"));
  rb_str_concat(str, rb_gsl_vector_to_s(r->den));
  return str;
}

static VALUE rb_gsl_rational_inspect(VALUE obj)
{
  VALUE str;
  str = rb_str_new2(rb_class2name(CLASS_OF(obj)));
  rb_str_concat(str, rb_str_new2("\n"));
  rb_str_concat(str, rb_gsl_rational_to_s(obj));
  return str;
}

static VALUE rb_gsl_rational_add(VALUE obj, VALUE other)
{
  gsl_rational *r = NULL, *r2 = NULL, *rnew = NULL;
  gsl_poly *p = NULL;
  int flag = 0;
  Data_Get_Struct(obj, gsl_rational, r);
  if (RATIONAL_P(other)) {
    Data_Get_Struct(other, gsl_rational, r2);
    rnew = gsl_rational_add(r, r2);
  } else {
    p = get_poly_get(other, &flag);
    rnew = gsl_rational_add_poly(r, p);
    if (flag == 1) gsl_vector_free(p);
  }
  return Data_Wrap_Struct(cgsl_rational, gsl_rational_mark, gsl_rational_free, rnew);
}

static VALUE rb_gsl_rational_deconv(VALUE obj)
{
  gsl_rational *r = NULL;
  Data_Get_Struct(obj, gsl_rational, r);
  return rb_gsl_poly_deconv(r->num, r->den);
}

static VALUE rb_gsl_rational_uminus(VALUE obj)
{
  gsl_rational *r = NULL, *rnew;
  gsl_poly *p = NULL, *ptmp;
  int flag = 0;
  size_t i;
  if (RATIONAL_P(obj)) {
    Data_Get_Struct(obj, gsl_rational, r);
    rnew = gsl_rational_new(r->pnum, r->pden);
    for (i = 0; i < rnew->pnum->size; i++) 
      gsl_vector_set(rnew->pnum, i, -gsl_vector_get(r->pnum, i));
    return Data_Wrap_Struct(cgsl_rational, gsl_rational_mark, gsl_rational_free, rnew);
  } else {
    if (POLY_P(obj)) {
      Data_Get_Struct(obj, gsl_vector, ptmp);
      p = make_vector_clone(ptmp);
    } else {
      p = get_poly_get(obj, &flag);
    }
    for (i = 0; i < p->size; i++) gsl_vector_set(p, i, -gsl_vector_get(p, i));
    return Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, p);
  }
}

static VALUE rb_gsl_rational_uplus(VALUE obj)
{
  return obj;
}

static VALUE rb_gsl_rational_sub(VALUE obj, VALUE other)
{
  return rb_gsl_rational_add(obj, rb_gsl_rational_uminus(other));
}

static VALUE rb_gsl_rational_mul(VALUE obj, VALUE other)
{
  gsl_rational *r = NULL, *r2 = NULL, *rnew = NULL;
  gsl_poly *p;
  int flag = 0;
  Data_Get_Struct(obj, gsl_rational, r);
  if (RATIONAL_P(other)) {
    Data_Get_Struct(other, gsl_rational, r2);
    rnew = gsl_rational_mul(r, r2);
  } else if (VECTOR_P(other)) {
    Data_Get_Struct(other, gsl_vector, p);
    rnew = gsl_rational_mul_poly(r, p);
  } else {
    p = get_poly_get(other, &flag);
    rnew = gsl_rational_mul_poly(r, p);
    gsl_vector_free(p);
  }
  return Data_Wrap_Struct(cgsl_rational, gsl_rational_mark, gsl_rational_free, rnew);
}

static VALUE rb_gsl_rational_div(VALUE obj, VALUE other)
{
  gsl_rational *r = NULL, *r2 = NULL, *rnew = NULL;
  gsl_poly *p;
  size_t i;
  Data_Get_Struct(obj, gsl_rational, r);
  if (RATIONAL_P(other)) {
    Data_Get_Struct(other, gsl_rational, r2);
    rnew = gsl_rational_div(r, r2);
  } else if (VECTOR_P(other)) {
    Data_Get_Struct(other, gsl_vector, p);
    rnew = gsl_rational_div_poly(r, p);
  } else {
    switch (TYPE(other)) {
    case T_ARRAY:
      p = gsl_vector_alloc(RARRAY_LEN(other));
      for (i = 0; i < p->size; i++)
	gsl_vector_set(p, i, NUM2DBL(rb_ary_entry(other, i)));
      rnew = gsl_rational_div_poly(r, p);
      gsl_vector_free(p);
      break;
    case T_FLOAT:
    case T_FIXNUM:
      rnew = gsl_rational_new(r->pnum, r->pden);
      gsl_vector_scale(rnew->pnum, 1.0/NUM2DBL(other));
      break;
    default:
      rb_raise(rb_eTypeError, "wrong argument type %s", 
	       rb_class2name(CLASS_OF(other)));
      break;
    }
  }
  return Data_Wrap_Struct(cgsl_rational, gsl_rational_mark, gsl_rational_free, rnew);
}

static VALUE rb_gsl_rational_inverse(VALUE obj)
{
  gsl_rational *r = NULL, *rnew = NULL;
  Data_Get_Struct(obj, gsl_rational, r);
  rnew = gsl_rational_inverse(r);
  return Data_Wrap_Struct(cgsl_rational, gsl_rational_mark, gsl_rational_free, rnew);
}

static VALUE rb_gsl_rational_den(VALUE obj)
{
  gsl_rational *r = NULL;
  Data_Get_Struct(obj, gsl_rational, r);
  return r->den;
}

static VALUE rb_gsl_rational_num(VALUE obj)
{
  gsl_rational *r = NULL;
  Data_Get_Struct(obj, gsl_rational, r);
  return r->num;
}

static VALUE rb_gsl_rational_coerce(VALUE obj, VALUE other)
{
  gsl_rational *r = NULL;
  gsl_poly *p = NULL, *ptmp = NULL;
  size_t i;
  switch (TYPE(other)) {
  case T_FLOAT:
  case T_FIXNUM:
    p = gsl_vector_alloc(1);
    gsl_vector_set(p, 0, NUM2DBL(other));
    break;
  case T_ARRAY:
    p = gsl_vector_alloc(RARRAY_LEN(other));
    for (i = 0; i < p->size; i++)
      gsl_vector_set(p, i, NUM2DBL(rb_ary_entry(other, i)));
    break;
  default:
    CHECK_VECTOR(other);
    Data_Get_Struct(other, gsl_vector, ptmp);
    p = make_vector_clone(ptmp);
    break;
  }
  ptmp = gsl_vector_alloc(1);
  gsl_vector_set(ptmp, 0, 1.0);
  r = gsl_rational_new2(p, ptmp);
  return rb_ary_new3(2, 
		     Data_Wrap_Struct(cgsl_rational, gsl_rational_mark, gsl_rational_free, r), obj);
}

static VALUE rb_gsl_rational_zero(VALUE obj)
{
  gsl_rational *r = NULL;
  Data_Get_Struct(obj, gsl_rational, r);
  return rb_gsl_poly_complex_solve(1, &(r->num), cgsl_rational);
}

static VALUE rb_gsl_rational_pole(VALUE obj)
{
  gsl_rational *r = NULL;
  Data_Get_Struct(obj, gsl_rational, r);
  return rb_gsl_poly_complex_solve(1, &(r->den), cgsl_rational);
}

static VALUE rb_gsl_poly_inverse(VALUE obj)
{
  gsl_poly *p = NULL, *ptmp = NULL;
  gsl_rational *r = NULL;
  Data_Get_Struct(obj, gsl_poly, p);
  ptmp = gsl_vector_alloc(1);
  gsl_vector_set(ptmp, 0, 1.0);
  r = gsl_rational_new(ptmp, p);
  gsl_vector_free(ptmp);
  return Data_Wrap_Struct(cgsl_rational, gsl_rational_mark, gsl_rational_free, r);
}

void Init_gsl_rational(VALUE module)
{
  cgsl_rational = rb_define_class_under(module, "Rational", cGSL_Object);
  rb_define_singleton_method(cgsl_rational, "new", rb_gsl_rational_new, -1);
  rb_define_singleton_method(cgsl_rational, "[]", rb_gsl_rational_new, -1);
  rb_define_singleton_method(cgsl_rational, "alloc", rb_gsl_rational_new, -1);

  rb_define_method(cgsl_rational, "den", rb_gsl_rational_den, 0);
  rb_define_method(cgsl_rational, "num", rb_gsl_rational_num, 0);

  rb_define_method(cgsl_rational, "print", rb_gsl_rational_print, 0);
  rb_define_method(cgsl_rational, "to_s", rb_gsl_rational_to_s, 0);
  rb_define_method(cgsl_rational, "inspect", rb_gsl_rational_inspect, 0);

  rb_define_method(cgsl_rational, "deconv", rb_gsl_rational_deconv, 0);

  rb_define_method(cgsl_rational, "-@", rb_gsl_rational_uminus, 0);
  rb_define_method(cgsl_rational, "+@", rb_gsl_rational_uplus, 0);

  rb_define_method(cgsl_rational, "add", rb_gsl_rational_add, 1);
  rb_define_alias(cgsl_rational, "+", "add");
  rb_define_method(cgsl_rational, "sub", rb_gsl_rational_sub, 1);
  rb_define_alias(cgsl_rational, "-", "sub");
  rb_define_method(cgsl_rational, "mul", rb_gsl_rational_mul, 1);
  rb_define_alias(cgsl_rational, "*", "mul");
  rb_define_method(cgsl_rational, "div", rb_gsl_rational_div, 1);
  rb_define_alias(cgsl_rational, "/", "div");
  rb_define_method(cgsl_rational, "inverse", rb_gsl_rational_inverse, 0);
  rb_define_alias(cgsl_rational, "inv", "inverse");

  /***** Methods of poly *****/
  rb_define_method(cgsl_poly, "inverse", rb_gsl_poly_inverse, 0);
  rb_define_alias(cgsl_poly, "inv", "inverse");
  rb_define_method(cgsl_poly, "make_rational", rb_gsl_poly_make_rational, 1);
  rb_define_alias(cgsl_poly, "/", "make_rational");

  rb_define_method(cgsl_rational, "coerce", rb_gsl_rational_coerce, 1);

  rb_define_method(cgsl_rational, "zero", rb_gsl_rational_zero, 0);
  rb_define_method(cgsl_rational, "pole", rb_gsl_rational_pole, 0);
}
