/*
  qrng.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_rng.h"
#include "rb_gsl_array.h"
#include <gsl/gsl_qrng.h>
#ifdef HAVE_QRNGEXTRA_QRNGEXTRA_H
#include <qrngextra/qrngextra.h>
#endif

enum rb_gsl_qrng_generator {
  GSL_QRNG_NIEDERREITER_2,
  GSL_QRNG_SOBOL,
  GSL_QRNG_HALTON,
  GSL_QRNG_REVERSEHALTON,
#ifdef HAVE_QRNGEXTRA_QRNGEXTRA_H
  GSL_QRNG_HDSOBOL,
#endif
};

static const gsl_qrng_type* get_gsl_qrng_type(VALUE t);

static const gsl_qrng_type* get_gsl_qrng_type(VALUE t)
{
  const gsl_qrng_type *T;
  char name[32];

  switch (TYPE(t)) {
  case T_STRING:
    strcpy(name, STR2CSTR(t));
    if (strstr(name, "niederreiter_2")) return T = gsl_qrng_niederreiter_2;
#ifdef HAVE_QRNGEXTRA_QRNGEXTRA_H
	else if (strstr(name, "hdsobol")) return T = qrngextra_hdsobol;
#endif
    else if (strstr(name, "sobol")) return T = gsl_qrng_sobol;
#ifdef GSL_1_11_LATER
    else if (strstr(name, "reversehalton")) return T = gsl_qrng_reversehalton;
    else if (strstr(name, "halton")) return T = gsl_qrng_halton;

#endif
    else rb_raise(rb_eArgError, "unknown type");
    break;
  case T_FIXNUM:
    switch (FIX2INT(t)) {
    case GSL_QRNG_NIEDERREITER_2: T = gsl_qrng_niederreiter_2; break;
    case GSL_QRNG_SOBOL: T = gsl_qrng_sobol; break;
#ifdef GSL_1_11_LATER
    case GSL_QRNG_HALTON: T = gsl_qrng_halton; break;
    case GSL_QRNG_REVERSEHALTON: T = gsl_qrng_reversehalton; break;
#endif
#ifdef HAVE_QRNGEXTRA_QRNGEXTRA_H
	case GSL_QRNG_HDSOBOL: T = qrngextra_hdsobol; break;
#endif
    default: 
      rb_raise(rb_eArgError, "unknown type");
    }
    break;
  default: 
    rb_raise(rb_eTypeError, "wrong argument type %s (String or Fixnum expected)",
	     rb_class2name(CLASS_OF(t)));
  }
  return T;
}

static VALUE rb_gsl_qrng_new(VALUE klass, VALUE t, VALUE dd)
{
  unsigned int d;
  gsl_qrng *q;
  const gsl_qrng_type *T;
  d = NUM2UINT(dd);
  T = get_gsl_qrng_type(t);
  q = gsl_qrng_alloc(T, d);
  return Data_Wrap_Struct(klass, 0, gsl_qrng_free, q);
}

#ifdef GSL_1_2_LATER
static VALUE rb_gsl_qrng_init(VALUE obj)
{
  gsl_qrng *q = NULL;
  Data_Get_Struct(obj, gsl_qrng, q);
  gsl_qrng_init(q);
  return obj;
}
#endif

static VALUE rb_gsl_qrng_name(VALUE obj)
{
  gsl_qrng *q = NULL;
  Data_Get_Struct(obj, gsl_qrng, q);
  return rb_str_new2(gsl_qrng_name(q));
}

static VALUE rb_gsl_qrng_size(VALUE obj)
{
  gsl_qrng *q = NULL;
  Data_Get_Struct(obj, gsl_qrng, q);
  return INT2FIX(gsl_qrng_size(q));
}

static VALUE rb_gsl_qrng_clone(VALUE obj)
{
  gsl_qrng *q = NULL, *q2 = NULL;
  Data_Get_Struct(obj, gsl_qrng, q);
  q2 = gsl_qrng_clone(q);
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_qrng_free, q2);
}

/* singleton */
static VALUE rb_gsl_qrng_memcpy(VALUE obj, VALUE dest, VALUE src)
{
  gsl_qrng *q = NULL, *q2 = NULL;
  Data_Get_Struct(dest, gsl_qrng, q);
  Data_Get_Struct(src, gsl_qrng, q2);
  gsl_qrng_memcpy(q, q2);
  return dest;
}

static VALUE rb_gsl_qrng_get(int argc, VALUE *argv, VALUE obj)
{
  gsl_qrng *q = NULL;
  gsl_vector *v;
  Data_Get_Struct(obj, gsl_qrng, q);
  if (argc == 0) {
    v = gsl_vector_alloc(q->dimension);
    gsl_qrng_get(q, v->data);
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
  } else {
    if (!rb_obj_is_kind_of(argv[0], cgsl_vector)) {
      rb_raise(rb_eArgError, "wrong type argument (GSL_Vector required)");
    }
    Data_Get_Struct(argv[0], gsl_vector, v);
    return INT2FIX(gsl_qrng_get(q, v->data));
  }
}

void Init_gsl_qrng(VALUE module)
{
  VALUE cgsl_qrng;
  cgsl_qrng = rb_define_class_under(module, "QRng", cGSL_Object);

  rb_define_singleton_method(cgsl_qrng, "new", rb_gsl_qrng_new, 2);
  rb_define_singleton_method(cgsl_qrng, "alloc", rb_gsl_qrng_new, 2);
#ifdef GSL_1_2_LATER
  rb_define_method(cgsl_qrng, "init", rb_gsl_qrng_init, 0);
#endif
  rb_define_method(cgsl_qrng, "name", rb_gsl_qrng_name, 0);
  rb_define_method(cgsl_qrng, "size", rb_gsl_qrng_size, 0);
  rb_define_method(cgsl_qrng, "clone", rb_gsl_qrng_clone, 0);
  rb_define_alias(cgsl_qrng, "duplicate", "clone");
  rb_define_singleton_method(cgsl_qrng, "memcpy", rb_gsl_qrng_memcpy, 2);

  rb_define_method(cgsl_qrng, "get", rb_gsl_qrng_get, -1);

  rb_define_const(cgsl_qrng, "NIEDERREITER_2", INT2FIX(GSL_QRNG_NIEDERREITER_2));
  rb_define_const(cgsl_qrng, "SOBOL", INT2FIX(GSL_QRNG_SOBOL));
#ifdef GSL_1_11_LATER
  rb_define_const(cgsl_qrng, "HALTON", INT2FIX(GSL_QRNG_HALTON));
  rb_define_const(cgsl_qrng, "REVERSEHALTON", INT2FIX(GSL_QRNG_REVERSEHALTON));
#endif
#ifdef HAVE_QRNGEXTRA_QRNGEXTRA_H
    rb_define_const(cgsl_qrng, "HDSOBOL", INT2FIX(GSL_QRNG_HDSOBOL));
#endif
}
