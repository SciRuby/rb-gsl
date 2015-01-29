/*
  blas1.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include <gsl/gsl_blas.h>
#include "include/rb_gsl_complex.h"
#include "include/rb_gsl_array.h"
#include "include/rb_gsl_common.h"

static int get_vector1(int argc, VALUE *argv, VALUE obj, gsl_vector **x);
static int get_vector_complex1(int argc, VALUE *argv, VALUE obj, gsl_vector_complex **x);
static int get_vector2(int argc, VALUE *argv, VALUE obj,
           gsl_vector **x, gsl_vector **y);
static int get_vector_complex2(int argc, VALUE *argv, VALUE obj,
           gsl_vector_complex **x, gsl_vector_complex **y);

static int get_vector1(int argc, VALUE *argv, VALUE obj, gsl_vector **x)
{
  int flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    Data_Get_Vector(argv[0], (*x));
    break;
  default:
    Data_Get_Vector(obj, (*x));
    flag = 1;
    break;
  }
  return flag;
}

static int get_vector_complex1(int argc, VALUE *argv, VALUE obj, gsl_vector_complex **x)
{
  int flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    CHECK_VECTOR_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector_complex, (*x));
    break;
  default:
    Data_Get_Struct(obj, gsl_vector_complex, (*x));
    flag = 1;
    break;
  }
  return flag;
}

static int get_vector2(int argc, VALUE *argv, VALUE obj,
           gsl_vector **x, gsl_vector **y)
{
  int flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    Data_Get_Vector(argv[0], (*x));
    Data_Get_Vector(argv[1], (*y));
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    Data_Get_Vector(obj, (*x));
    Data_Get_Vector(argv[0], (*y));
    flag = 1;
    break;
  }
  return flag;
}


static int get_vector_complex2(int argc, VALUE *argv, VALUE obj,
           gsl_vector_complex **x, gsl_vector_complex **y)
{
  int flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    CHECK_VECTOR_COMPLEX(argv[0]);
    CHECK_VECTOR_COMPLEX(argv[1]);
    Data_Get_Struct(argv[0], gsl_vector_complex, (*x));
    Data_Get_Struct(argv[1], gsl_vector_complex, (*y));
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    CHECK_VECTOR_COMPLEX(argv[0]);
    Data_Get_Struct(obj, gsl_vector_complex, (*x));
    Data_Get_Struct(argv[0], gsl_vector_complex, (*y));
    flag = 1;
    break;
  }
  return flag;
}

static VALUE rb_gsl_blas_ddot(int argc, VALUE *argv, VALUE obj)
{
  double r;
  // local variable "status" declared and set, but never used
  //int status;
  gsl_vector *x = NULL, *y = NULL;
  get_vector2(argc, argv, obj, &x, &y);
  /*status =*/ gsl_blas_ddot(x, y, &r);
  return rb_float_new(r);
}

static VALUE rb_gsl_blas_zdotu(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex *r;
  // local variable "status" declared and set, but never used
  //int status;
  gsl_vector_complex *x = NULL, *y = NULL;
  get_vector_complex2(argc, argv, obj, &x, &y);
  r = ALLOC(gsl_complex);
  /*status =*/ gsl_blas_zdotu(x, y, r);
  return Data_Wrap_Struct(cgsl_complex, 0, free, r);
}

static VALUE rb_gsl_blas_zdotc(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex *r;
  // local variable "status" declared and set, but never used
  //int status;
  gsl_vector_complex *x = NULL, *y = NULL;
  get_vector_complex2(argc, argv, obj, &x, &y);
  r = ALLOC(gsl_complex);
  /*status =*/ gsl_blas_zdotc(x, y, r);
  return Data_Wrap_Struct(cgsl_complex, 0, free, r);
}

static VALUE rb_gsl_blas_dnrm2(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *x = NULL;
  get_vector1(argc, argv, obj, &x);
  return rb_float_new(gsl_blas_dnrm2(x));
}

static VALUE rb_gsl_blas_dnrm(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *x = NULL;
  double a;
  get_vector1(argc, argv, obj, &x);
  a = gsl_blas_dnrm2(x);
  return rb_float_new(a*a);
}

static VALUE rb_gsl_blas_dznrm2(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *x = NULL;
  get_vector_complex1(argc, argv, obj, &x);
  return rb_float_new(gsl_blas_dznrm2(x));
}

static VALUE rb_gsl_blas_dasum(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *x = NULL;
  get_vector1(argc, argv, obj, &x);
  return rb_float_new(gsl_blas_dasum(x));
}

static VALUE rb_gsl_blas_dzasum(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *x = NULL;
  get_vector_complex1(argc, argv, obj, &x);
  return rb_float_new(gsl_blas_dzasum(x));
}

static VALUE rb_gsl_blas_idamax(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *x = NULL;
  get_vector1(argc, argv, obj, &x);
  return INT2FIX(gsl_blas_idamax(x));
}

static VALUE rb_gsl_blas_izamax(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *x = NULL;
  get_vector_complex1(argc, argv, obj, &x);
  return INT2FIX(gsl_blas_izamax(x));
}

static VALUE rb_gsl_blas_dswap(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *x = NULL, *y = NULL;
  get_vector2(argc, argv, obj, &x, &y);
  return INT2FIX(gsl_blas_dswap(x, y));
}

static VALUE rb_gsl_blas_zswap(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *x = NULL, *y = NULL;
  get_vector_complex2(argc, argv, obj, &x, &y);
  return INT2FIX(gsl_blas_zswap(x, y));
}

static VALUE rb_gsl_blas_dcopy(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *x = NULL, *y = NULL;
  get_vector2(argc, argv, obj, &x, &y);
  return INT2FIX(gsl_blas_dcopy(x, y));
}

static VALUE rb_gsl_blas_zcopy(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *x = NULL, *y = NULL;
  get_vector_complex2(argc, argv, obj, &x, &y);
  return INT2FIX(gsl_blas_zcopy(x, y));
}

static VALUE rb_gsl_blas_daxpy(int argc, VALUE *argv, VALUE obj)
{
  double a;
  gsl_vector *x = NULL, *y = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    get_vector2(argc-1, argv+1, obj, &x, &y);
    Need_Float(argv[0]);
    //    a = RFLOAT(argv[0])->value;
    a = NUM2DBL(argv[0]);
    break;
  default:
    Data_Get_Struct(obj, gsl_vector, x);
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    Need_Float(argv[0]);
    //    a = RFLOAT(argv[0])->value;
    a = NUM2DBL(argv[0]);
    Data_Get_Vector(argv[1], y);
    break;
  }
  gsl_blas_daxpy(a, x, y);
  return argv[argc-1];
}

static VALUE rb_gsl_blas_daxpy2(int argc, VALUE *argv, VALUE obj)
{
  double a;
  gsl_vector *x = NULL, *y = NULL, *y2 = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    get_vector2(argc-1, argv+1, obj, &x, &y);
    Need_Float(argv[0]);
    //    a = RFLOAT(argv[0])->value;
    a = NUM2DBL(argv[0]);
    break;
  default:
    Data_Get_Struct(obj, gsl_vector, x);
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    Need_Float(argv[0]);
    CHECK_VECTOR(argv[1]);
    //    a = RFLOAT(argv[0])->value;
    a = NUM2DBL(argv[0]);
    Data_Get_Struct(argv[1], gsl_vector, y);
    break;
  }
  y2 = gsl_vector_alloc(y->size);
  gsl_vector_memcpy(y2, y);
  gsl_blas_daxpy(a, x, y2);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, y2);
}

static VALUE rb_gsl_blas_zaxpy(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex *a = NULL;
  gsl_vector_complex *x = NULL, *y = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    CHECK_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_complex, a);
    get_vector_complex2(argc-1, argv+1, obj, &x, &y);
    break;
  default:
    Data_Get_Struct(obj, gsl_vector_complex, x);
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    CHECK_COMPLEX(argv[0]);
    CHECK_VECTOR_COMPLEX(argv[1]);
    Data_Get_Struct(argv[0], gsl_complex, a);
    Data_Get_Struct(argv[1], gsl_vector_complex, y);
    break;
  }

  gsl_blas_zaxpy(*a, x, y);
  return argv[argc-1];
}

static VALUE rb_gsl_blas_zaxpy2(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex *a = NULL;
  gsl_vector_complex *x = NULL, *y = NULL, *y2 = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    get_vector_complex2(argc-1, argv+1, obj, &x, &y);
    CHECK_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_complex, a);
    break;
  default:
    Data_Get_Struct(obj, gsl_vector_complex, x);
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    CHECK_COMPLEX(argv[0]);
    CHECK_VECTOR_COMPLEX(argv[1]);
    Data_Get_Struct(argv[0], gsl_complex, a);
    Data_Get_Struct(argv[1], gsl_vector_complex, y);
    break;
  }
  y2 = gsl_vector_complex_alloc(y->size);
  gsl_vector_complex_memcpy(y2, y);
  gsl_blas_zaxpy(*a, x, y2);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, y2);
}

static VALUE rb_gsl_blas_dscal(int argc, VALUE *argv, VALUE obj)
{
  double a;
  gsl_vector *x = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    Need_Float(argv[0]);
    CHECK_VECTOR(argv[1]);
    //    a = RFLOAT(argv[0])->value;
    a = NUM2DBL(argv[0]);
    Data_Get_Struct(argv[1], gsl_vector, x);
    gsl_blas_dscal(a, x);
    return argv[1];
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    Need_Float(argv[0]);
    //    a = RFLOAT(argv[0])->value;
    a = NUM2DBL(argv[0]);
    Data_Get_Struct(obj, gsl_vector, x);
    gsl_blas_dscal(a, x);
    return obj;
    break;
  }
  return Qnil; /* never reach here */
}

static VALUE rb_gsl_blas_dscal2(int argc, VALUE *argv, VALUE obj)
{
  double a;
  gsl_vector *x = NULL, *xnew = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    Need_Float(argv[0]);
    CHECK_VECTOR(argv[1]);
    a = NUM2DBL(argv[0]);
    Data_Get_Struct(argv[1], gsl_vector, x);
    break;
  default:
    Data_Get_Struct(obj, gsl_vector, x);
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    Need_Float(argv[0]);
    a = NUM2DBL(argv[0]);
    break;
  }
  xnew = gsl_vector_alloc(x->size);
  gsl_vector_memcpy(xnew, x);
  gsl_blas_dscal(a, xnew);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, xnew);
}

static VALUE rb_gsl_blas_zdscal(int argc, VALUE *argv, VALUE obj)
{
  double a;
  gsl_vector_complex *x = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    Need_Float(argv[0]);
    CHECK_VECTOR_COMPLEX(argv[1]);
    //    a = RFLOAT(argv[0])->value;
    a = NUM2DBL(argv[0]);
    Data_Get_Struct(argv[1], gsl_vector_complex, x);
    gsl_blas_zdscal(a, x);
    return argv[1];
    break;
  default:
    Data_Get_Struct(obj, gsl_vector_complex, x);
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    Need_Float(argv[0]);
    a = NUM2DBL(argv[0]);
    gsl_blas_zdscal(a, x);
    return obj;
    break;
  }
}

static VALUE rb_gsl_blas_zdscal2(int argc, VALUE *argv, VALUE obj)
{
  double a;
  gsl_vector_complex *x = NULL, *xnew = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    Need_Float(argv[0]);
    CHECK_VECTOR_COMPLEX(argv[1]);
    a = NUM2DBL(argv[0]);
    Data_Get_Struct(argv[1], gsl_vector_complex, x);
    break;
  default:
    Data_Get_Struct(obj, gsl_vector_complex, x);
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    Need_Float(argv[0]);
    a = NUM2DBL(argv[0]);
    break;
  }
  xnew = gsl_vector_complex_alloc(x->size);
  gsl_vector_complex_memcpy(xnew, x);
  gsl_blas_zdscal(a, xnew);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, xnew);
}

static VALUE rb_gsl_blas_zscal(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex *a = NULL;
  gsl_vector_complex *x = NULL;
  CHECK_COMPLEX(argv[0]);
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    CHECK_VECTOR_COMPLEX(argv[1]);
    Data_Get_Struct(argv[0], gsl_complex, a);
    Data_Get_Struct(argv[1], gsl_vector_complex, x);
    gsl_blas_zscal(*a, x);
    return argv[1];
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    Data_Get_Struct(obj, gsl_vector_complex, x);
    Data_Get_Struct(argv[0], gsl_complex, a);
    gsl_blas_zscal(*a, x);
    return obj;
    break;
  }
}

static VALUE rb_gsl_blas_zscal2(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex *a = NULL;
  gsl_vector_complex *x = NULL, *xnew = NULL;
  CHECK_COMPLEX(argv[0]);
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    CHECK_VECTOR_COMPLEX(argv[1]);
    Data_Get_Struct(argv[0], gsl_complex, a);
    Data_Get_Struct(argv[1], gsl_vector_complex, x);
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    Data_Get_Struct(obj, gsl_vector_complex, x);
    Data_Get_Struct(argv[0], gsl_complex, a);
    break;
  }
  xnew = gsl_vector_complex_alloc(x->size);
  gsl_vector_complex_memcpy(xnew, x);
  gsl_blas_zscal(*a, xnew);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, xnew);
}

static VALUE rb_gsl_blas_drot(VALUE obj, VALUE xx, VALUE yy, VALUE cc, VALUE ss)
{
  gsl_vector *x = NULL, *y = NULL;
  double c, s;
  CHECK_VECTOR(xx);
  CHECK_VECTOR(yy);
  Need_Float(cc);
  Need_Float(ss);
  Data_Get_Struct(xx, gsl_vector, x);
  Data_Get_Struct(yy, gsl_vector, y);
  c = NUM2DBL(cc);
  s = NUM2DBL(ss);
  gsl_blas_drot(x, y, c, s);
  return rb_ary_new3(2, xx, yy);
}

static VALUE rb_gsl_blas_drot2(VALUE obj, VALUE xx, VALUE yy, VALUE cc, VALUE ss)
{
  gsl_vector *x = NULL, *y = NULL, *xnew = NULL, *ynew = NULL;
  double c, s;
  CHECK_VECTOR(xx);
  CHECK_VECTOR(yy);
  Need_Float(cc);
  Need_Float(ss);
  Data_Get_Struct(xx, gsl_vector, x);
  Data_Get_Struct(yy, gsl_vector, y);
  c = NUM2DBL(cc);
  s = NUM2DBL(ss);
  xnew = gsl_vector_alloc(x->size);
  ynew = gsl_vector_alloc(y->size);
  gsl_vector_memcpy(xnew, x);
  gsl_vector_memcpy(ynew, y);
  gsl_blas_drot(xnew, ynew, c, s);
  return rb_ary_new3(2, Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, xnew),
         Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, ynew));
}

static VALUE rb_gsl_blas_drotm(VALUE obj, VALUE xx, VALUE yy, VALUE PP)
{
  gsl_vector *x = NULL, *y = NULL, *p = NULL;
  int flag = 0, i;
  CHECK_VECTOR(xx);
  CHECK_VECTOR(yy);
  Data_Get_Struct(xx, gsl_vector, x);
  Data_Get_Struct(yy, gsl_vector, y);
  if (rb_obj_is_kind_of(PP, cgsl_vector)) {
    Data_Get_Struct(PP, gsl_vector, p);
  } else {
    if (TYPE(PP) != T_ARRAY) rb_raise(rb_eTypeError, "wrong argument type %s (Array of Vector expected", rb_class2name(CLASS_OF(PP)));
    //    p = gsl_vector_alloc(RARRAY(PP)->len);
    p = gsl_vector_alloc(RARRAY_LEN(PP));
    for (i = 0; i < RARRAY_LEN(PP); i++) {
      gsl_vector_set(p, i, rb_ary_entry(PP, i));
    }
    flag = 1;
  }
  gsl_blas_drotm(x, y, p->data);
  if (flag == 1) gsl_vector_free(p);
  return rb_ary_new3(2, xx, yy);
}

static VALUE rb_gsl_blas_drotm2(VALUE obj, VALUE xx, VALUE yy, VALUE PP)
{
  gsl_vector *x = NULL, *y = NULL, *p = NULL, *xnew = NULL, *ynew = NULL;
  int flag = 0, i;
  CHECK_VECTOR(xx);
  CHECK_VECTOR(yy);
  Data_Get_Struct(xx, gsl_vector, x);
  Data_Get_Struct(yy, gsl_vector, y);
  if (rb_obj_is_kind_of(PP, cgsl_vector)) {
    Data_Get_Struct(PP, gsl_vector, p);
  } else {
    if (TYPE(PP) != T_ARRAY) rb_raise(rb_eTypeError, "wrong argument type %s (Array of Vector expected", rb_class2name(CLASS_OF(PP)));
    //    p = gsl_vector_alloc(RARRAY(PP)->len);
    p = gsl_vector_alloc(RARRAY_LEN(PP));
    for (i = 0; i < RARRAY_LEN(PP); i++) {
      gsl_vector_set(p, i, rb_ary_entry(PP, i));
    }
    flag = 1;
  }
  xnew = gsl_vector_alloc(x->size);
  ynew = gsl_vector_alloc(y->size);
  gsl_vector_memcpy(xnew, x);
  gsl_vector_memcpy(ynew, y);
  gsl_blas_drotm(xnew, ynew, p->data);
  if (flag == 1) gsl_vector_free(p);
  return rb_ary_new3(2, Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, xnew),
         Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, ynew));
}

void Init_gsl_blas1(VALUE module)
{
  rb_define_module_function(module, "ddot", rb_gsl_blas_ddot, -1);
  rb_define_method(cgsl_vector, "blas_ddot", rb_gsl_blas_ddot, -1);
  rb_define_alias(cgsl_vector, "ddot", "blas_ddot");
  /*  rb_define_alias(cgsl_vector, "dot", "blas_ddot");*/

  rb_define_module_function(module, "zdotu", rb_gsl_blas_zdotu, -1);
  rb_define_method(cgsl_vector_complex, "blas_zdotu", rb_gsl_blas_zdotu, -1);
  rb_define_alias(cgsl_vector_complex, "zdotu", "blas_zdotu");
  rb_define_alias(cgsl_vector_complex, "dotu", "blas_zdotu");

  rb_define_module_function(module, "zdotc", rb_gsl_blas_zdotc, -1);
  rb_define_method(cgsl_vector_complex, "blas_zdotc", rb_gsl_blas_zdotc, -1);
  rb_define_alias(cgsl_vector_complex, "zdotc", "blas_zdotc");
  rb_define_alias(cgsl_vector_complex, "dotc", "blas_zdotc");

  rb_define_module_function(module, "dnrm2", rb_gsl_blas_dnrm2, -1);
  rb_define_method(cgsl_vector, "blas_dnrm2", rb_gsl_blas_dnrm2, -1);
  rb_define_alias(cgsl_vector, "dnrm2", "blas_dnrm2");
  rb_define_alias(cgsl_vector, "nrm2", "blas_dnrm2");
  rb_define_alias(cgsl_vector, "norm", "blas_dnrm2");

  rb_define_module_function(module, "dnrm", rb_gsl_blas_dnrm, -1);
  rb_define_method(cgsl_vector, "blas_dnrm", rb_gsl_blas_dnrm, -1);
  rb_define_alias(cgsl_vector, "dnrm", "blas_dnrm");
  rb_define_alias(cgsl_vector, "nrm", "blas_dnrm");
  rb_define_alias(cgsl_vector, "sumsq", "blas_dnrm");

  rb_define_module_function(module, "dznrm2", rb_gsl_blas_dznrm2, -1);
  rb_define_method(cgsl_vector_complex, "blas_dznrm2", rb_gsl_blas_dznrm2, -1);
  rb_define_alias(cgsl_vector_complex, "dznrm2", "blas_dznrm2");
  rb_define_alias(cgsl_vector_complex, "nrm2", "blas_dznrm2");

  rb_define_module_function(module, "dasum", rb_gsl_blas_dasum, -1);
  rb_define_method(cgsl_vector, "blas_dasum", rb_gsl_blas_dasum, -1);
  rb_define_alias(cgsl_vector, "dasum", "blas_dasum");
  rb_define_alias(cgsl_vector, "asum", "blas_dasum");

  rb_define_module_function(module, "dzasum", rb_gsl_blas_dzasum, -1);
  rb_define_method(cgsl_vector_complex, "blas_dzasum", rb_gsl_blas_dzasum, -1);
  rb_define_alias(cgsl_vector_complex, "dzasum", "blas_dzasum");
  rb_define_alias(cgsl_vector_complex, "asum", "blas_dzasum");

  rb_define_module_function(module, "idamax", rb_gsl_blas_idamax, -1);
  rb_define_method(cgsl_vector, "blas_idamax", rb_gsl_blas_idamax, -1);
  rb_define_alias(cgsl_vector, "idamax", "blas_idamax");

  rb_define_module_function(module, "izamax", rb_gsl_blas_izamax, -1);
  rb_define_method(cgsl_vector_complex, "blas_izamax", rb_gsl_blas_izamax, -1);
  rb_define_alias(cgsl_vector_complex, "izamax", "blas_izamax");

  rb_define_module_function(module, "dswap", rb_gsl_blas_dswap, -1);
  rb_define_method(cgsl_vector, "blas_dswap", rb_gsl_blas_dswap, -1);
  rb_define_alias(cgsl_vector, "dswap", "blas_dswap");
  rb_define_alias(cgsl_vector, "swap", "blas_dswap");

  rb_define_module_function(module, "zswap", rb_gsl_blas_zswap, -1);
  rb_define_method(cgsl_vector_complex, "blas_zswap", rb_gsl_blas_zswap, -1);
  rb_define_alias(cgsl_vector_complex, "zswap", "blas_zswap");
  rb_define_alias(cgsl_vector_complex, "swap", "blas_zswap");

  rb_define_module_function(module, "dcopy", rb_gsl_blas_dcopy, -1);
  rb_define_method(cgsl_vector, "blas_dcopy", rb_gsl_blas_dcopy, -1);
  rb_define_alias(cgsl_vector, "dcopy", "blas_dcopy");
  rb_define_alias(cgsl_vector, "copy", "blas_dcopy");

  rb_define_module_function(module, "zcopy", rb_gsl_blas_zcopy, -1);
  rb_define_method(cgsl_vector_complex, "blas_zcopy", rb_gsl_blas_zcopy, -1);
  rb_define_alias(cgsl_vector_complex, "zcopy", "blas_zcopy");
  rb_define_alias(cgsl_vector_complex, "copy", "blas_zcopy");

  rb_define_module_function(module, "daxpy!", rb_gsl_blas_daxpy, -1);
  rb_define_method(cgsl_vector, "blas_daxpy!", rb_gsl_blas_daxpy, -1);
  rb_define_alias(cgsl_vector, "daxpy!", "blas_daxpy!");
  rb_define_alias(cgsl_vector, "axpy!", "blas_daxpy!");

  rb_define_module_function(module, "daxpy", rb_gsl_blas_daxpy2, -1);
  rb_define_method(cgsl_vector, "blas_daxpy", rb_gsl_blas_daxpy2, -1);
  rb_define_alias(cgsl_vector, "daxpy", "blas_daxpy");
  rb_define_alias(cgsl_vector, "axpy", "blas_daxpy");

  rb_define_module_function(module, "zaxpy!", rb_gsl_blas_zaxpy, -1);
  rb_define_method(cgsl_vector_complex, "blas_zaxpy!", rb_gsl_blas_zaxpy, -1);
  rb_define_alias(cgsl_vector_complex, "zaxpy!", "blas_zaxpy!");
  rb_define_alias(cgsl_vector_complex, "axpy!", "blas_zaxpy!");

  rb_define_module_function(module, "zaxpy", rb_gsl_blas_zaxpy2, -1);
  rb_define_method(cgsl_vector_complex, "blas_zaxpy", rb_gsl_blas_zaxpy2, -1);
  rb_define_alias(cgsl_vector_complex, "zaxpy", "blas_zaxpy");
  rb_define_alias(cgsl_vector_complex, "axpy", "blas_zaxpy");

  rb_define_module_function(module, "dscal!", rb_gsl_blas_dscal, -1);
  rb_define_method(cgsl_vector, "blas_dscal!", rb_gsl_blas_dscal, -1);
  rb_define_alias(cgsl_vector, "dscal!", "blas_dscal!");
  rb_define_alias(cgsl_vector, "scal!", "blas_dscal!");

  rb_define_module_function(module, "dscal", rb_gsl_blas_dscal2, -1);
  rb_define_method(cgsl_vector, "blas_dscal", rb_gsl_blas_dscal2, -1);
  rb_define_alias(cgsl_vector, "dscal", "blas_dscal");
  rb_define_alias(cgsl_vector, "scal", "blas_dscal");

  rb_define_module_function(module, "zdscal!", rb_gsl_blas_zdscal, -1);
  rb_define_method(cgsl_vector_complex, "blas_zdscal!", rb_gsl_blas_zdscal, -1);
  rb_define_alias(cgsl_vector_complex, "zdscal!", "blas_zdscal!");
  rb_define_alias(cgsl_vector_complex, "scal!", "blas_zdscal!");

  rb_define_module_function(module, "zdscal", rb_gsl_blas_zdscal2, -1);
  rb_define_method(cgsl_vector_complex, "blas_zdscal", rb_gsl_blas_zdscal2, -1);
  rb_define_alias(cgsl_vector_complex, "zdscal", "blas_zdscal");
  rb_define_alias(cgsl_vector_complex, "scal", "blas_zdscal");

  rb_define_module_function(module, "zscal!", rb_gsl_blas_zscal, -1);
  rb_define_method(cgsl_vector_complex, "blas_zscal!", rb_gsl_blas_zscal, -1);
  rb_define_alias(cgsl_vector_complex, "zscal!", "blas_zscal!");

  rb_define_module_function(module, "zscal2", rb_gsl_blas_zscal2, -1);
  rb_define_method(cgsl_vector_complex, "blas_zscal2", rb_gsl_blas_zscal2, -1);
  rb_define_alias(cgsl_vector_complex, "zscal2", "blas_zscal2");
  rb_define_alias(cgsl_vector_complex, "scal2", "blas_zscal2");

  rb_define_module_function(module, "drot!", rb_gsl_blas_drot, 4);
  rb_define_module_function(module, "drot", rb_gsl_blas_drot2, 4);

  rb_define_module_function(module, "drotm!", rb_gsl_blas_drotm, 3);
  rb_define_module_function(module, "drotm", rb_gsl_blas_drotm2, 3);
}
