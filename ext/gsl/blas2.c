/*
 blas2.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include <gsl/gsl_blas.h>
#include "include/rb_gsl_common.h"
#include "include/rb_gsl_complex.h"
#include "include/rb_gsl_array.h"

static void define_const(VALUE module);

static void define_const(VALUE module)
{
  rb_define_const(module, "CblasRowMajor", INT2FIX(CblasRowMajor));
  rb_define_const(module, "CblasColMajor", INT2FIX(CblasColMajor));
  rb_define_const(module, "RowMajor", INT2FIX(CblasRowMajor));
  rb_define_const(module, "ColMajor", INT2FIX(CblasColMajor));

  rb_define_const(module, "CblasNoTrans", INT2FIX(CblasNoTrans));
  rb_define_const(module, "CblasTrans", INT2FIX(CblasTrans));
  rb_define_const(module, "CblasConjTrans", INT2FIX(CblasConjTrans));
  rb_define_const(module, "NoTrans", INT2FIX(CblasNoTrans));
  rb_define_const(module, "Trans", INT2FIX(CblasTrans));
  rb_define_const(module, "ConjTrans", INT2FIX(CblasConjTrans));

  rb_define_const(module, "CblasUpper", INT2FIX(CblasUpper));
  rb_define_const(module, "CblasLower", INT2FIX(CblasLower));
  rb_define_const(module, "Upper", INT2FIX(CblasUpper));
  rb_define_const(module, "Lower", INT2FIX(CblasLower));

  rb_define_const(module, "CblasNonUnit", INT2FIX(CblasNonUnit));
  rb_define_const(module, "CblasUnit", INT2FIX(CblasUnit));
  rb_define_const(module, "NonUnit", INT2FIX(CblasNonUnit));
  rb_define_const(module, "Unit", INT2FIX(CblasUnit));

  rb_define_const(module, "CblasLeft", INT2FIX(CblasLeft));
  rb_define_const(module, "CblasRight", INT2FIX(CblasRight));
  rb_define_const(module, "Left", INT2FIX(CblasLeft));
  rb_define_const(module, "Right", INT2FIX(CblasRight));
}

static VALUE rb_gsl_blas_dgemv(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *x = NULL, *y = NULL;
  double a, b;
  int type, istart, flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 4)",
         argc);
    CHECK_FIXNUM(argv[0]);
    Need_Float(argv[1]);
    CHECK_MATRIX(argv[2]);
    CHECK_VECTOR(argv[3]);
    type = FIX2INT(argv[0]);
    a = NUM2DBL(argv[1]);
    Data_Get_Struct(argv[2], gsl_matrix, A);
    Data_Get_Struct(argv[3], gsl_vector, x);
    istart = 4;
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix, A);
    if (argc < 3) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 3)",
         argc);
    CHECK_FIXNUM(argv[0]);
    Need_Float(argv[1]);
    CHECK_VECTOR(argv[2]);
    type = FIX2INT(argv[0]);
    a = NUM2DBL(argv[1]);
    Data_Get_Struct(argv[2], gsl_vector, x);
    istart = 3;
    break;
  }
  switch (argc - istart) {
  case 2:
    Need_Float(argv[istart]);
    CHECK_VECTOR(argv[istart+1]);
    b = NUM2DBL(argv[istart]);
    Data_Get_Struct(argv[istart+1], gsl_vector, y);
    break;
  case 0:
    b = 0.0;
    y = gsl_vector_alloc(x->size);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  gsl_blas_dgemv(type, a, A, x, b, y);
  if (flag == 1) return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, y);
  else return argv[argc-1];

}

static VALUE rb_gsl_blas_dgemv2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *x = NULL, *y, *ynew;
  double a, b;
  int type, istart, flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 4)",
         argc);
    CHECK_FIXNUM(argv[0]);
    Need_Float(argv[1]);
    CHECK_MATRIX(argv[2]);
    CHECK_VECTOR(argv[3]);
    type = FIX2INT(argv[0]);
    a = NUM2DBL(argv[1]);
    Data_Get_Struct(argv[2], gsl_matrix, A);
    Data_Get_Struct(argv[3], gsl_vector, x);
    istart = 4;
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix, A);
    if (argc < 3) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 3)",
         argc);
    CHECK_FIXNUM(argv[0]);
    Need_Float(argv[1]);
    CHECK_VECTOR(argv[2]);
    type = FIX2INT(argv[0]);
    a = NUM2DBL(argv[1]);
    Data_Get_Struct(argv[2], gsl_vector, x);
    istart = 3;
    break;
  }
  switch (argc - istart) {
  case 2:
    Need_Float(argv[istart]);
    CHECK_VECTOR(argv[istart+1]);
    b = NUM2DBL(argv[istart]);
    Data_Get_Struct(argv[istart+1], gsl_vector, y);
    break;
  case 0:
    b = 0.0;
    y = gsl_vector_alloc(x->size);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  ynew = gsl_vector_alloc(y->size);
  gsl_vector_memcpy(ynew, y);
  gsl_blas_dgemv(type, a, A, x, b, ynew);
  if (flag == 1) gsl_vector_free(y);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, ynew);
}

static VALUE rb_gsl_blas_zgemv(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *x = NULL, *y = NULL;
  gsl_complex *a, *b, z;
  int type, istart, flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 4)",
         argc);
    CHECK_FIXNUM(argv[0]);
    CHECK_COMPLEX(argv[1]);
    CHECK_MATRIX_COMPLEX(argv[2]);
    CHECK_VECTOR_COMPLEX(argv[3]);
    type = FIX2INT(argv[0]);
    Data_Get_Struct(argv[1], gsl_complex, a);
    Data_Get_Struct(argv[2], gsl_matrix_complex, A);
    Data_Get_Struct(argv[3], gsl_vector_complex, x);
    istart = 4;
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix_complex, A);
    if (argc < 3) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 3)",
         argc);
    CHECK_FIXNUM(argv[0]);
    CHECK_COMPLEX(argv[1]);
    CHECK_VECTOR_COMPLEX(argv[2]);
    type = FIX2INT(argv[0]);
    Data_Get_Struct(argv[1], gsl_complex, a);
    Data_Get_Struct(argv[2], gsl_vector_complex, x);
    istart = 3;
    break;
  }
  switch (argc - istart) {
  case 2:
    CHECK_COMPLEX(argv[istart]);
    CHECK_VECTOR_COMPLEX(argv[istart+1]);
    Data_Get_Struct(argv[istart], gsl_complex, b);
    Data_Get_Struct(argv[istart+1], gsl_vector_complex, y);
    break;
  case 0:
    z = gsl_complex_rect(0.0, 0.0);
    b = &z;
    y = gsl_vector_complex_alloc(x->size);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  gsl_blas_zgemv(type, *a, A, x, *b, y);
  if (flag == 1) return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, y);
  else return argv[argc-1];

}


static VALUE rb_gsl_blas_zgemv2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *x = NULL, *y, *ynew;
  gsl_complex *a, *b, z;
  int type, istart, flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 4)",
         argc);
    CHECK_FIXNUM(argv[0]);
    CHECK_COMPLEX(argv[1]);
    CHECK_MATRIX_COMPLEX(argv[2]);
    CHECK_VECTOR_COMPLEX(argv[3]);
    type = FIX2INT(argv[0]);
    Data_Get_Struct(argv[1], gsl_complex, a);
    Data_Get_Struct(argv[2], gsl_matrix_complex, A);
    Data_Get_Struct(argv[3], gsl_vector_complex, x);
    istart = 4;
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix_complex, A);
    if (argc < 3) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 3)",
         argc);
    CHECK_FIXNUM(argv[0]);
    CHECK_COMPLEX(argv[1]);
    CHECK_VECTOR_COMPLEX(argv[2]);
    type = FIX2INT(argv[0]);
    Data_Get_Struct(argv[1], gsl_complex, a);
    Data_Get_Struct(argv[2], gsl_vector_complex, x);
    istart = 3;
    break;
  }
  switch (argc - istart) {
  case 2:
    CHECK_COMPLEX(argv[istart]);
    CHECK_VECTOR_COMPLEX(argv[istart+1]);
    Data_Get_Struct(argv[istart], gsl_complex, b);
    Data_Get_Struct(argv[istart+1], gsl_vector_complex, y);
    break;
  case 0:
    z = gsl_complex_rect(0.0, 0.0);
    b = &z;
    y = gsl_vector_complex_alloc(x->size);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  ynew = gsl_vector_complex_alloc(y->size);
  gsl_vector_complex_memcpy(ynew, y);
  gsl_blas_zgemv(type, *a, A, x, *b, ynew);
  if (flag == 1) gsl_vector_complex_free(y);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, ynew);
}

static VALUE rb_gsl_blas_dtrmv(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *x = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 5) rb_raise(rb_eArgError, "wrong number of arguments (%d for 5)",
          argc);
    CHECK_MATRIX(argv[3]);
    CHECK_VECTOR(argv[4]);
    Data_Get_Struct(argv[3], gsl_matrix, A);
    Data_Get_Struct(argv[4], gsl_vector, x);
    break;
  default:
    if (argc != 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for 4)",
          argc);
    CHECK_VECTOR(argv[3]);
    Data_Get_Struct(obj, gsl_matrix, A);
    Data_Get_Struct(argv[3], gsl_vector, x);
    break;
  }
  CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
  gsl_blas_dtrmv(FIX2INT(argv[0]), FIX2INT(argv[1]), FIX2INT(argv[2]),
     A, x);
  return argv[argc-1];
}


static VALUE rb_gsl_blas_dtrmv2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *x = NULL, *xnew;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 5) rb_raise(rb_eArgError, "wrong number of arguments (%d for 5)",
          argc);
    CHECK_MATRIX(argv[3]);
    CHECK_VECTOR(argv[4]);
    Data_Get_Struct(argv[3], gsl_matrix, A);
    Data_Get_Struct(argv[4], gsl_vector, x);
    break;
  default:
    if (argc != 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for 4)",
          argc);
    CHECK_VECTOR(argv[3]);
    Data_Get_Struct(obj, gsl_matrix, A);
    Data_Get_Struct(argv[3], gsl_vector, x);
    break;
  }
  CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
  xnew = gsl_vector_alloc(x->size);
  gsl_vector_memcpy(xnew, x);
  gsl_blas_dtrmv(FIX2INT(argv[0]), FIX2INT(argv[1]), FIX2INT(argv[2]),
     A, xnew);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, xnew);
}


static VALUE rb_gsl_blas_ztrmv(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *x = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 5) rb_raise(rb_eArgError, "wrong number of arguments (%d for 5)",
          argc);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_VECTOR_COMPLEX(argv[4]);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_vector_complex, x);
    break;
  default:
    if (argc != 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for 4)",
          argc);
    CHECK_VECTOR_COMPLEX(argv[3]);
    Data_Get_Struct(obj, gsl_matrix_complex, A);
    Data_Get_Struct(argv[3], gsl_vector_complex, x);
    break;
  }
  CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
  gsl_blas_ztrmv(FIX2INT(argv[0]), FIX2INT(argv[1]), FIX2INT(argv[2]),
     A, x);
  return argv[argc-1];
}

static VALUE rb_gsl_blas_ztrmv2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *x = NULL, *xnew;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 5) rb_raise(rb_eArgError, "wrong number of arguments (%d for 5)",
          argc);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_VECTOR_COMPLEX(argv[4]);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_vector_complex, x);
    break;
  default:
    if (argc != 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for 4)",
          argc);
    CHECK_VECTOR_COMPLEX(argv[3]);
    Data_Get_Struct(obj, gsl_matrix_complex, A);
    Data_Get_Struct(argv[3], gsl_vector_complex, x);
    break;
  }
  CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
  xnew = gsl_vector_complex_alloc(x->size);
  gsl_vector_complex_memcpy(xnew, x);
  gsl_blas_ztrmv(FIX2INT(argv[0]), FIX2INT(argv[1]), FIX2INT(argv[2]),
     A, xnew);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, xnew);
}

static VALUE rb_gsl_blas_dtrsv(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *x = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 5) rb_raise(rb_eArgError, "wrong number of arguments (%d for 5)",
          argc);
    CHECK_MATRIX(argv[3]);
    CHECK_VECTOR(argv[4]);
    Data_Get_Struct(argv[3], gsl_matrix, A);
    Data_Get_Struct(argv[4], gsl_vector, x);
    break;
  default:
    if (argc != 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for 4)",
          argc);
    CHECK_VECTOR(argv[3]);
    Data_Get_Struct(obj, gsl_matrix, A);
    Data_Get_Struct(argv[3], gsl_vector, x);
    break;
  }
  CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
  gsl_blas_dtrsv(FIX2INT(argv[0]), FIX2INT(argv[1]), FIX2INT(argv[2]),
     A, x);
  return argv[argc-1];
}

static VALUE rb_gsl_blas_dtrsv2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *x = NULL, *xnew;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 5) rb_raise(rb_eArgError, "wrong number of arguments (%d for 5)",
          argc);
    CHECK_MATRIX(argv[3]);
    CHECK_VECTOR(argv[4]);
    Data_Get_Struct(argv[3], gsl_matrix, A);
    Data_Get_Struct(argv[4], gsl_vector, x);
    break;
  default:
    if (argc != 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for 4)",
          argc);
    CHECK_VECTOR(argv[3]);
    Data_Get_Struct(obj, gsl_matrix, A);
    Data_Get_Struct(argv[3], gsl_vector, x);
    break;
  }
  CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
  xnew = gsl_vector_alloc(x->size);
  gsl_vector_memcpy(xnew, x);
  gsl_blas_dtrsv(FIX2INT(argv[0]), FIX2INT(argv[1]), FIX2INT(argv[2]),
     A, xnew);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, xnew);
}

static VALUE rb_gsl_blas_ztrsv(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *x = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 5) rb_raise(rb_eArgError, "wrong number of arguments (%d for 5)",
          argc);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_VECTOR_COMPLEX(argv[4]);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_vector_complex, x);
    break;
  default:
    if (argc != 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for 4)",
          argc);
    CHECK_VECTOR_COMPLEX(argv[3]);
    Data_Get_Struct(obj, gsl_matrix_complex, A);
    Data_Get_Struct(argv[3], gsl_vector_complex, x);
    break;
  }
  CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
  gsl_blas_ztrsv(FIX2INT(argv[0]), FIX2INT(argv[1]), FIX2INT(argv[2]),
     A, x);
  return argv[argc-1];
}

static VALUE rb_gsl_blas_ztrsv2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *x = NULL, *xnew;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 5) rb_raise(rb_eArgError, "wrong number of arguments (%d for 5)",
          argc);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_VECTOR_COMPLEX(argv[4]);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_vector_complex, x);
    break;
  default:
    if (argc != 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for 4)",
          argc);
    CHECK_VECTOR_COMPLEX(argv[3]);
    Data_Get_Struct(obj, gsl_matrix_complex, A);
    Data_Get_Struct(argv[3], gsl_vector_complex, x);
    break;
  }
  CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
  xnew = gsl_vector_complex_alloc(x->size);
  gsl_vector_complex_memcpy(xnew, x);
  gsl_blas_ztrsv(FIX2INT(argv[0]), FIX2INT(argv[1]), FIX2INT(argv[2]),
     A, xnew);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, xnew);
}

static VALUE rb_gsl_blas_dsymv(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *x = NULL, *y = NULL;
  double a, b;
  int type, istart, flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 4)",
         argc);
    CHECK_FIXNUM(argv[0]);
    Need_Float(argv[1]);
    CHECK_MATRIX(argv[2]);
    CHECK_VECTOR(argv[3]);
    type = FIX2INT(argv[0]);
    a = NUM2DBL(argv[1]);
    Data_Get_Struct(argv[2], gsl_matrix, A);
    Data_Get_Struct(argv[3], gsl_vector, x);
    istart = 4;
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix, A);
    if (argc < 3) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 3)",
         argc);
    CHECK_FIXNUM(argv[0]);
    Need_Float(argv[1]);
    CHECK_VECTOR(argv[2]);
    type = FIX2INT(argv[0]);
    a = NUM2DBL(argv[1]);
    Data_Get_Struct(argv[2], gsl_vector, x);
    istart = 3;
    break;
  }
  switch (argc - istart) {
  case 2:
    Need_Float(argv[istart]);
    CHECK_VECTOR(argv[istart+1]);
    b = NUM2DBL(argv[istart]);
    Data_Get_Struct(argv[istart+1], gsl_vector, y);
    break;
  case 0:
    b = 0.0;
    y = gsl_vector_alloc(x->size);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  gsl_blas_dsymv(type, a, A, x, b, y);
  if (flag == 1) return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, y);
  else return argv[argc-1];

}

static VALUE rb_gsl_blas_dsymv2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *x = NULL, *y, *ynew;
  double a, b;
  int type, istart, flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 4)",
         argc);
    CHECK_FIXNUM(argv[0]);
    Need_Float(argv[1]);
    CHECK_MATRIX(argv[2]);
    CHECK_VECTOR(argv[3]);
    type = FIX2INT(argv[0]);
    a = NUM2DBL(argv[1]);
    Data_Get_Struct(argv[2], gsl_matrix, A);
    Data_Get_Struct(argv[3], gsl_vector, x);
    istart = 4;
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix, A);
    if (argc < 3) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 3)",
         argc);
    CHECK_FIXNUM(argv[0]);
    Need_Float(argv[1]);
    CHECK_VECTOR(argv[2]);
    type = FIX2INT(argv[0]);
    a = NUM2DBL(argv[1]);
    Data_Get_Struct(argv[2], gsl_vector, x);
    istart = 3;
    break;
  }
  switch (argc - istart) {
  case 2:
    Need_Float(argv[istart]);
    CHECK_VECTOR(argv[istart+1]);
    b = NUM2DBL(argv[istart]);
    Data_Get_Struct(argv[istart+1], gsl_vector, y);
    break;
  case 0:
    b = 0.0;
    y = gsl_vector_alloc(x->size);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  ynew = gsl_vector_alloc(y->size);
  gsl_vector_memcpy(ynew, y);
  gsl_blas_dsymv(type, a, A, x, b, ynew);
  if (flag == 1) gsl_vector_free(y);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, ynew);
}

static VALUE rb_gsl_blas_zhemv(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *x = NULL, *y = NULL;
  gsl_complex *a, *b, z;
  int type, istart, flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 4)",
         argc);
    CHECK_FIXNUM(argv[0]);
    CHECK_COMPLEX(argv[1]);
    CHECK_MATRIX_COMPLEX(argv[2]);
    CHECK_VECTOR_COMPLEX(argv[3]);
    type = FIX2INT(argv[0]);
    Data_Get_Struct(argv[1], gsl_complex, a);
    Data_Get_Struct(argv[2], gsl_matrix_complex, A);
    Data_Get_Struct(argv[3], gsl_vector_complex, x);
    istart = 4;
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix_complex, A);
    if (argc < 3) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 3)",
         argc);
    CHECK_FIXNUM(argv[0]);
    CHECK_COMPLEX(argv[1]);
    CHECK_VECTOR_COMPLEX(argv[2]);
    type = FIX2INT(argv[0]);
    Data_Get_Struct(argv[1], gsl_complex, a);
    Data_Get_Struct(argv[2], gsl_vector_complex, x);
    istart = 3;
    break;
  }
  switch (argc - istart) {
  case 2:
    CHECK_COMPLEX(argv[istart]);
    CHECK_VECTOR_COMPLEX(argv[istart+1]);
    Data_Get_Struct(argv[istart], gsl_complex, b);
    Data_Get_Struct(argv[istart+1], gsl_vector_complex, y);
    break;
  case 0:
    z = gsl_complex_rect(0.0, 0.0);
    b = &z;
    y = gsl_vector_complex_alloc(x->size);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  gsl_blas_zhemv(type, *a, A, x, *b, y);
  if (flag == 1) return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, y);
  else return argv[argc-1];

}

static VALUE rb_gsl_blas_zhemv2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *x = NULL, *y, *ynew;
  gsl_complex *a, *b, z;
  int type, istart, flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 4) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 4)",
         argc);
    CHECK_FIXNUM(argv[0]);
    CHECK_COMPLEX(argv[1]);
    CHECK_MATRIX_COMPLEX(argv[2]);
    CHECK_VECTOR_COMPLEX(argv[3]);
    type = FIX2INT(argv[0]);
    Data_Get_Struct(argv[1], gsl_complex, a);
    Data_Get_Struct(argv[2], gsl_matrix_complex, A);
    Data_Get_Struct(argv[3], gsl_vector_complex, x);
    istart = 4;
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix_complex, A);
    if (argc < 3) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 3)",
         argc);
    CHECK_FIXNUM(argv[0]);
    CHECK_COMPLEX(argv[1]);
    CHECK_VECTOR_COMPLEX(argv[2]);
    type = FIX2INT(argv[0]);
    Data_Get_Struct(argv[1], gsl_complex, a);
    Data_Get_Struct(argv[2], gsl_vector_complex, x);
    istart = 3;
    break;
  }
  switch (argc - istart) {
  case 2:
    CHECK_COMPLEX(argv[istart]);
    CHECK_VECTOR_COMPLEX(argv[istart+1]);
    Data_Get_Struct(argv[istart], gsl_complex, b);
    Data_Get_Struct(argv[istart+1], gsl_vector_complex, y);
    break;
  case 0:
    z = gsl_complex_rect(0.0, 0.0);
    b = &z;
    y = gsl_vector_complex_alloc(x->size);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  ynew = gsl_vector_complex_alloc(y->size);
  gsl_vector_complex_memcpy(ynew, y);
  gsl_blas_zhemv(type, *a, A, x, *b, ynew);
  if (flag == 1) gsl_vector_complex_free(y);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, ynew);
}

static VALUE rb_gsl_blas_dger(VALUE obj, VALUE aa, VALUE xx, VALUE yy, VALUE AA)
{
  gsl_matrix *A = NULL;
  gsl_vector *x = NULL, *y = NULL;
  double a;
  Need_Float(aa);
  CHECK_VECTOR(xx);  CHECK_VECTOR(yy);
  CHECK_MATRIX(AA);
  a = NUM2DBL(aa);
  Data_Get_Struct(xx, gsl_vector, x);
  Data_Get_Struct(yy, gsl_vector, y);
  Data_Get_Struct(AA, gsl_matrix, A);
  gsl_blas_dger(a, x, y, A);
  return AA;
}

static VALUE rb_gsl_blas_dger2(VALUE obj, VALUE aa, VALUE xx, VALUE yy, VALUE AA)
{
  gsl_matrix *A = NULL, *Anew = NULL;
  gsl_vector *x = NULL, *y = NULL;
  double a;
  Need_Float(aa);
  CHECK_VECTOR(xx);  CHECK_VECTOR(yy);
  CHECK_MATRIX(AA);
  a = NUM2DBL(aa);
  Data_Get_Struct(xx, gsl_vector, x);
  Data_Get_Struct(yy, gsl_vector, y);
  Data_Get_Struct(AA, gsl_matrix, A);
  Anew = gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix_memcpy(Anew, A);
  gsl_blas_dger(a, x, y, Anew);
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Anew);
}

static VALUE rb_gsl_blas_zgeru(VALUE obj, VALUE aa, VALUE xx, VALUE yy, VALUE AA)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *x = NULL, *y = NULL;
  gsl_complex *a;
  CHECK_COMPLEX(aa);
  CHECK_VECTOR_COMPLEX(xx);CHECK_VECTOR_COMPLEX(yy);
  CHECK_MATRIX_COMPLEX(AA);
  Data_Get_Struct(aa, gsl_complex, a);
  Data_Get_Struct(xx, gsl_vector_complex, x);
  Data_Get_Struct(yy, gsl_vector_complex, y);
  Data_Get_Struct(AA, gsl_matrix_complex, A);
  gsl_blas_zgeru(*a, x, y, A);
  return AA;
}

static VALUE rb_gsl_blas_zgeru2(VALUE obj, VALUE aa, VALUE xx, VALUE yy, VALUE AA)
{
  gsl_matrix_complex *A = NULL, *Anew = NULL;
  gsl_vector_complex *x = NULL, *y = NULL;
  gsl_complex *a;
  CHECK_COMPLEX(aa);
  CHECK_VECTOR_COMPLEX(xx);CHECK_VECTOR_COMPLEX(yy);
  CHECK_MATRIX_COMPLEX(AA);
  Data_Get_Struct(aa, gsl_complex, a);
  Data_Get_Struct(xx, gsl_vector_complex, x);
  Data_Get_Struct(yy, gsl_vector_complex, y);
  Data_Get_Struct(AA, gsl_matrix_complex, A);
  Anew = gsl_matrix_complex_alloc(A->size1, A->size2);
  gsl_matrix_complex_memcpy(Anew, A);
  gsl_blas_zgeru(*a, x, y, Anew);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, Anew);
}

static VALUE rb_gsl_blas_zgerc(VALUE obj, VALUE aa, VALUE xx, VALUE yy, VALUE AA)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *x = NULL, *y = NULL;
  gsl_complex *a;
  CHECK_COMPLEX(aa);
  CHECK_VECTOR_COMPLEX(xx);CHECK_VECTOR_COMPLEX(yy);
  CHECK_MATRIX_COMPLEX(AA);
  Data_Get_Struct(aa, gsl_complex, a);
  Data_Get_Struct(xx, gsl_vector_complex, x);
  Data_Get_Struct(yy, gsl_vector_complex, y);
  Data_Get_Struct(AA, gsl_matrix_complex, A);
  gsl_blas_zgerc(*a, x, y, A);
  return AA;
}

static VALUE rb_gsl_blas_zgerc2(VALUE obj, VALUE aa, VALUE xx, VALUE yy, VALUE AA)
{
  gsl_matrix_complex *A = NULL, *Anew = NULL;
  gsl_vector_complex *x = NULL, *y = NULL;
  gsl_complex *a;
  CHECK_COMPLEX(aa);
  CHECK_VECTOR_COMPLEX(xx);CHECK_VECTOR_COMPLEX(yy);
  CHECK_MATRIX_COMPLEX(AA);
  Data_Get_Struct(aa, gsl_complex, a);
  Data_Get_Struct(xx, gsl_vector_complex, x);
  Data_Get_Struct(yy, gsl_vector_complex, y);
  Data_Get_Struct(AA, gsl_matrix_complex, A);
  Anew = gsl_matrix_complex_alloc(A->size1, A->size2);
  gsl_matrix_complex_memcpy(Anew, A);
  gsl_blas_zgerc(*a, x, y, Anew);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, Anew);
}

static VALUE rb_gsl_blas_dsyr(VALUE obj, VALUE tt, VALUE aa, VALUE xx, VALUE AA)
{
  gsl_matrix *A = NULL;
  gsl_vector *x = NULL;
  double a;
  CHECK_FIXNUM(tt);
  Need_Float(aa);
  CHECK_VECTOR(xx);  CHECK_MATRIX(AA);
  a = NUM2DBL(aa);
  Data_Get_Struct(xx, gsl_vector, x);
  Data_Get_Struct(AA, gsl_matrix, A);
  gsl_blas_dsyr(FIX2INT(tt), a, x, A);
  return AA;
}

static VALUE rb_gsl_blas_dsyr_a(VALUE obj, VALUE tt, VALUE aa, VALUE xx, VALUE AA)
{
  gsl_matrix *A = NULL, *Anew = NULL;
  gsl_vector *x = NULL;
  double a;
  CHECK_FIXNUM(tt);
  Need_Float(aa);
  CHECK_VECTOR(xx);  CHECK_MATRIX(AA);
  a = NUM2DBL(aa);
  Data_Get_Struct(xx, gsl_vector, x);
  Data_Get_Struct(AA, gsl_matrix, A);
  Anew = gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix_memcpy(Anew, A);
  gsl_blas_dsyr(FIX2INT(tt), a, x, Anew);
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Anew);
}

static VALUE rb_gsl_blas_zher(VALUE obj, VALUE tt, VALUE aa, VALUE xx, VALUE AA)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *x = NULL;
  double a;
  CHECK_FIXNUM(tt);
  Need_Float(aa);
  CHECK_VECTOR_COMPLEX(xx);  CHECK_MATRIX_COMPLEX(AA);
  a = NUM2DBL(aa);
  Data_Get_Struct(xx, gsl_vector_complex, x);
  Data_Get_Struct(AA, gsl_matrix_complex, A);
  gsl_blas_zher(FIX2INT(tt), a, x, A);
  return AA;
}

static VALUE rb_gsl_blas_zher_a(VALUE obj, VALUE tt, VALUE aa, VALUE xx, VALUE AA)
{
  gsl_matrix_complex  *A = NULL, *Anew = NULL;
  gsl_vector_complex  *x = NULL;
  double a;
  CHECK_FIXNUM(tt);
  Need_Float(aa);
  CHECK_VECTOR_COMPLEX(xx);  CHECK_MATRIX_COMPLEX(AA);
  a = NUM2DBL(aa);
  Data_Get_Struct(xx, gsl_vector_complex, x);
  Data_Get_Struct(AA, gsl_matrix_complex, A);
  Anew = gsl_matrix_complex_alloc(A->size1, A->size2);
  gsl_matrix_complex_memcpy(Anew, A);
  gsl_blas_zher(FIX2INT(tt), a, x, Anew);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, Anew);
}

static VALUE rb_gsl_blas_dsyr2(VALUE obj, VALUE tt, VALUE aa, VALUE xx,
            VALUE yy, VALUE AA)
{
  gsl_matrix *A = NULL;
  gsl_vector *x = NULL, *y = NULL;
  double a;
  CHECK_FIXNUM(tt);
  Need_Float(aa);
  CHECK_VECTOR(xx); CHECK_VECTOR(yy);  CHECK_MATRIX(AA);
  a = NUM2DBL(aa);
  Data_Get_Struct(xx, gsl_vector, x);
  Data_Get_Struct(yy, gsl_vector, y);
  Data_Get_Struct(AA, gsl_matrix, A);
  gsl_blas_dsyr2(FIX2INT(tt), a, x, y, A);
  return AA;
}

static VALUE rb_gsl_blas_dsyr2_a(VALUE obj, VALUE tt, VALUE aa, VALUE xx,
         VALUE yy, VALUE AA)
{
  gsl_matrix *A = NULL, *Anew = NULL;
  gsl_vector *x = NULL, *y = NULL;
  double a;
  CHECK_FIXNUM(tt);
  Need_Float(aa);
  CHECK_VECTOR(xx); CHECK_VECTOR(yy);  CHECK_MATRIX(AA);
  a = NUM2DBL(aa);
  Data_Get_Struct(xx, gsl_vector, x);
  Data_Get_Struct(yy, gsl_vector, y);
  Data_Get_Struct(AA, gsl_matrix, A);
  Anew = gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix_memcpy(Anew, A);
  gsl_blas_dsyr2(FIX2INT(tt), a, x, y, Anew);
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Anew);
}

static VALUE rb_gsl_blas_zher2(VALUE obj, VALUE tt, VALUE aa, VALUE xx,
            VALUE yy, VALUE AA)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *x = NULL, *y = NULL;
  gsl_complex *a;
  CHECK_FIXNUM(tt);
  CHECK_COMPLEX(aa);
  CHECK_VECTOR_COMPLEX(xx); CHECK_VECTOR_COMPLEX(yy);  CHECK_MATRIX_COMPLEX(AA);
  Data_Get_Struct(aa, gsl_complex, a);
  Data_Get_Struct(xx, gsl_vector_complex, x);
  Data_Get_Struct(yy, gsl_vector_complex, y);
  Data_Get_Struct(AA, gsl_matrix_complex, A);
  gsl_blas_zher2(FIX2INT(tt), *a, x, y, A);
  return AA;
}

static VALUE rb_gsl_blas_zher2_a(VALUE obj, VALUE tt, VALUE aa, VALUE xx,
         VALUE yy, VALUE AA)
{
  gsl_matrix_complex *A = NULL, *Anew = NULL;
  gsl_vector_complex *x = NULL, *y = NULL;
  gsl_complex *a;
  CHECK_FIXNUM(tt);
  CHECK_COMPLEX(aa);
  CHECK_VECTOR_COMPLEX(xx); CHECK_VECTOR_COMPLEX(yy);  CHECK_MATRIX_COMPLEX(AA);
  Data_Get_Struct(aa, gsl_complex, a);
  Data_Get_Struct(xx, gsl_vector_complex, x);
  Data_Get_Struct(yy, gsl_vector_complex, y);
  Data_Get_Struct(AA, gsl_matrix_complex, A);
  Anew = gsl_matrix_complex_alloc(A->size1, A->size2);
  gsl_matrix_complex_memcpy(Anew, A);
  gsl_blas_zher2(FIX2INT(tt), *a, x, y, Anew);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, Anew);
}

void Init_gsl_blas2(VALUE module)
{
  define_const(module);

  rb_define_module_function(module, "dgemv!", rb_gsl_blas_dgemv, -1);
  rb_define_method(cgsl_matrix, "blas_dgemv!", rb_gsl_blas_dgemv, -1);
  rb_define_alias(cgsl_matrix, "dgemv!", "blas_dgemv!");
  rb_define_alias(cgsl_matrix, "gemv!", "blas_dgemv!");

  rb_define_module_function(module, "dgemv", rb_gsl_blas_dgemv2, -1);
  rb_define_method(cgsl_matrix, "blas_dgemv", rb_gsl_blas_dgemv2, -1);
  rb_define_alias(cgsl_matrix, "dgemv", "blas_dgemv");
  rb_define_alias(cgsl_matrix, "gemv", "blas_dgemv");

  rb_define_module_function(module, "zgemv!", rb_gsl_blas_zgemv, -1);
  rb_define_method(cgsl_matrix_complex, "blas_zgemv!", rb_gsl_blas_zgemv, -1);
  rb_define_alias(cgsl_matrix_complex, "zgemv!", "blas_zgemv!");
  rb_define_alias(cgsl_matrix_complex, "gemv!", "blas_zgemv!");

  rb_define_module_function(module, "zgemv", rb_gsl_blas_zgemv2, -1);
  rb_define_method(cgsl_matrix_complex, "blas_zgemv", rb_gsl_blas_zgemv2, -1);
  rb_define_alias(cgsl_matrix_complex, "zgemv", "blas_zgemv");
  rb_define_alias(cgsl_matrix_complex, "gemv", "blas_zgemv");

  rb_define_module_function(module, "dtrmv!", rb_gsl_blas_dtrmv, -1);
  rb_define_method(cgsl_matrix, "blas_dtrmv!", rb_gsl_blas_dtrmv, -1);
  rb_define_alias(cgsl_matrix, "dtrmv!", "blas_dtrmv!");
  rb_define_alias(cgsl_matrix, "trmv!", "blas_dtrmv!");

  rb_define_module_function(module, "dtrmv", rb_gsl_blas_dtrmv2, -1);
  rb_define_method(cgsl_matrix, "blas_dtrmv", rb_gsl_blas_dtrmv2, -1);
  rb_define_alias(cgsl_matrix, "dtrmv", "blas_dtrmv");
  rb_define_alias(cgsl_matrix, "trmv", "blas_dtrmv");

  rb_define_module_function(module, "ztrmv!", rb_gsl_blas_ztrmv, -1);
  rb_define_method(cgsl_matrix_complex, "blas_ztrmv!", rb_gsl_blas_ztrmv, -1);
  rb_define_alias(cgsl_matrix_complex, "ztrmv!", "blas_ztrmv!");

  rb_define_module_function(module, "ztrmv", rb_gsl_blas_ztrmv2, -1);
  rb_define_method(cgsl_matrix_complex, "blas_ztrmv", rb_gsl_blas_ztrmv2, -1);
  rb_define_alias(cgsl_matrix_complex, "ztrmv", "blas_ztrmv");
  rb_define_alias(cgsl_matrix_complex, "trmv", "blas_ztrmv");

  rb_define_module_function(module, "dtrsv!", rb_gsl_blas_dtrsv, -1);
  rb_define_method(cgsl_matrix, "blas_dtrsv!", rb_gsl_blas_dtrsv, -1);
  rb_define_alias(cgsl_matrix, "dtrsv!", "blas_dtrsv!");
  rb_define_alias(cgsl_matrix, "trsv!", "blas_dtrsv!");

  rb_define_module_function(module, "dtrsv", rb_gsl_blas_dtrsv2, -1);
  rb_define_method(cgsl_matrix, "blas_dtrsv", rb_gsl_blas_dtrsv2, -1);
  rb_define_alias(cgsl_matrix, "dtrsv", "blas_dtrsv");
  rb_define_alias(cgsl_matrix, "trsv", "blas_dtrsv");

  rb_define_module_function(module, "ztrsv!", rb_gsl_blas_ztrsv, -1);
  rb_define_method(cgsl_matrix_complex, "blas_ztrsv!", rb_gsl_blas_ztrsv, -1);
  rb_define_alias(cgsl_matrix_complex, "ztrsv!", "blas_ztrsv!");
  rb_define_alias(cgsl_matrix_complex, "trsv!", "blas_ztrsv!");

  rb_define_module_function(module, "ztrsv", rb_gsl_blas_ztrsv2, -1);
  rb_define_method(cgsl_matrix_complex, "blas_ztrsv", rb_gsl_blas_ztrsv2, -1);
  rb_define_alias(cgsl_matrix_complex, "ztrsv", "blas_ztrsv");
  rb_define_alias(cgsl_matrix_complex, "trsv", "blas_ztrsv");

  rb_define_module_function(module, "dsymv!", rb_gsl_blas_dsymv, -1);
  rb_define_method(cgsl_matrix, "blas_dsymv!", rb_gsl_blas_dsymv, -1);
  rb_define_alias(cgsl_matrix, "dsymv!", "blas_dsymv!");
  rb_define_alias(cgsl_matrix, "symv!", "blas_dsymv!");

  rb_define_module_function(module, "dsymv", rb_gsl_blas_dsymv2, -1);
  rb_define_method(cgsl_matrix, "blas_dsymv", rb_gsl_blas_dsymv2, -1);
  rb_define_alias(cgsl_matrix, "dsymv", "blas_dsymv");
  rb_define_alias(cgsl_matrix, "symv", "blas_dsymv");

  rb_define_module_function(module, "zhemv!", rb_gsl_blas_zhemv, -1);
  rb_define_method(cgsl_matrix_complex, "blas_zhemv!", rb_gsl_blas_zhemv, -1);
  rb_define_alias(cgsl_matrix_complex, "zhemv!", "blas_zhemv!");
  rb_define_alias(cgsl_matrix_complex, "symv!", "blas_zhemv!");

  rb_define_module_function(module, "zhemv", rb_gsl_blas_zhemv2, -1);
  rb_define_method(cgsl_matrix_complex, "blas_zhemv", rb_gsl_blas_zhemv2, -1);
  rb_define_alias(cgsl_matrix_complex, "zhemv", "blas_zhemv");
  rb_define_alias(cgsl_matrix_complex, "symv", "blas_zhemv");

  rb_define_module_function(module, "dger!", rb_gsl_blas_dger, 4);
  rb_define_module_function(module, "dger", rb_gsl_blas_dger2, 4);
  rb_define_module_function(module, "zgeru!", rb_gsl_blas_zgeru, 4);
  rb_define_module_function(module, "zgeru", rb_gsl_blas_zgeru2, 4);
  rb_define_module_function(module, "zgerc!", rb_gsl_blas_zgerc, 4);
  rb_define_module_function(module, "zgerc", rb_gsl_blas_zgerc2, 4);
  rb_define_module_function(module, "dsyr!", rb_gsl_blas_dsyr, 4);
  rb_define_module_function(module, "dsyr", rb_gsl_blas_dsyr_a, 4);
  rb_define_module_function(module, "zher!", rb_gsl_blas_zher, 4);
  rb_define_module_function(module, "zher", rb_gsl_blas_zher_a, 4);

  rb_define_module_function(module, "dsyr2!", rb_gsl_blas_dsyr2, 4);
  rb_define_module_function(module, "dsyr2", rb_gsl_blas_dsyr2_a, 4);
  rb_define_module_function(module, "zher2!", rb_gsl_blas_zher2, 4);
  rb_define_module_function(module, "zher2", rb_gsl_blas_zher2_a, 4);
}
