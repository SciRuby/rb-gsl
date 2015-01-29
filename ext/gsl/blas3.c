/*
  blas3.c
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

static VALUE rb_gsl_blas_dgemm(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *B = NULL, *C = NULL;
  double alpha, beta;
  CBLAS_TRANSPOSE_t TransA, TransB;
  int flag = 0;
  switch (argc) {
  case 2:
    CHECK_MATRIX(argv[0]);
    CHECK_MATRIX(argv[1]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    Data_Get_Struct(argv[1], gsl_matrix, B);
    C = gsl_matrix_calloc(A->size1, B->size2);
    alpha = 1.0;
    beta = 0.0;
    TransA = CblasNoTrans;  TransB = CblasNoTrans;
    flag = 1;
    break;
  case 5:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    Need_Float(argv[2]);
    CHECK_MATRIX(argv[3]);
    CHECK_MATRIX(argv[4]);
    TransA = FIX2INT(argv[0]);
    TransB = FIX2INT(argv[1]);
    alpha = NUM2DBL(argv[2]);
    Data_Get_Struct(argv[3], gsl_matrix, A);
    Data_Get_Struct(argv[4], gsl_matrix, B);
    C = gsl_matrix_calloc(A->size1, B->size2);
    beta = 0.0;
    flag = 1;
    break;
  case 6:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    Need_Float(argv[2]);
    CHECK_MATRIX(argv[3]);
    CHECK_MATRIX(argv[4]);
    Need_Float(argv[5]);
    TransA = FIX2INT(argv[0]);
    TransB = FIX2INT(argv[1]);
    alpha = NUM2DBL(argv[2]);
    Data_Get_Struct(argv[3], gsl_matrix, A);
    Data_Get_Struct(argv[4], gsl_matrix, B);
    beta = NUM2DBL(argv[5]);
    C = gsl_matrix_calloc(A->size1, B->size2);
    flag = 1;
    break;
  case 7:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    Need_Float(argv[2]);
    CHECK_MATRIX(argv[3]);
    CHECK_MATRIX(argv[4]);
    Need_Float(argv[5]);
    CHECK_MATRIX(argv[6]);
    TransA = FIX2INT(argv[0]);
    TransB = FIX2INT(argv[1]);
    alpha = NUM2DBL(argv[2]);
    Data_Get_Struct(argv[3], gsl_matrix, A);
    Data_Get_Struct(argv[4], gsl_matrix, B);
    beta = NUM2DBL(argv[5]);
    Data_Get_Struct(argv[6], gsl_matrix, C);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2, 5, 6, or 7)", argc);
    break;
  }
  gsl_blas_dgemm(TransA, TransB, alpha, A, B, beta, C);
  if (flag == 1) return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, C);
  else return argv[6];
}

static VALUE rb_gsl_blas_zgemm(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL, *B = NULL, *C = NULL;
  gsl_complex alpha, beta, *pa = &alpha, *pb = &beta;
  CBLAS_TRANSPOSE_t TransA, TransB;
  int flag = 0;
  alpha.dat[0] = 1.0; alpha.dat[1] = 0.0;
  beta.dat[0] = 0.0; beta.dat[1] = 0.0;
  switch (argc) {
  case 2:
    CHECK_MATRIX_COMPLEX(argv[0]);
    CHECK_MATRIX_COMPLEX(argv[1]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, A);
    Data_Get_Struct(argv[1], gsl_matrix_complex, B);
    C = gsl_matrix_complex_calloc(A->size1, B->size2);
    TransA = CblasNoTrans;  TransB = CblasNoTrans;
    flag = 1;
    break;
  case 5:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    CHECK_COMPLEX(argv[2]);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_MATRIX_COMPLEX(argv[4]);
    TransA = FIX2INT(argv[0]);
    TransB = FIX2INT(argv[1]);
    Data_Get_Struct(argv[2], gsl_complex, pa);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_matrix_complex, B);
    C = gsl_matrix_complex_calloc(A->size1, B->size2);
    flag = 1;
    break;
  case 6:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    CHECK_COMPLEX(argv[2]);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_MATRIX_COMPLEX(argv[4]);
    CHECK_COMPLEX(argv[5]);
    TransA = FIX2INT(argv[0]);
    TransB = FIX2INT(argv[1]);
    Data_Get_Struct(argv[2], gsl_complex, pa);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_matrix_complex, B);
    Data_Get_Struct(argv[5], gsl_complex, pb);
    C = gsl_matrix_complex_calloc(A->size1, B->size2);
    flag = 1;
    break;
  case 7:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    CHECK_COMPLEX(argv[2]);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_MATRIX_COMPLEX(argv[4]);
    CHECK_COMPLEX(argv[5]);
    CHECK_MATRIX_COMPLEX(argv[6]);
    TransA = FIX2INT(argv[0]);
    TransB = FIX2INT(argv[1]);
    Data_Get_Struct(argv[2], gsl_complex, pa);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_matrix_complex, B);
    Data_Get_Struct(argv[5], gsl_complex, pb);
    Data_Get_Struct(argv[6], gsl_matrix_complex, C);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 7)", argc);
    break;
  }
  gsl_blas_zgemm(TransA, TransB, alpha, A, B, beta, C);
  if (flag == 1) return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, C);
  else return argv[6];
}

static VALUE rb_gsl_blas_dsymm(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *B = NULL, *C = NULL;
  double alpha, beta;
  CBLAS_SIDE_t Side;
  CBLAS_UPLO_t Uplo;
  int flag = 0;
  switch (argc) {
  case 2:
    CHECK_MATRIX(argv[0]);
    CHECK_MATRIX(argv[1]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    Data_Get_Struct(argv[1], gsl_matrix, B);
    C = gsl_matrix_calloc(A->size1, B->size2);
    alpha = 1.0;
    beta = 0.0;
    Side = CblasLeft;  Uplo = CblasUpper;
    flag = 1;
    break;
  case 5:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    Need_Float(argv[2]);
    CHECK_MATRIX(argv[3]);
    CHECK_MATRIX(argv[4]);
    Side = FIX2INT(argv[0]);
    Uplo = FIX2INT(argv[1]);
    alpha = NUM2DBL(argv[2]);
    Data_Get_Struct(argv[3], gsl_matrix, A);
    Data_Get_Struct(argv[4], gsl_matrix, B);
    C = gsl_matrix_calloc(A->size1, B->size2);
    beta = 0.0;
    flag = 1;
    break;
  case 6:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    Need_Float(argv[2]);
    CHECK_MATRIX(argv[3]);
    CHECK_MATRIX(argv[4]);
    Need_Float(argv[5]);
    CHECK_MATRIX(argv[6]);
    Side = FIX2INT(argv[0]);
    Uplo = FIX2INT(argv[1]);
    alpha = NUM2DBL(argv[2]);
    Data_Get_Struct(argv[3], gsl_matrix, A);
    Data_Get_Struct(argv[4], gsl_matrix, B);
    beta = NUM2DBL(argv[5]);
    C = gsl_matrix_calloc(A->size1, B->size2);
    flag = 1;
    break;
  case 7:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    Need_Float(argv[2]);
    CHECK_MATRIX(argv[3]);
    CHECK_MATRIX(argv[4]);
    Need_Float(argv[5]);
    CHECK_MATRIX(argv[6]);
    Side = FIX2INT(argv[0]);
    Uplo = FIX2INT(argv[1]);
    alpha = NUM2DBL(argv[2]);
    Data_Get_Struct(argv[3], gsl_matrix, A);
    Data_Get_Struct(argv[4], gsl_matrix, B);
    beta = NUM2DBL(argv[5]);
    Data_Get_Struct(argv[6], gsl_matrix, C);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 7)", argc);
    break;
  }
  gsl_blas_dsymm(Side, Uplo, alpha, A, B, beta, C);
  if (flag == 1) return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, C);
  else return argv[6];
}

static VALUE rb_gsl_blas_zsymm(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL, *B = NULL, *C = NULL;
  gsl_complex alpha, beta, *pa = &alpha, *pb = &beta;
  CBLAS_SIDE_t Side;
  CBLAS_UPLO_t Uplo;
  int flag = 0;
  alpha = gsl_complex_rect(1.0, 0.0);
  beta = gsl_complex_rect(0.0, 0.0);
  switch (argc) {
  case 2:
    CHECK_MATRIX_COMPLEX(argv[0]);
    CHECK_MATRIX_COMPLEX(argv[1]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, A);
    Data_Get_Struct(argv[1], gsl_matrix_complex, B);
    C = gsl_matrix_complex_calloc(A->size1, B->size2);
    Side = CblasLeft;  Uplo = CblasUpper;
    flag = 1;
    break;
  case 5:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    CHECK_COMPLEX(argv[2]);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_MATRIX_COMPLEX(argv[4]);
    Side = FIX2INT(argv[0]);
    Uplo = FIX2INT(argv[1]);
    Data_Get_Struct(argv[2], gsl_complex, pa);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_matrix_complex, B);
    C = gsl_matrix_complex_calloc(A->size1, B->size2);
    flag = 1;
    break;
  case 6:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    CHECK_COMPLEX(argv[2]);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_MATRIX_COMPLEX(argv[4]);
    CHECK_COMPLEX(argv[5]);
    CHECK_MATRIX_COMPLEX(argv[6]);
    Side = FIX2INT(argv[0]);
    Uplo = FIX2INT(argv[1]);
    Data_Get_Struct(argv[2], gsl_complex, pa);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_matrix_complex, B);
    Data_Get_Struct(argv[5], gsl_complex, pb);
    C = gsl_matrix_complex_calloc(A->size1, B->size2);
    flag = 1;
    break;
  case 7:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    CHECK_COMPLEX(argv[2]);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_MATRIX_COMPLEX(argv[4]);
    CHECK_COMPLEX(argv[5]);
    CHECK_MATRIX_COMPLEX(argv[6]);
    Side = FIX2INT(argv[0]);
    Uplo = FIX2INT(argv[1]);
    Data_Get_Struct(argv[2], gsl_complex, pa);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_matrix_complex, B);
    Data_Get_Struct(argv[5], gsl_complex, pb);
    Data_Get_Struct(argv[6], gsl_matrix_complex, C);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 7)", argc);
    break;
  }
  gsl_blas_zsymm(Side, Uplo, alpha, A, B, beta, C);
  if (flag == 1) return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, C);
  else return argv[6];
}

static VALUE rb_gsl_blas_zhemm(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL, *B = NULL, *C = NULL;
  gsl_complex alpha, beta, *pa = &alpha, *pb = &beta;
  CBLAS_SIDE_t Side;
  CBLAS_UPLO_t Uplo;
  int flag = 0;
  alpha = gsl_complex_rect(1.0, 0.0);
  beta = gsl_complex_rect(0.0, 0.0);
  switch (argc) {
  case 2:
    CHECK_MATRIX_COMPLEX(argv[0]);
    CHECK_MATRIX_COMPLEX(argv[1]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, A);
    Data_Get_Struct(argv[1], gsl_matrix_complex, B);
    C = gsl_matrix_complex_calloc(A->size1, B->size2);
    Side = CblasLeft;  Uplo = CblasUpper;
    flag = 1;
    break;
  case 5:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    CHECK_COMPLEX(argv[2]);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_MATRIX_COMPLEX(argv[4]);
    Side = FIX2INT(argv[0]);
    Uplo = FIX2INT(argv[1]);
    Data_Get_Struct(argv[2], gsl_complex, pa);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_matrix_complex, B);
    C = gsl_matrix_complex_calloc(A->size1, B->size2);
    flag = 1;
    break;
  case 6:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    CHECK_COMPLEX(argv[2]);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_MATRIX_COMPLEX(argv[4]);
    CHECK_COMPLEX(argv[5]);
    CHECK_MATRIX_COMPLEX(argv[6]);
    Side = FIX2INT(argv[0]);
    Uplo = FIX2INT(argv[1]);
    Data_Get_Struct(argv[2], gsl_complex, pa);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_matrix_complex, B);
    Data_Get_Struct(argv[5], gsl_complex, pb);
    C = gsl_matrix_complex_calloc(A->size1, B->size2);
    flag = 1;
    break;
  case 7:
    CHECK_FIXNUM(argv[0]);
    CHECK_FIXNUM(argv[1]);
    CHECK_COMPLEX(argv[2]);
    CHECK_MATRIX_COMPLEX(argv[3]);
    CHECK_MATRIX_COMPLEX(argv[4]);
    CHECK_COMPLEX(argv[5]);
    CHECK_MATRIX_COMPLEX(argv[6]);
    Side = FIX2INT(argv[0]);
    Uplo = FIX2INT(argv[1]);
    Data_Get_Struct(argv[2], gsl_complex, pa);
    Data_Get_Struct(argv[3], gsl_matrix_complex, A);
    Data_Get_Struct(argv[4], gsl_matrix_complex, B);
    Data_Get_Struct(argv[5], gsl_complex, pb);
    Data_Get_Struct(argv[6], gsl_matrix_complex, C);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 7)", argc);
    break;
  }
  gsl_blas_zhemm(Side, Uplo, alpha, A, B, beta, C);
  if (flag == 1) return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, C);
  else return argv[6];
}

static VALUE rb_gsl_blas_dtrmm(VALUE obj, VALUE s, VALUE u, VALUE ta,
             VALUE d, VALUE a, VALUE aa, VALUE bb)
{
  gsl_matrix *A = NULL, *B = NULL;
  double alpha;
  CBLAS_SIDE_t Side;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t TransA;
  CBLAS_DIAG_t Diag;
  CHECK_FIXNUM(s);  CHECK_FIXNUM(u);  CHECK_FIXNUM(ta);  CHECK_FIXNUM(d);
  Need_Float(a);
  CHECK_MATRIX(aa);  CHECK_MATRIX(bb);
  Side = FIX2INT(s);
  Uplo = FIX2INT(u);
  TransA = FIX2INT(ta);
  Diag = FIX2INT(d);
  alpha = NUM2DBL(a);
  Data_Get_Struct(aa, gsl_matrix, A);
  Data_Get_Struct(bb, gsl_matrix, B);
  gsl_blas_dtrmm(Side, Uplo, TransA, Diag, alpha, A, B);
  return bb;
}

static VALUE rb_gsl_blas_dtrmm2(VALUE obj, VALUE s, VALUE u, VALUE ta,
             VALUE d, VALUE a, VALUE aa, VALUE bb)
{
  gsl_matrix *A = NULL, *B = NULL, *Bnew = NULL;
  double alpha;
  CBLAS_SIDE_t Side;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t TransA;
  CBLAS_DIAG_t Diag;
  CHECK_FIXNUM(s);  CHECK_FIXNUM(u);  CHECK_FIXNUM(ta);  CHECK_FIXNUM(d);
  Need_Float(a);
  CHECK_MATRIX(aa);  CHECK_MATRIX(bb);
  Side = FIX2INT(s);
  Uplo = FIX2INT(u);
  TransA = FIX2INT(ta);
  Diag = FIX2INT(d);
  alpha = NUM2DBL(a);
  Data_Get_Struct(aa, gsl_matrix, A);
  Data_Get_Struct(bb, gsl_matrix, B);
  Bnew = gsl_matrix_alloc(B->size1, B->size2);
  gsl_matrix_memcpy(Bnew, B);
  gsl_blas_dtrmm(Side, Uplo, TransA, Diag, alpha, A, Bnew);
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Bnew);
}

static VALUE rb_gsl_blas_ztrmm(VALUE obj, VALUE s, VALUE u, VALUE ta,
             VALUE d, VALUE a, VALUE aa, VALUE bb)
{
  gsl_matrix_complex *A = NULL, *B = NULL;
  gsl_complex alpha, *pa = &alpha;
  CBLAS_SIDE_t Side;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t TransA;
  CBLAS_DIAG_t Diag;
  CHECK_FIXNUM(s);  CHECK_FIXNUM(u);  CHECK_FIXNUM(ta);  CHECK_FIXNUM(d);
  CHECK_COMPLEX(a);
  CHECK_MATRIX_COMPLEX(aa);
  CHECK_MATRIX_COMPLEX(bb);
  Side = FIX2INT(s);
  Uplo = FIX2INT(u);
  TransA = FIX2INT(ta);
  Diag = FIX2INT(d);
  Data_Get_Struct(a, gsl_complex, pa);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  Data_Get_Struct(bb, gsl_matrix_complex, B);
  gsl_blas_ztrmm(Side, Uplo, TransA, Diag, alpha, A, B);
  return bb;
}

static VALUE rb_gsl_blas_ztrmm2(VALUE obj, VALUE s, VALUE u, VALUE ta,
             VALUE d, VALUE a, VALUE aa, VALUE bb)
{
  gsl_matrix_complex *A = NULL, *B = NULL, *Bnew = NULL;
  gsl_complex alpha, *pa = &alpha;
  CBLAS_SIDE_t Side;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t TransA;
  CBLAS_DIAG_t Diag;
  CHECK_FIXNUM(s);  CHECK_FIXNUM(u);  CHECK_FIXNUM(ta);  CHECK_FIXNUM(d);
  CHECK_COMPLEX(a);
  CHECK_MATRIX_COMPLEX(aa);  CHECK_MATRIX_COMPLEX(bb);
  Side = FIX2INT(s);
  Uplo = FIX2INT(u);
  TransA = FIX2INT(ta);
  Diag = FIX2INT(d);
  Data_Get_Struct(a, gsl_complex, pa);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  Data_Get_Struct(bb, gsl_matrix_complex, B);
  Bnew = gsl_matrix_complex_alloc(B->size1, B->size2);
  gsl_matrix_complex_memcpy(Bnew, B);
  gsl_blas_ztrmm(Side, Uplo, TransA, Diag, alpha, A, Bnew);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, Bnew);
}

static VALUE rb_gsl_blas_dtrsm(VALUE obj, VALUE s, VALUE u, VALUE ta,
             VALUE d, VALUE a, VALUE aa, VALUE bb)
{
  gsl_matrix *A = NULL, *B = NULL;
  double alpha;
  CBLAS_SIDE_t Side;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t TransA;
  CBLAS_DIAG_t Diag;
  CHECK_FIXNUM(s);  CHECK_FIXNUM(u);  CHECK_FIXNUM(ta);  CHECK_FIXNUM(d);
  Need_Float(a);
  CHECK_MATRIX(aa);  CHECK_MATRIX(bb);
  Side = FIX2INT(s);
  Uplo = FIX2INT(u);
  TransA = FIX2INT(ta);
  Diag = FIX2INT(d);
  alpha = NUM2DBL(a);
  Data_Get_Struct(aa, gsl_matrix, A);
  Data_Get_Struct(bb, gsl_matrix, B);
  gsl_blas_dtrsm(Side, Uplo, TransA, Diag, alpha, A, B);
  return bb;
}

static VALUE rb_gsl_blas_dtrsm2(VALUE obj, VALUE s, VALUE u, VALUE ta,
             VALUE d, VALUE a, VALUE aa, VALUE bb)
{
  gsl_matrix *A = NULL, *B = NULL, *Bnew = NULL;
  double alpha;
  CBLAS_SIDE_t Side;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t TransA;
  CBLAS_DIAG_t Diag;
  CHECK_FIXNUM(s);  CHECK_FIXNUM(u);  CHECK_FIXNUM(ta);  CHECK_FIXNUM(d);
  Need_Float(a);
  CHECK_MATRIX(aa);  CHECK_MATRIX(bb);
  Side = FIX2INT(s);
  Uplo = FIX2INT(u);
  TransA = FIX2INT(ta);
  Diag = FIX2INT(d);
  alpha = NUM2DBL(a);
  Data_Get_Struct(aa, gsl_matrix, A);
  Data_Get_Struct(bb, gsl_matrix, B);
  Bnew = gsl_matrix_alloc(B->size1, B->size2);
  gsl_matrix_memcpy(Bnew, B);
  gsl_blas_dtrsm(Side, Uplo, TransA, Diag, alpha, A, Bnew);
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Bnew);
}

static VALUE rb_gsl_blas_ztrsm(VALUE obj, VALUE s, VALUE u, VALUE ta,
             VALUE d, VALUE a, VALUE aa, VALUE bb)
{
  gsl_matrix_complex *A = NULL, *B = NULL;
  gsl_complex alpha, *pa = &alpha;
  CBLAS_SIDE_t Side;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t TransA;
  CBLAS_DIAG_t Diag;
  CHECK_FIXNUM(s);  CHECK_FIXNUM(u);  CHECK_FIXNUM(ta);  CHECK_FIXNUM(d);
  CHECK_COMPLEX(a);
  CHECK_MATRIX_COMPLEX(aa);
  CHECK_MATRIX_COMPLEX(bb);
  Side = FIX2INT(s);
  Uplo = FIX2INT(u);
  TransA = FIX2INT(ta);
  Diag = FIX2INT(d);
  Data_Get_Struct(a, gsl_complex, pa);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  Data_Get_Struct(bb, gsl_matrix_complex, B);
  gsl_blas_ztrsm(Side, Uplo, TransA, Diag, alpha, A, B);
  return bb;
}

static VALUE rb_gsl_blas_ztrsm2(VALUE obj, VALUE s, VALUE u, VALUE ta,
             VALUE d, VALUE a, VALUE aa, VALUE bb)
{
  gsl_matrix_complex *A = NULL, *B = NULL, *Bnew = NULL;
  gsl_complex alpha, *pa = &alpha;
  CBLAS_SIDE_t Side;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t TransA;
  CBLAS_DIAG_t Diag;
  CHECK_FIXNUM(s);  CHECK_FIXNUM(u);  CHECK_FIXNUM(ta);  CHECK_FIXNUM(d);
  CHECK_COMPLEX(a);
  CHECK_MATRIX_COMPLEX(aa);
  CHECK_MATRIX_COMPLEX(bb);
  Side = FIX2INT(s);
  Uplo = FIX2INT(u);
  TransA = FIX2INT(ta);
  Diag = FIX2INT(d);
  Data_Get_Struct(a, gsl_complex, pa);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  Data_Get_Struct(bb, gsl_matrix_complex, B);
  Bnew = gsl_matrix_complex_alloc(B->size1, B->size2);
  gsl_matrix_complex_memcpy(Bnew, B);
  gsl_blas_ztrsm(Side, Uplo, TransA, Diag, alpha, A, Bnew);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, Bnew);
}

static VALUE rb_gsl_blas_dsyrk(VALUE obj, VALUE u, VALUE t, VALUE a, VALUE aa,
             VALUE b, VALUE cc)
{
  gsl_matrix *A = NULL, *C = NULL;
  double alpha, beta;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t Trans;
  CHECK_FIXNUM(u);  CHECK_FIXNUM(t);
  Need_Float(a);  Need_Float(b);
  CHECK_MATRIX(aa);  CHECK_MATRIX(cc);
  Uplo = FIX2INT(u);
  Trans = FIX2INT(t);
  alpha = NUM2DBL(a);
  beta = NUM2DBL(b);
  Data_Get_Struct(aa, gsl_matrix, A);
  Data_Get_Struct(cc, gsl_matrix, C);
  gsl_blas_dsyrk(Uplo, Trans, alpha, A, beta, C);
  return cc;
}

static VALUE rb_gsl_blas_dsyrk2(VALUE obj, VALUE u, VALUE t, VALUE a, VALUE aa,
             VALUE b, VALUE cc)
{
  gsl_matrix *A = NULL, *C = NULL, *Cnew = NULL;
  double alpha, beta;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t Trans;
  CHECK_FIXNUM(u);  CHECK_FIXNUM(t);
  Need_Float(a);  Need_Float(b);
  CHECK_MATRIX(aa);  CHECK_MATRIX(cc);
  Uplo = FIX2INT(u);
  Trans = FIX2INT(t);
  alpha = NUM2DBL(a);
  beta = NUM2DBL(b);
  Data_Get_Struct(aa, gsl_matrix, A);
  Data_Get_Struct(cc, gsl_matrix, C);
  Cnew = gsl_matrix_alloc(C->size1, C->size2);
  gsl_matrix_memcpy(Cnew, C);
  gsl_blas_dsyrk(Uplo, Trans, alpha, A, beta, Cnew);
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Cnew);
}

static VALUE rb_gsl_blas_zsyrk(VALUE obj, VALUE u, VALUE t, VALUE a, VALUE aa,
             VALUE b, VALUE cc)
{
  gsl_matrix_complex *A = NULL, *C = NULL;
  gsl_complex alpha, beta, *pa = &alpha, *pb = &beta;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t Trans;
  CHECK_FIXNUM(u);  CHECK_FIXNUM(t);
  CHECK_COMPLEX(a);  CHECK_COMPLEX(b);
  CHECK_MATRIX_COMPLEX(aa);  CHECK_MATRIX_COMPLEX(cc);
  Uplo = FIX2INT(u);
  Trans = FIX2INT(t);
  Data_Get_Struct(a, gsl_complex, pa);
  Data_Get_Struct(b, gsl_complex, pb);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  Data_Get_Struct(cc, gsl_matrix_complex, C);
  gsl_blas_zsyrk(Uplo, Trans, alpha, A, beta, C);
  return cc;
}

static VALUE rb_gsl_blas_zsyrk2(VALUE obj, VALUE u, VALUE t, VALUE a, VALUE aa,
             VALUE b, VALUE cc)
{
  gsl_matrix_complex *A = NULL, *C = NULL, *Cnew = NULL;
  gsl_complex alpha, beta, *pa = &alpha, *pb = &beta;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t Trans;
  CHECK_FIXNUM(u);  CHECK_FIXNUM(t);
  CHECK_COMPLEX(a);  CHECK_COMPLEX(b);
  CHECK_MATRIX_COMPLEX(aa);  CHECK_MATRIX_COMPLEX(cc);
  Uplo = FIX2INT(u);
  Trans = FIX2INT(t);
  Data_Get_Struct(a, gsl_complex, pa);
  Data_Get_Struct(b, gsl_complex, pb);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  Data_Get_Struct(cc, gsl_matrix_complex, C);
  Cnew = gsl_matrix_complex_alloc(C->size1, C->size2);
  gsl_matrix_complex_memcpy(Cnew, C);
  gsl_blas_zsyrk(Uplo, Trans, alpha, A, beta, Cnew);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, Cnew);
}

static VALUE rb_gsl_blas_zherk(VALUE obj, VALUE u, VALUE t, VALUE a, VALUE aa,
             VALUE b, VALUE cc)
{
  gsl_matrix_complex *A = NULL, *C = NULL;
  double alpha, beta;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t Trans;
  CHECK_FIXNUM(u);  CHECK_FIXNUM(t);
  Need_Float(a);  Need_Float(b);
  CHECK_MATRIX_COMPLEX(aa);
  CHECK_MATRIX_COMPLEX(cc);
  Uplo = FIX2INT(u);
  Trans = FIX2INT(t);
  alpha = NUM2DBL(a);
  beta = NUM2DBL(b);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  Data_Get_Struct(cc, gsl_matrix_complex, C);
  gsl_blas_zherk(Uplo, Trans, alpha, A, beta, C);
  return cc;
}

static VALUE rb_gsl_blas_zherk2(VALUE obj, VALUE u, VALUE t, VALUE a, VALUE aa,
             VALUE b, VALUE cc)
{
  gsl_matrix_complex *A = NULL, *C = NULL, *Cnew = NULL;
  double alpha, beta;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t Trans;
  CHECK_FIXNUM(u);  CHECK_FIXNUM(t);
  Need_Float(a);  Need_Float(b);
  CHECK_MATRIX_COMPLEX(aa);  CHECK_MATRIX_COMPLEX(cc);
  Uplo = FIX2INT(u);
  Trans = FIX2INT(t);
  alpha = NUM2DBL(a);
  beta = NUM2DBL(b);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  Data_Get_Struct(cc, gsl_matrix_complex, C);
  Cnew = gsl_matrix_complex_alloc(C->size1, C->size2);
  gsl_matrix_complex_memcpy(Cnew, C);
  gsl_blas_zherk(Uplo, Trans, alpha, A, beta, Cnew);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, Cnew);
}

static VALUE rb_gsl_blas_dsyr2k(VALUE obj, VALUE u, VALUE t, VALUE a, VALUE aa,
        VALUE bb, VALUE b, VALUE cc)
{
  gsl_matrix *A = NULL, *B = NULL, *C = NULL;
  double alpha, beta;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t Trans;
  CHECK_FIXNUM(u);  CHECK_FIXNUM(t);
  Need_Float(a);  Need_Float(b);
  CHECK_MATRIX(aa);  CHECK_MATRIX(bb);  CHECK_MATRIX(cc);
  Uplo = FIX2INT(u);
  Trans = FIX2INT(t);
  alpha = NUM2DBL(a);
  Data_Get_Struct(aa, gsl_matrix, A);
  Data_Get_Struct(bb, gsl_matrix, B);
  beta = NUM2DBL(b);
  Data_Get_Struct(cc, gsl_matrix, C);
  gsl_blas_dsyr2k(Uplo, Trans, alpha, A, B, beta, C);
  return cc;
}

static VALUE rb_gsl_blas_dsyr2k2(VALUE obj, VALUE u, VALUE t, VALUE a, VALUE aa,
        VALUE bb, VALUE b, VALUE cc)
{
  gsl_matrix *A = NULL, *B = NULL, *C = NULL, *Cnew = NULL;
  double alpha, beta;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t Trans;
  CHECK_FIXNUM(u);  CHECK_FIXNUM(t);
  Need_Float(a);  Need_Float(b);
  CHECK_MATRIX(aa);  CHECK_MATRIX(bb);  CHECK_MATRIX(cc);
  Uplo = FIX2INT(u);
  Trans = FIX2INT(t);
  alpha = NUM2DBL(a);
  Data_Get_Struct(aa, gsl_matrix, A);
  Data_Get_Struct(bb, gsl_matrix, B);
  beta = NUM2DBL(b);
  Data_Get_Struct(cc, gsl_matrix, C);
  Cnew = gsl_matrix_alloc(C->size1, C->size2);
  gsl_matrix_memcpy(Cnew, C);
  gsl_blas_dsyr2k(Uplo, Trans, alpha, A, B, beta, Cnew);
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Cnew);
}

static VALUE rb_gsl_blas_zsyr2k(VALUE obj, VALUE u, VALUE t, VALUE a, VALUE aa,
        VALUE bb, VALUE b, VALUE cc)
{
  gsl_matrix_complex *A = NULL, *B = NULL, *C = NULL;
  gsl_complex alpha, beta, *pa = &alpha, *pb = &beta;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t Trans;
  CHECK_FIXNUM(u);  CHECK_FIXNUM(t);
  CHECK_COMPLEX(a);  CHECK_COMPLEX(b);
  CHECK_MATRIX_COMPLEX(aa);  CHECK_MATRIX_COMPLEX(bb);  CHECK_MATRIX_COMPLEX(cc);
  Uplo = FIX2INT(u);
  Trans = FIX2INT(t);
  Data_Get_Struct(a, gsl_complex, pa);
  Data_Get_Struct(b, gsl_complex, pb);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  Data_Get_Struct(bb, gsl_matrix_complex, B);
  Data_Get_Struct(cc, gsl_matrix_complex, C);
  gsl_blas_zsyr2k(Uplo, Trans, alpha, A, B, beta, C);
  return cc;
}

static VALUE rb_gsl_blas_zsyr2k2(VALUE obj, VALUE u, VALUE t, VALUE a, VALUE aa,
        VALUE bb, VALUE b, VALUE cc)
{
  gsl_matrix_complex *A = NULL, *B = NULL, *C = NULL, *Cnew = NULL;
  gsl_complex alpha, beta, *pa = &alpha, *pb = &beta;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t Trans;
  CHECK_FIXNUM(u);  CHECK_FIXNUM(t);
  CHECK_COMPLEX(a);  CHECK_COMPLEX(b);
  CHECK_MATRIX_COMPLEX(aa);  CHECK_MATRIX_COMPLEX(bb);  CHECK_MATRIX_COMPLEX(cc);
  Uplo = FIX2INT(u);
  Trans = FIX2INT(t);
  Data_Get_Struct(a, gsl_complex, pa);
  Data_Get_Struct(b, gsl_complex, pb);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  Data_Get_Struct(bb, gsl_matrix_complex, B);
  Data_Get_Struct(cc, gsl_matrix_complex, C);
  Cnew = gsl_matrix_complex_alloc(C->size1, C->size2);
  gsl_matrix_complex_memcpy(Cnew, C);
  gsl_blas_zsyr2k(Uplo, Trans, alpha, A, B, beta, Cnew);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, Cnew);
}

static VALUE rb_gsl_blas_zher2k(VALUE obj, VALUE u, VALUE t, VALUE a, VALUE aa,
        VALUE bb, VALUE b, VALUE cc)
{
  gsl_matrix_complex *A = NULL, *B = NULL, *C = NULL;
  gsl_complex alpha, *pa = &alpha;
  double beta;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t Trans;
  CHECK_FIXNUM(u);  CHECK_FIXNUM(t);
  CHECK_COMPLEX(a);
  Need_Float(b);
  CHECK_MATRIX_COMPLEX(aa);  CHECK_MATRIX_COMPLEX(bb);  CHECK_MATRIX_COMPLEX(cc);
  Uplo = FIX2INT(u);
  Trans = FIX2INT(t);
  Data_Get_Struct(a, gsl_complex, pa);
  beta = NUM2DBL(b);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  Data_Get_Struct(bb, gsl_matrix_complex, B);
  Data_Get_Struct(cc, gsl_matrix_complex, C);
  gsl_blas_zher2k(Uplo, Trans, alpha, A, B, beta, C);
  return cc;
}

static VALUE rb_gsl_blas_zher2k2(VALUE obj, VALUE u, VALUE t, VALUE a, VALUE aa,
        VALUE bb, VALUE b, VALUE cc)
{
  gsl_matrix_complex *A = NULL, *B = NULL, *C = NULL, *Cnew = NULL;
  gsl_complex alpha, *pa = &alpha;
  double beta;
  CBLAS_UPLO_t Uplo;
  CBLAS_TRANSPOSE_t Trans;
  CHECK_FIXNUM(u);  CHECK_FIXNUM(t);
  CHECK_COMPLEX(a);
  Need_Float(b);
  CHECK_MATRIX_COMPLEX(aa);  CHECK_MATRIX_COMPLEX(bb);  CHECK_MATRIX_COMPLEX(cc);
  Uplo = FIX2INT(u);
  Trans = FIX2INT(t);
  Data_Get_Struct(a, gsl_complex, pa);
  beta = NUM2DBL(b);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  Data_Get_Struct(bb, gsl_matrix_complex, B);
  Data_Get_Struct(cc, gsl_matrix_complex, C);
  Cnew = gsl_matrix_complex_alloc(C->size1, C->size2);
  gsl_matrix_complex_memcpy(Cnew, C);
  gsl_blas_zher2k(Uplo, Trans, alpha, A, B, beta, Cnew);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, Cnew);
}

void Init_gsl_blas3(VALUE module)
{
  rb_define_module_function(module, "dgemm", rb_gsl_blas_dgemm, -1);
  rb_define_module_function(module, "zgemm", rb_gsl_blas_zgemm, -1);

  rb_define_module_function(module, "dsymm", rb_gsl_blas_dsymm, -1);
  rb_define_module_function(module, "zsymm", rb_gsl_blas_zsymm, -1);
  rb_define_module_function(module, "zhemm", rb_gsl_blas_zhemm, -1);

  rb_define_module_function(module, "dtrmm!", rb_gsl_blas_dtrmm, 7);
  rb_define_module_function(module, "dtrmm", rb_gsl_blas_dtrmm2, 7);
  rb_define_module_function(module, "ztrmm!", rb_gsl_blas_ztrmm, 7);
  rb_define_module_function(module, "ztrmm", rb_gsl_blas_ztrmm2, 7);

  rb_define_module_function(module, "dtrsm!", rb_gsl_blas_dtrsm, 7);
  rb_define_module_function(module, "dtrsm", rb_gsl_blas_dtrsm2, 7);
  rb_define_module_function(module, "ztrsm!", rb_gsl_blas_ztrsm, 7);
  rb_define_module_function(module, "ztrsm", rb_gsl_blas_ztrsm2, 7);

  rb_define_module_function(module, "dsyrk!", rb_gsl_blas_dsyrk, 6);
  rb_define_module_function(module, "dsyrk", rb_gsl_blas_dsyrk2, 6);
  rb_define_module_function(module, "zsyrk!", rb_gsl_blas_zsyrk, 6);
  rb_define_module_function(module, "zsyrk", rb_gsl_blas_zsyrk2, 6);
  rb_define_module_function(module, "zherk!", rb_gsl_blas_zherk, 6);
  rb_define_module_function(module, "zherk", rb_gsl_blas_zherk2, 6);

  rb_define_module_function(module, "dsyr2k!", rb_gsl_blas_dsyr2k, 7);
  rb_define_module_function(module, "dsyr2k", rb_gsl_blas_dsyr2k2, 7);
  rb_define_module_function(module, "zsyr2k!", rb_gsl_blas_zsyr2k, 7);
  rb_define_module_function(module, "zsyr2k", rb_gsl_blas_zsyr2k2, 7);
  rb_define_module_function(module, "zher2k!", rb_gsl_blas_zher2k, 7);
  rb_define_module_function(module, "zher2k", rb_gsl_blas_zher2k2, 7);
}
