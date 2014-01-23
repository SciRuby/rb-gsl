/*
  linalg.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_config.h"
#include <gsl/gsl_math.h>
#include "rb_gsl_array.h"
#include "rb_gsl_common.h"
#include "rb_gsl_linalg.h"

static VALUE cgsl_matrix_LU;
static VALUE cgsl_matrix_QR;
static VALUE cgsl_matrix_QRPT;
static VALUE cgsl_vector_tau;
static VALUE cgsl_matrix_Q;
static VALUE cgsl_matrix_R;

static VALUE cgsl_matrix_LQ;
static VALUE cgsl_matrix_PTLQ;
static VALUE cgsl_matrix_L;

static VALUE cgsl_matrix_SV;
static VALUE cgsl_matrix_U;
static VALUE cgsl_matrix_V;
static VALUE cgsl_vector_S;
static VALUE cgsl_matrix_C;

enum {
  LINALG_DECOMP,
  LINALG_DECOMP_BANG,
};

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_LU_decomp_narray(int argc, VALUE *argv, VALUE obj,
					    int flag);
#endif

static VALUE rb_gsl_linalg_LU_decomposition(int argc, VALUE *argv, VALUE obj, int flag)
{
  gsl_matrix *mtmp = NULL, *m = NULL;
  gsl_permutation *p = NULL;
  int signum, itmp;
  size_t size;
  VALUE objp, objm, omatrix;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0])) 
      return rb_gsl_linalg_LU_decomp_narray(argc, argv, obj, flag);
#endif
    if (MATRIX_COMPLEX_P(argv[0]))
      return rb_gsl_linalg_complex_LU_decomp2(argc, argv, obj);
    omatrix = argv[0];
    itmp = 1;
    break;
  default:
    if (MATRIX_COMPLEX_P(obj))
      return rb_gsl_linalg_complex_LU_decomp2(argc, argv, obj);
    omatrix = obj;
    itmp = 0;
    break;
  }
  CHECK_MATRIX(omatrix);
  Data_Get_Struct(omatrix, gsl_matrix, mtmp);
  if (flag == LINALG_DECOMP_BANG) {
    m = mtmp;
    RBASIC(omatrix)->klass = cgsl_matrix_LU;
    objm = omatrix;
  } else {
    m = make_matrix_clone(mtmp);
    objm = Data_Wrap_Struct(cgsl_matrix_LU, 0, gsl_matrix_free, m);
  }
  size = m->size1;
  switch (argc-itmp) {
  case 0:
    p = gsl_permutation_alloc(size);
    gsl_linalg_LU_decomp(m, p, &signum);
    objp = Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p);
    if (flag == LINALG_DECOMP_BANG) return rb_ary_new3(2, objp, INT2FIX(signum));
    else return rb_ary_new3(3, objm, objp, INT2FIX(signum));
    break;
  case 1:
    CHECK_PERMUTATION(argv[itmp]);
    Data_Get_Struct(argv[itmp], gsl_permutation, p);
    gsl_linalg_LU_decomp(m, p, &signum);
    if (flag == LINALG_DECOMP_BANG) return INT2FIX(signum);
    else return rb_ary_new3(2, objm, INT2FIX(signum));
    break;
  default:
    rb_raise(rb_eArgError, "Usage: LU_decomp() or LU_decomp(permutation)");
    break;
  }
  return Qnil; /* never reach here */
}

static VALUE rb_gsl_linalg_LU_decomp(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_LU_decomposition(argc, argv, obj, LINALG_DECOMP);
}

static VALUE rb_gsl_linalg_LU_decomp_bang(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_LU_decomposition(argc, argv, obj, LINALG_DECOMP_BANG);
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_LU_decomp_narray(int argc, VALUE *argv, VALUE obj,
					    int flag)
{
  struct NARRAY *na, *na2;
  VALUE m;
  gsl_matrix_view mv;
  gsl_permutation *p;
  int signum;

  if (argc != 1) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
  GetNArray(argv[0], na);
  if (na->rank < 2) rb_raise(rb_eRuntimeError, "rank >= 2 required");
  if (na->shape[0] != na->shape[1])
    rb_raise(rb_eRuntimeError, "square matrix required");
  if (flag == LINALG_DECOMP) {
    m = na_make_object(NA_DFLOAT, 2, na->shape, CLASS_OF(argv[0]));
    GetNArray(m, na2);
    memcpy((double*)na2->ptr, (double*)na->ptr, sizeof(double)*na2->total);
    mv = gsl_matrix_view_array((double*)na2->ptr, na->shape[1], na->shape[0]);
  } else {
    mv = gsl_matrix_view_array((double*)na->ptr, na->shape[1], na->shape[0]);
  }
  p = gsl_permutation_alloc(mv.matrix.size1);
  gsl_linalg_LU_decomp(&mv.matrix, p, &signum);
  if (flag == LINALG_DECOMP) {
    return rb_ary_new3(3, m, 
		       Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p),
		       INT2FIX(signum));
  } else {
    return rb_ary_new3(3, argv[0], 
		Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p),
		INT2FIX(signum));
  }
  
}
#endif

static gsl_matrix* get_matrix(VALUE obj, VALUE klass,int *flagm);
static gsl_permutation* get_permutation(VALUE obj,  size_t size, int *flagp);
static gsl_vector* get_vector2(VALUE obj,  int *flagv);

static gsl_matrix* get_matrix(VALUE obj, VALUE klass, int *flagm)
{
  gsl_matrix *mtmp = NULL, *m = NULL;
#ifdef HAVE_NARRAY_H
  gsl_matrix_view mv;
  struct NARRAY *na;
#endif
  if (CLASS_OF(obj) == klass) {
    Data_Get_Struct(obj, gsl_matrix, m);
    *flagm = 0;
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(obj)) {
    GetNArray(obj, na);
    mv = gsl_matrix_view_array((double*)na->ptr, na->shape[1], na->shape[0]);
    m = &mv.matrix;
    *flagm = -1;
#endif
  } else {
    CHECK_MATRIX(obj);
    Data_Get_Struct(obj, gsl_matrix, mtmp);
    m = make_matrix_clone(mtmp);
    *flagm = 1;
  }
  return m;
}

static gsl_permutation* get_permutation(VALUE obj, size_t size, int *flagp)
{
  gsl_permutation *p = NULL;
  if (CLASS_OF(obj) == cgsl_permutation) {
    Data_Get_Struct(obj, gsl_permutation, p);
    *flagp = 0;
  } else {
    p = gsl_permutation_alloc(size);
    *flagp = 1;
  }
  return p;
}

static gsl_vector* get_vector2(VALUE obj, int *flagv)
{
  gsl_vector *v = NULL;
#ifdef HAVE_NARRAY_H
  gsl_vector_view vv;
  struct NARRAY *na;
#endif
  if (TYPE(obj) == T_ARRAY) {
    v = make_cvector_from_rarray(obj);
    *flagv = 1;
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(obj)) {
    GetNArray(obj, na);
    vv = gsl_vector_view_array((double*) na->ptr, na->total);
    v = &vv.vector;
    *flagv = -1;
#endif
  } else {
    CHECK_VECTOR(obj);
    Data_Get_Struct(obj, gsl_vector, v);
    *flagv = 0;
  }
  return v;
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_LU_solve_narray(int argc, VALUE *argv, VALUE obj);
#endif

VALUE rb_gsl_linalg_LU_solve(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m = NULL;
  gsl_permutation *p = NULL;
  gsl_vector *b = NULL, *x = NULL;
  int signum, flagm = 0, flagp = 0, flagb = 0, flagx = 0, itmp;
  size_t size;
  VALUE bb;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 2 || argc > 4) 
      rb_raise(rb_eArgError, "Usage: solve(m, b), solve(m, b, x), solve(lu, p, b), solve(lu, p, b, x)");
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0])) 
      return rb_gsl_linalg_LU_solve_narray(argc, argv, obj);
#endif
    m = get_matrix(argv[0], cgsl_matrix_LU, &flagm);
    itmp = 1;
    break;
  default:
    if (argc < 1 || argc > 3) 
      rb_raise(rb_eArgError, "Usage: LU_solve(b), LU_solve(p, b), LU_solve(b, x), solve(p, b, x)");
    
    m = get_matrix(obj, cgsl_matrix_LU, &flagm);
    itmp = 0;
    break;
  }
  size = m->size1;

  p = get_permutation(argv[itmp], size, &flagp);
  if (flagp == 1 && flagm == 0) rb_raise(rb_eArgError, "permutation must be given");
  if (flagp == 0) itmp++;

  bb = argv[itmp];
  b = get_vector2(argv[itmp], &flagb);
  itmp++;

  if (itmp == argc) {
    x = gsl_vector_alloc(size);
    flagx = 1;
  } else {
    CHECK_VECTOR(argv[itmp]);
    Data_Get_Struct(argv[itmp], gsl_vector, x);
  }
  if (flagm == 1) gsl_linalg_LU_decomp(m, p, &signum);
  gsl_linalg_LU_solve(m, p, b, x);
  if (flagm == 1) gsl_matrix_free(m);
  if (flagp == 1) gsl_permutation_free(p);
  if (flagb == 1) gsl_vector_free(b);
  if (flagx == 1) return Data_Wrap_Struct(VECTOR_ROW_COL(bb), 0, gsl_vector_free, x);
  else return argv[argc-1];
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_LU_solve_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na, *b;
  VALUE ret;
  gsl_permutation *p;
  gsl_matrix_view mv;
  gsl_vector_view bv, xv;
  double *x;
  int shape[1];
  if (argc < 3)
    rb_raise(rb_eArgError, 
	     "wrong number of arguments %d(NArray, GSL::Permutation and NArray expected", 
	     argc);
  GetNArray(argv[0], na);
  mv = gsl_matrix_view_array((double*) na->ptr, na->shape[1], na->shape[0]);
  CHECK_PERMUTATION(argv[1]);
  Data_Get_Struct(argv[1], gsl_permutation, p);
  GetNArray(argv[2], b);
  bv = gsl_vector_view_array((double*) b->ptr, b->total);
  if (argc == 3) {
    shape[0] = b->total;
    ret = na_make_object(NA_DFLOAT, 1, shape, CLASS_OF(argv[0]));
  } else {
    ret = argv[3];
  }
  x = NA_PTR_TYPE(ret,double*);
  xv = gsl_vector_view_array(x, b->total);
  gsl_linalg_LU_solve(&mv.matrix, p, &bv.vector, &xv.vector);
  return ret;
}
#endif

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_LU_svx_narray(int argc, VALUE *argv, VALUE obj);
#endif

/* bb must be Vector, it is replaced by the root of the system */
static VALUE rb_gsl_linalg_LU_svx(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m = NULL;
  gsl_permutation *p = NULL;
  gsl_vector *b = NULL;
  int signum, flagm = 0, flagp = 0, flagb = 0, itmp;
  size_t size;

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 2 || argc > 3) 
      rb_raise(rb_eArgError, "Usage: svx(m, b), svx(lu, p, b)");
#ifdef HAVE_NARRAY_H
  if (NA_IsNArray(argv[0])) 
    return rb_gsl_linalg_LU_svx_narray(argc, argv, obj);
#endif
    m = get_matrix(argv[0], cgsl_matrix_LU, &flagm);
    itmp = 1;
    break;
  default:
   if (argc < 1 || argc > 2) 
     rb_raise(rb_eArgError, "Usage: LU_svx(b), LU_svx(p, b)");
    m = get_matrix(obj, cgsl_matrix_LU, &flagm);
    itmp = 0;
    break;
  }
  size = m->size1;
  p = get_permutation(argv[itmp], size, &flagp);
  if (flagp == 1 && flagm == 0) rb_raise(rb_eArgError, "permutation must be given");
  if (flagp == 0) itmp++;
  CHECK_VECTOR(argv[itmp]);
  b = get_vector2(argv[itmp], &flagb);
  if (flagm == 1) gsl_linalg_LU_decomp(m, p, &signum);
  gsl_linalg_LU_svx(m, p, b);
  if (flagm == 1) gsl_matrix_free(m);
  if (flagp == 1) gsl_permutation_free(p);
  return argv[itmp];
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_LU_svx_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na, *b;
  gsl_permutation *p;
  gsl_matrix_view mv;
  gsl_vector_view bv;
  if (argc != 3)
    rb_raise(rb_eArgError, 
	     "wrong number of arguments %d(NArray, GSL::Permutation and NArray expected", 
	     argc);
  GetNArray(argv[0], na);
  mv = gsl_matrix_view_array((double*) na->ptr, na->shape[1], na->shape[0]);
  CHECK_PERMUTATION(argv[1]);
  Data_Get_Struct(argv[1], gsl_permutation, p);
  GetNArray(argv[2], b);
  bv = gsl_vector_view_array((double*) b->ptr, b->total);
  gsl_linalg_LU_svx(&mv.matrix, p, &bv.vector);
  return argv[2];
}
#endif

/* singleton */
static VALUE rb_gsl_linalg_LU_refine(VALUE obj, VALUE vm,
				     VALUE lu, VALUE pp, VALUE bb,
				     VALUE xx)
{
  gsl_matrix *m = NULL, *mlu = NULL;
  gsl_permutation *p = NULL;
  gsl_vector *b = NULL, *x = NULL, *r = NULL;
  int flagb = 0;
  VALUE vr;
  CHECK_MATRIX(vm);  CHECK_MATRIX(lu);
  CHECK_PERMUTATION(pp);  CHECK_VECTOR(xx);
  Data_Get_Struct(vm, gsl_matrix, m);
  Data_Get_Struct(lu, gsl_matrix, mlu);
  Data_Get_Struct(pp, gsl_permutation, p);
  if (TYPE(bb) == T_ARRAY) {
    b = make_cvector_from_rarray(bb);
    flagb = 1;
  } else {
    CHECK_VECTOR(bb);
    Data_Get_Struct(bb, gsl_vector, b);
  }
  Data_Get_Struct(xx, gsl_vector, x);
  r = gsl_vector_alloc(m->size1);
  gsl_linalg_LU_refine(m, mlu, p, b, x, r);
  vr = Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, r);
  if (flagb == 1) gsl_vector_free(b);
  return rb_ary_new3(2, xx, vr);
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_LU_invert_narray(int argc, VALUE *argv, VALUE obj);
#endif

static VALUE rb_gsl_linalg_LU_invert(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m = NULL, *inverse = NULL;
  gsl_permutation *p = NULL;
  int signum, flagm = 0, flagp = 0, itmp;
  size_t size;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0]))
      return rb_gsl_linalg_LU_invert_narray(argc, argv, obj);
#endif
    m = get_matrix(argv[0], cgsl_matrix_LU, &flagm);
    itmp = 1;
    break;
  default:
    m = get_matrix(obj, cgsl_matrix_LU, &flagm);
    itmp = 0;
  }
  size = m->size1;

  if (argc == itmp) {
    p = gsl_permutation_alloc(size);
    flagp = 1;
  } else {
    CHECK_PERMUTATION(argv[itmp]);
    p = get_permutation(argv[itmp], size, &flagp);
  }
  if (flagp == 1 && flagm == 0) rb_raise(rb_eArgError, "permutation must be given");
  if (flagp == 0) itmp++;

  if (flagm == 1 || flagp == 1) {
    gsl_linalg_LU_decomp(m, p, &signum);  
  }

  if (argc-1 == itmp) {
    CHECK_MATRIX(argv[itmp]);
    Data_Get_Struct(argv[itmp], gsl_matrix, inverse);
  } else {
    inverse = gsl_matrix_alloc(size, size);
  }
  gsl_linalg_LU_invert(m, p, inverse);
  if (flagm == 1) gsl_matrix_free(m);
  if (flagp == 1) gsl_permutation_free(p);
  if (argc-1 == itmp) return argv[itmp];
  else return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, inverse);

}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_LU_invert_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na;
  VALUE inv;
  gsl_permutation *p;
  gsl_matrix_view mv1, mv2;
  if (argc != 2) {
    rb_raise(rb_eArgError, "Usage: LU.invert(lu, perm)");
  }
  CHECK_PERMUTATION(argv[1]);
  GetNArray(argv[0], na);
  inv = na_make_object(NA_DFLOAT, 2, na->shape, CLASS_OF(argv[0]));
  mv1 = gsl_matrix_view_array((double*)na->ptr, na->shape[1], na->shape[0]);
  mv2 = gsl_matrix_view_array(NA_PTR_TYPE(inv, double*), na->shape[1], na->shape[0]);
  CHECK_PERMUTATION(argv[1]);
  Data_Get_Struct(argv[1], gsl_permutation, p);
  gsl_linalg_LU_invert(&mv1.matrix, p, &mv2.matrix);
  return inv;
}
static VALUE rb_gsl_linalg_LU_det_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na;
  gsl_matrix_view mv;
  int signum = 1;
  switch (argc) {
  case 2:
    signum = FIX2INT(argv[1]);
    /* no break */
  case 1:
    GetNArray(argv[0], na);
    mv = gsl_matrix_view_array((double*)na->ptr, na->shape[1], na->shape[0]);
    break;
  default:
    rb_raise(rb_eArgError, "Usage: LU.det(lu, perm)");
    break;
  }
  return rb_float_new(gsl_linalg_LU_det(&mv.matrix, signum));
}
static VALUE rb_gsl_linalg_LU_lndet_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na;
  gsl_matrix_view mv;
  switch (argc) {
  case 1:
    GetNArray(argv[0], na);
    mv = gsl_matrix_view_array((double*)na->ptr, na->shape[1], na->shape[0]);
    break;
  default:
    rb_raise(rb_eArgError, "Usage: LU.lndet(lu)");
    break;
  }
  return rb_float_new(gsl_linalg_LU_lndet(&mv.matrix));
}

#endif

static VALUE rb_gsl_linalg_LU_det(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m = NULL;
  gsl_permutation *p = NULL;
  int flagm = 0, flagp = 0, itmp, sign;
  size_t size;
  double det;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0]))
      return rb_gsl_linalg_LU_det_narray(argc, argv, obj);
#endif

    m = get_matrix(argv[0], cgsl_matrix_LU, &flagm);
    itmp = 1;
    break;
  default:
    m = get_matrix(obj, cgsl_matrix_LU, &flagm);
    itmp = 0;
    break;
  }
  size = m->size1;
  if (flagm == 0) {
    if (argc-itmp == 1)  sign = FIX2INT(argv[itmp]);
    else sign = 1;
  } else {
    if (argc-itmp >= 2) {
      get_permutation(argv[itmp], size, &flagp);
    } else {
      p = gsl_permutation_alloc(size);
      flagp = 1;
    }
  } 
  if (flagm == 1) gsl_linalg_LU_decomp(m, p, &sign);  
  det = gsl_linalg_LU_det(m, sign);
  if (flagm == 1) gsl_matrix_free(m);
  if (flagp == 1) gsl_permutation_free(p);
  return rb_float_new(det);
}

static VALUE rb_gsl_linalg_LU_lndet(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m = NULL;
  gsl_permutation *p = NULL;
  int flagm = 0, sign;
  double lndet;

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0]))
      return rb_gsl_linalg_LU_lndet_narray(argc, argv, obj);
#endif

    m = get_matrix(argv[0], cgsl_matrix_LU, &flagm);
    break;
  default:
    m = get_matrix(obj, cgsl_matrix_LU, &flagm);
    break;
  }
  if (flagm == 1) {
    p = gsl_permutation_alloc(m->size1);
    gsl_linalg_LU_decomp(m, p, &sign);  
  }
  lndet = gsl_linalg_LU_lndet(m);
  if (flagm == 1) {
    gsl_matrix_free(m);
    gsl_permutation_free(p);
  }
  return rb_float_new(lndet);
}

static VALUE rb_gsl_linalg_LU_sgndet(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m = NULL;
  gsl_permutation *p = NULL;
  int flagm = 0, sign, signdet, itmp;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    m = get_matrix(argv[0], cgsl_matrix_LU, &flagm);
    itmp = 1;
    break;
  default:
    m = get_matrix(obj, cgsl_matrix_LU, &flagm);
    itmp = 0;
    break;
  }
  if (flagm == 1) {
    p = gsl_permutation_alloc(m->size1);
    gsl_linalg_LU_decomp(m, p, &sign);  
  } else {
    if (argc-itmp != 1) rb_raise(rb_eArgError, "sign must be given");
    sign = FIX2INT(argv[itmp]);
  }
  signdet = gsl_linalg_LU_sgndet(m, sign);
  if (flagm == 1) {
    gsl_matrix_free(m);
    gsl_permutation_free(p);
  }
  return INT2FIX(signdet);
}

#ifdef GSL_1_6_LATER
int gsl_linalg_LQ_solve_T(const gsl_matrix*, const gsl_vector*, const gsl_vector*, gsl_vector*);
int gsl_linalg_LQ_svx_T (const gsl_matrix*, const gsl_vector*, gsl_vector*);
int gsl_linalg_LQ_lssolve_T(const gsl_matrix * LQ, const gsl_vector * tau, 
                           const gsl_vector * b, gsl_vector * x, 
                           gsl_vector * residual);
int
gsl_linalg_LQ_Lsolve_T (const gsl_matrix * LQ, const gsl_vector * b, gsl_vector* x);
int
gsl_linalg_LQ_Lsvx_T (const gsl_matrix * LQ, gsl_vector * x);
int
gsl_linalg_L_solve_T (const gsl_matrix * L, const gsl_vector * b, gsl_vector * x);


#endif

enum {
  LINALG_QR_DECOMP,
  LINALG_QR_DECOMP_BANG,
  LINALG_LQ_DECOMP,
  LINALG_LQ_DECOMP_BANG,
  LINALG_QR_SOLVE,
  LINALG_LQ_SOLVE,
  LINALG_QR_QTvec,
  LINALG_QR_Qvec,
  LINALG_LQ_vecQ,
  LINALG_LQ_vecQT,
  LINALG_QR_RSOLVE,
  LINALG_LQ_LSOLVE,
  LINALG_QR_RSVX,
  LINALG_LQ_LSVX,
  LINALG_R_SOLVE,
  LINALG_R_SVX,
  LINALG_L_SOLVE,
  LINALG_L_SVX,
  LINALG_QR_UNPACK,
  LINALG_LQ_UNPACK,
};

static VALUE rb_gsl_linalg_QR_LQ_decomposition(int argc, VALUE *argv, VALUE obj,
					       int flag)
{
  gsl_matrix *m = NULL, *mtmp = NULL;
  gsl_vector *tau = NULL;
  int (*fdecomp)(gsl_matrix *, gsl_vector *);
  int itmp, status;
  VALUE vtau, mdecomp, omatrix;

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments.");
    omatrix = argv[0];
    itmp = 1;
    break;
  default:
    omatrix = obj;
    itmp = 0;
    break;
  }
  CHECK_MATRIX(omatrix);
  Data_Get_Struct(omatrix, gsl_matrix, mtmp);    

  switch (flag) {
  case LINALG_QR_DECOMP:
    fdecomp = &gsl_linalg_QR_decomp;
    m = make_matrix_clone(mtmp);
    mdecomp = Data_Wrap_Struct(cgsl_matrix_QR, 0, gsl_matrix_free, m);
    break;
  case LINALG_QR_DECOMP_BANG:
    fdecomp = &gsl_linalg_QR_decomp;
    m = mtmp;
    mdecomp = omatrix;
    RBASIC(mdecomp)->klass = cgsl_matrix_QR;
    break;
#ifdef GSL_1_6_LATER
  case LINALG_LQ_DECOMP:
    fdecomp = &gsl_linalg_LQ_decomp;
    m = make_matrix_clone(mtmp);
    mdecomp = Data_Wrap_Struct(cgsl_matrix_LQ, 0, gsl_matrix_free, m);
    break;
  case LINALG_LQ_DECOMP_BANG:
    fdecomp = &gsl_linalg_LQ_decomp;
    m = mtmp;
    mdecomp = omatrix;
    RBASIC(mdecomp)->klass = cgsl_matrix_LQ;
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  switch (argc - itmp) {
  case 0:
    tau = gsl_vector_alloc(GSL_MIN(mtmp->size1, mtmp->size2));
    break;
  case 1:
    CHECK_VECTOR(argv[itmp]);
    Data_Get_Struct(argv[itmp], gsl_vector, tau);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  status = (*fdecomp)(m, tau);
  switch (flag) {
  case LINALG_QR_DECOMP:
  case LINALG_LQ_DECOMP:
    if (argc == itmp) {
      vtau = Data_Wrap_Struct(cgsl_vector_tau, 0, gsl_vector_free, tau);
      return rb_ary_new3(2, mdecomp, vtau);
    } else {
      RBASIC(argv[itmp])->klass = cgsl_vector_tau;
      return mdecomp;
    }
    break;
  case LINALG_QR_DECOMP_BANG:
  case LINALG_LQ_DECOMP_BANG:
   if (argc == itmp) {
      return Data_Wrap_Struct(cgsl_vector_tau, 0, gsl_vector_free, tau);
    } else {
      RBASIC(argv[itmp])->klass = cgsl_vector_tau;
      return INT2FIX(status);
    }
    break;
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  return Qnil;
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_QR_decomp_narray(int argc, VALUE *argv, VALUE obj);
#endif

static VALUE rb_gsl_linalg_QR_decomp(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_NARRAY_H
  if (argc >= 1 && NA_IsNArray(argv[0]))
    return rb_gsl_linalg_QR_decomp_narray(argc, argv, obj);
#endif
  return rb_gsl_linalg_QR_LQ_decomposition(argc, argv, obj, LINALG_QR_DECOMP);
}

static VALUE rb_gsl_linalg_QR_decomp_bang(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QR_LQ_decomposition(argc, argv, obj, LINALG_QR_DECOMP_BANG);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_LQ_decomp(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QR_LQ_decomposition(argc, argv, obj, LINALG_LQ_DECOMP);
}

static VALUE rb_gsl_linalg_LQ_decomp_bang(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QR_LQ_decomposition(argc, argv, obj, LINALG_LQ_DECOMP_BANG);
}
#endif

static VALUE rb_gsl_linalg_QR_LQ_solve(int argc, VALUE *argv, VALUE obj, int flag)
{
  gsl_matrix *m = NULL;
  gsl_vector *b = NULL, *x = NULL, *tau = NULL;
  VALUE omatrix;
  int flagm = 0, flagt = 0, flagb = 0, flagx = 0, itmp;
  size_t size;
  int (*fdecomp)(gsl_matrix*, gsl_vector*);
  int (*fsolve)(const gsl_matrix*, const gsl_vector*, const gsl_vector*, gsl_vector*);

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments.");
    omatrix = argv[0];
    itmp = 1;
    break;
  default:
    omatrix = obj;
    itmp = 0;
    break;
  }
  if (argc-itmp < 1 || argc-itmp > 3)
    rb_raise(rb_eArgError,  "wrong number of arguments");
  CHECK_MATRIX(omatrix);
  switch (flag) {
  case LINALG_QR_SOLVE:
    m = get_matrix(omatrix, cgsl_matrix_QR, &flagm);
    fdecomp = &gsl_linalg_QR_decomp;
    fsolve = &gsl_linalg_QR_solve;
    break;
#ifdef GSL_1_6_LATER
  case LINALG_LQ_SOLVE:
    m = get_matrix(omatrix, cgsl_matrix_LQ, &flagm);
    fdecomp = &gsl_linalg_LQ_decomp;
    fsolve = &gsl_linalg_LQ_solve_T;
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operatioin");
    break;
  }
  size = m->size1;
  if (flagm == 0) { /* the matrix given is already decomped */
    if (CLASS_OF(argv[itmp]) != cgsl_vector_tau)
      rb_raise(rb_eArgError, "tau vector must be given");
    Data_Get_Struct(argv[itmp], gsl_vector, tau);
    flagt = 0;
    itmp++;
  } else {
    if (CLASS_OF(argv[itmp]) == cgsl_vector_tau) {
      Data_Get_Struct(argv[itmp], gsl_vector, tau);
      flagt = 0;
      itmp++;
    } else {
      tau = gsl_vector_alloc(size);
      flagt = 1;
    }
  }
  b = get_vector2(argv[itmp], &flagb);
  itmp++;
  if (itmp == argc) {
    x = gsl_vector_alloc(m->size1);
    flagx = 1;
  } else {
    CHECK_VECTOR(argv[itmp]);
    Data_Get_Struct(argv[itmp], gsl_vector, x);
    flagx = 0;
  }
  if (flagm == 1) (*fdecomp)(m, tau);
  (*fsolve)(m, tau, b, x);
  if (flagm == 1) gsl_matrix_free(m);
  if (flagt == 1) gsl_vector_free(tau);
  if (flagb == 1) gsl_vector_free(b);
  if (flagx == 1) return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
  else return argv[itmp];
}

static VALUE rb_gsl_linalg_QR_LQ_svx(int argc, VALUE *argv, VALUE obj, int flag)
{
  gsl_matrix *m = NULL;
  gsl_vector *b = NULL, *tau = NULL;
  VALUE omatrix;
  int flagm = 0, flagt = 0, flagb = 0, itmp;
  size_t size;
  int (*fdecomp)(gsl_matrix*, gsl_vector*);
  int (*fsvx)(const gsl_matrix*, const gsl_vector*, gsl_vector*);

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments.");
    omatrix = argv[0];
    itmp = 1;
    break;
  default:
    omatrix = obj;
    itmp = 0;
    break;
  }
  if (argc-itmp < 1 || argc-itmp > 2) 
    rb_raise(rb_eArgError,  "wrong number of arguments");
  CHECK_MATRIX(omatrix);
  switch (flag) {
  case LINALG_QR_SOLVE:
    m = get_matrix(omatrix, cgsl_matrix_QR, &flagm);
    fdecomp = &gsl_linalg_QR_decomp;
    fsvx = &gsl_linalg_QR_svx;
    break;
#ifdef GSL_1_6_LATER
  case LINALG_LQ_SOLVE:
    m = get_matrix(omatrix, cgsl_matrix_LQ, &flagm);
    fdecomp = &gsl_linalg_LQ_decomp;
    fsvx = &gsl_linalg_LQ_svx_T;
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operatioin");
    break;
  }
  size = m->size1;
  if (flagm == 0) { /* the matrix given is already decomped */
    if (CLASS_OF(argv[itmp]) != cgsl_vector_tau)
      rb_raise(rb_eArgError, "tau vector must be given");
    Data_Get_Struct(argv[itmp], gsl_vector, tau);
    flagt = 0;
    itmp++;
  } else {
    if (CLASS_OF(argv[itmp]) == cgsl_vector_tau) {
      Data_Get_Struct(argv[itmp], gsl_vector, tau);
      flagt = 0;
      itmp++;
    } else {
      tau = gsl_vector_alloc(size);
      flagt = 1;
    }
  }
  b = get_vector2(argv[itmp], &flagb);
  if (flagm == 1 && flagt == 1) (*fdecomp)(m, tau);
  (*fsvx)(m, tau, b);
  if (flagm == 1) gsl_matrix_free(m);
  if (flagt == 1) gsl_vector_free(tau);
  return argv[itmp];
}

static VALUE rb_gsl_linalg_QR_LQ_lssolve(int argc, VALUE *argv, VALUE obj, int flag)
{
  gsl_matrix *m = NULL;
  gsl_vector *b = NULL, *x = NULL, *tau = NULL, *r = NULL;
  VALUE omatrix;
  int flagm = 0, flagt = 0, flagb = 0, itmp, status;
  size_t size;
  int (*fdecomp)(gsl_matrix*, gsl_vector*);
  int (*flssolve)(const gsl_matrix*, const gsl_vector*, const gsl_vector*, gsl_vector*,
		  gsl_vector*);

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments.");
    omatrix = argv[0];
    itmp = 1;
    break;
  default:
    omatrix = obj;
    itmp = 0;
    break;
  }
  if (argc-itmp < 1 || argc-itmp > 4)
    rb_raise(rb_eArgError,  "wrong number of arguments");
  CHECK_MATRIX(omatrix);
  switch (flag) {
  case LINALG_QR_SOLVE:
    m = get_matrix(omatrix, cgsl_matrix_QR, &flagm);
    fdecomp = &gsl_linalg_QR_decomp;
    flssolve = &gsl_linalg_QR_lssolve;
    break;
#ifdef GSL_1_6_LATER
  case LINALG_LQ_SOLVE:
    m = get_matrix(omatrix, cgsl_matrix_LQ, &flagm);
    fdecomp = &gsl_linalg_LQ_decomp;
    flssolve = &gsl_linalg_LQ_lssolve_T;
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operatioin");
    break;
  }
  size = m->size1;
  if (flagm == 0) { /* the matrix given is already decomped */
    if (CLASS_OF(argv[itmp]) != cgsl_vector_tau)
      rb_raise(rb_eArgError, "tau vector must be given");
    Data_Get_Struct(argv[itmp], gsl_vector, tau);
    flagt = 0;
    itmp++;
  } else {
    if (CLASS_OF(argv[itmp]) == cgsl_vector_tau) {
      Data_Get_Struct(argv[itmp], gsl_vector, tau);
      flagt = 0;
      itmp++;
    } else {
      tau = gsl_vector_alloc(size);
      flagt = 1;
    }
  }
  b = get_vector2(argv[itmp], &flagb);
  itmp++;
  switch (argc - itmp) {
  case 2:
    CHECK_VECTOR(argv[argc-2]);
    Data_Get_Struct(argv[argc-2], gsl_vector, x);
    CHECK_VECTOR(argv[argc-1]);
    Data_Get_Struct(argv[argc-1], gsl_vector, r);
    break;
  case 1:
    CHECK_VECTOR(argv[argc-1]);
    Data_Get_Struct(argv[argc-1], gsl_vector, x);
    r = gsl_vector_alloc(x->size);
    break;
  case 0:
    x = gsl_vector_alloc(m->size1);
    r = gsl_vector_alloc(m->size1);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  if (flagm == 1) (*fdecomp)(m, tau);
  status = (*flssolve)(m, tau, b, x, r);
  if (flagm == 1) gsl_matrix_free(m);
  if (flagt == 1) gsl_vector_free(tau);
  if (flagb == 1) gsl_vector_free(b);

  switch (argc - itmp) {
  case 2:
    return INT2FIX(status);
    break;
  case 1:
    return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, r);
    break;
  default:
    return rb_ary_new3(2, Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x),
		       Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, r));
  }
  return Qnil;
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_QR_solve_narray(int argc, VALUE *argv, VALUE obj);
static VALUE rb_gsl_linalg_QR_svx_narray(int argc, VALUE *argv, VALUE obj);
#endif

static VALUE rb_gsl_linalg_QR_solve(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_NARRAY_H
  if (argc == 3 && NA_IsNArray(argv[0]))
    return rb_gsl_linalg_QR_solve_narray(argc, argv, obj);
#endif
  return rb_gsl_linalg_QR_LQ_solve(argc, argv, obj, LINALG_QR_SOLVE);
}

static VALUE rb_gsl_linalg_QR_svx(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_NARRAY_H
  if (argc == 2 && NA_IsNArray(argv[0]))
    return rb_gsl_linalg_QR_svx_narray(argc, argv, obj);
#endif
  return rb_gsl_linalg_QR_LQ_svx(argc, argv, obj, LINALG_QR_SOLVE);
}

static VALUE rb_gsl_linalg_QR_lssolve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QR_LQ_lssolve(argc, argv, obj, LINALG_QR_SOLVE);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_LQ_solve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QR_LQ_solve(argc, argv, obj, LINALG_LQ_SOLVE);
}

static VALUE rb_gsl_linalg_LQ_svx(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QR_LQ_svx(argc, argv, obj, LINALG_LQ_SOLVE);
}

static VALUE rb_gsl_linalg_LQ_lssolve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QR_LQ_lssolve(argc, argv, obj, LINALG_LQ_SOLVE);
}
#endif

static VALUE rb_gsl_linalg_QRLQ_QTvec(int argc, VALUE *argv, VALUE obj,
				      int flag)
{
  gsl_matrix *QR = NULL;
  gsl_vector *tau = NULL, *v = NULL;
  VALUE ret;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 3) rb_raise(rb_eArgError, 
			    "wrong number of arguments (%d for 3)", argc); 
    CHECK_MATRIX(argv[0]); CHECK_VECTOR(argv[1]); CHECK_VECTOR(argv[2]);
    Data_Get_Struct(argv[0], gsl_matrix, QR);
    Data_Get_Struct(argv[1], gsl_vector, tau);
    Data_Get_Struct(argv[2], gsl_vector, v);
    ret = argv[2];
    break;
  default:
    if (argc != 2) rb_raise(rb_eArgError, 
			    "wrong number of arguments (%d for 2)", argc); 
    CHECK_VECTOR(argv[2]); CHECK_VECTOR(argv[1]);
    Data_Get_Struct(obj, gsl_matrix, QR);
    Data_Get_Struct(argv[0], gsl_vector, tau);
    Data_Get_Struct(argv[1], gsl_vector, v);
    ret = argv[1];
    break;
  }
  switch (flag) {
  case LINALG_QR_QTvec:
    gsl_linalg_QR_QTvec(QR, tau, v);
    break;
  case LINALG_QR_Qvec:
    gsl_linalg_QR_Qvec(QR, tau, v);
    break;
#ifdef GSL_1_6_LATER
  case LINALG_LQ_vecQ:
    gsl_linalg_LQ_vecQ(QR, tau, v);
    break;
  case LINALG_LQ_vecQT:
    gsl_linalg_LQ_vecQT(QR, tau, v);
    break;
#endif
  default:
    break;
  }
  return ret;
}

static VALUE rb_gsl_linalg_QR_QTvec(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQ_QTvec(argc, argv, obj, LINALG_QR_QTvec);
}

static VALUE rb_gsl_linalg_QR_Qvec(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQ_QTvec(argc, argv, obj, LINALG_QR_Qvec);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_LQ_vecQT(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQ_QTvec(argc, argv, obj, LINALG_LQ_vecQT);
}

static VALUE rb_gsl_linalg_LQ_vecQ(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQ_QTvec(argc, argv, obj, LINALG_LQ_vecQ);
}
#endif

static VALUE rb_gsl_linalg_QRLQ_unpack(int argc, VALUE *argv, VALUE obj,
				       int flag)
{
  gsl_matrix *QR = NULL, *Q = NULL, *R = NULL;
  gsl_vector *tau = NULL;
  int itmp;
  VALUE vtmp, vQ, vR, klass;
  switch (flag) {
  case LINALG_QR_UNPACK:
    klass = cgsl_matrix_QR;
    break;
  case LINALG_LQ_UNPACK:
    klass = cgsl_matrix_LQ;
    break;
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, 
			    "wrong number of arguments (%d for 2)", argc);
    vtmp = argv[0];
    itmp = 1;
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, 
			    "wrong number of arguments (%d for 1)", argc);
    vtmp = obj;
    itmp = 0;
    break;
  }
  CHECK_MATRIX(vtmp);
  if (CLASS_OF(vtmp) != klass) {
    rb_raise(rb_eTypeError, "not a QR matrix");
  }
  Data_Get_Struct(vtmp, gsl_matrix, QR);
  if (CLASS_OF(argv[itmp]) != cgsl_vector_tau)
    rb_raise(rb_eTypeError, "tau vector must be given.");
  Data_Get_Struct(argv[itmp], gsl_vector, tau);
  Q = gsl_matrix_alloc(QR->size1, QR->size1);
  R = gsl_matrix_alloc(QR->size1, QR->size2);
  switch (flag) {
  case LINALG_QR_UNPACK:
    gsl_linalg_QR_unpack(QR, tau, Q, R);
    vQ = Data_Wrap_Struct(cgsl_matrix_Q, 0, gsl_matrix_free, Q);
    vR = Data_Wrap_Struct(cgsl_matrix_R, 0, gsl_matrix_free, R);
    break;
#ifdef GSL_1_6_LATER
  case LINALG_LQ_UNPACK:
  gsl_linalg_LQ_unpack(QR, tau, Q, R);
  vQ = Data_Wrap_Struct(cgsl_matrix_L, 0, gsl_matrix_free, Q);
  vR = Data_Wrap_Struct(cgsl_matrix_Q, 0, gsl_matrix_free, R);
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  return rb_ary_new3(2, vQ, vR);
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_QR_unpack_narray(int argc, VALUE *argv, VALUE obj);
#endif

static VALUE rb_gsl_linalg_QR_unpack(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_NARRAY_H
  if (argc == 2 && NA_IsNArray(argv[0]))
    return rb_gsl_linalg_QR_unpack_narray(argc, argv, obj);
#endif
  return rb_gsl_linalg_QRLQ_unpack(argc, argv, obj, LINALG_QR_UNPACK);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_LQ_unpack(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQ_unpack(argc, argv, obj, LINALG_LQ_UNPACK);
}
#endif

/* singleton */
static VALUE rb_gsl_linalg_QRLQ_QRLQsolve(int argc, VALUE *argv, VALUE obj,
					  int flag)
{
  gsl_matrix *Q = NULL, *R = NULL;
  gsl_vector *b = NULL, *x = NULL;
  int (*fsolve)(gsl_matrix*, gsl_matrix *, const gsl_vector*, gsl_vector *);
  int flagb = 0;
  VALUE retval;
  switch (argc) {
  case 3:
    CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
    Data_Get_Struct(argv[0], gsl_matrix, Q);
    Data_Get_Struct(argv[1], gsl_matrix, R);
    x = gsl_vector_alloc(Q->size1);
    retval = Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
    break;
  case 4:
    CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
    CHECK_VECTOR(argv[3]);
    Data_Get_Struct(argv[0], gsl_matrix, Q);
    Data_Get_Struct(argv[1], gsl_matrix, R);
    Data_Get_Struct(argv[3], gsl_vector, x);
    retval = argv[3];
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 3 or 4)", argc);
    break;
  }
  switch (flag) {
  case LINALG_QR_DECOMP:
    if (CLASS_OF(argv[0]) != cgsl_matrix_Q) 
      rb_raise(rb_eTypeError, "not a Q matrix");
    if (CLASS_OF(argv[1]) != cgsl_matrix_R) 
      rb_raise(rb_eTypeError, "not a R matrix");
    fsolve = &gsl_linalg_QR_QRsolve;
    break;
#ifdef GSL_1_6_LATER
  case LINALG_LQ_DECOMP:
    /*    if (CLASS_OF(argv[0]) != cgsl_matrix_L) 
      rb_raise(rb_eTypeError, "not a L matrix");
    if (CLASS_OF(argv[1]) != cgsl_matrix_Q) 
    rb_raise(rb_eTypeError, "not a Q matrix");*/
    fsolve = &gsl_linalg_LQ_LQsolve;
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  if (TYPE(argv[2]) == T_ARRAY) {
    b = make_cvector_from_rarray(argv[2]);
    flagb = 1;
  } else {
    CHECK_VECTOR(argv[2]);
    Data_Get_Struct(argv[2], gsl_vector, b);
  }
  (*fsolve)(Q, R, b, x);
  if (flagb == 1) gsl_vector_free(b);
  return retval;
}

/*****/
static VALUE rb_gsl_linalg_QRLQ_RLsolve(int argc, VALUE *argv, VALUE obj,
					  int flag)
{
  gsl_matrix *QR = NULL, *mtmp;
  gsl_vector *b = NULL, *x = NULL, *tau = NULL;
  size_t istart;
  int (*fsolve)(const gsl_matrix*, const gsl_vector*, gsl_vector *);
  int flagb = 0, flagq = 0;
  VALUE omatrix,retval;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError,  "too few arguments");
    omatrix = argv[0];
    istart = 1;
    break;
  default:
    omatrix = obj;
    istart = 0;
    break;
  }
  CHECK_MATRIX(omatrix);
  Data_Get_Struct(omatrix, gsl_matrix, mtmp);
  switch (argc - istart) {
  case 1:
    x = gsl_vector_alloc(mtmp->size1);
    retval = Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
    break;
  case 2:
    Data_Get_Struct(argv[istart+1], gsl_vector, x);
    retval = argv[istart+1];
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 3 or 4)", argc);
    break;
  }
  QR = mtmp; flagq = 0;
  switch (flag) {
  case LINALG_QR_RSOLVE:
    if (CLASS_OF(omatrix) != cgsl_matrix_QR) {
      QR = make_matrix_clone(mtmp);
      tau = gsl_vector_alloc(QR->size1);
      gsl_linalg_QR_decomp(QR, tau);
      flagq = 1;
    }
    fsolve = &gsl_linalg_QR_Rsolve;
    break;
  case LINALG_R_SOLVE:
    if (CLASS_OF(omatrix) != cgsl_matrix_QR) {
      QR = make_matrix_clone(mtmp);
      tau = gsl_vector_alloc(QR->size1);
      gsl_linalg_QR_decomp(QR, tau);
      flagq = 1;
    }
    fsolve = &gsl_linalg_R_solve;
    break;
#ifdef GSL_1_6_LATER
  case LINALG_LQ_LSOLVE:
    if (CLASS_OF(omatrix) != cgsl_matrix_LQ) {
      QR = make_matrix_clone(mtmp);
      tau = gsl_vector_alloc(QR->size1);
      gsl_linalg_LQ_decomp(QR, tau);
      flagq = 1;
    }
    fsolve = &gsl_linalg_LQ_Lsolve_T;
    break;
  case LINALG_L_SOLVE:
    if (CLASS_OF(omatrix) != cgsl_matrix_LQ) {
      QR = make_matrix_clone(mtmp);
      tau = gsl_vector_alloc(QR->size1);
      gsl_linalg_LQ_decomp(QR, tau);
      flagq = 1;
    }
    fsolve = &gsl_linalg_L_solve_T;
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  if (TYPE(argv[istart]) == T_ARRAY) {
    b = make_cvector_from_rarray(argv[istart]);
    flagb = 1;
  } else {
    CHECK_VECTOR(argv[istart]);
    Data_Get_Struct(argv[istart], gsl_vector, b);
  }
  (*fsolve)(QR, b, x);
  if (flagb == 1) gsl_vector_free(b);
  if (flagq == 1) {
    gsl_matrix_free(QR);
    gsl_vector_free(tau);
  }
  return retval;
}

static VALUE rb_gsl_linalg_QRLQ_RLsvx(int argc, VALUE *argv, VALUE obj,
					  int flag)
{
  gsl_matrix *QR = NULL, *mtmp;
  gsl_vector *x = NULL, *tau = NULL;
  size_t istart;
  int (*fsolve)(const gsl_matrix*, gsl_vector *);
  int flagq = 0;
  VALUE omatrix,retval;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError,  "too few arguments");
    omatrix = argv[0];
    istart = 1;
    break;
  default:
    omatrix = obj;
    istart = 0;
    break;
  }
  CHECK_MATRIX(omatrix);
  Data_Get_Struct(omatrix, gsl_matrix, mtmp);
  switch (argc - istart) {
  case 0:
    x = gsl_vector_alloc(mtmp->size1);
    retval = Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
    break;
  case 1:
    Data_Get_Struct(argv[istart+1], gsl_vector, x);
    retval = argv[istart+1];
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 3 or 4)", argc);
    break;
  }
  QR = mtmp; flagq = 0;
  switch (flag) {
  case LINALG_QR_RSVX:
    if (CLASS_OF(omatrix) != cgsl_matrix_QR) {
      QR = make_matrix_clone(mtmp);
      tau = gsl_vector_alloc(QR->size1);
      gsl_linalg_QR_decomp(QR, tau);
      flagq = 1;
    }
    fsolve = &gsl_linalg_QR_Rsvx;
    break;
    /*
  case LINALG_R_SVX:
    if (CLASS_OF(omatrix) != cgsl_matrix_QR) {
      QR = make_matrix_clone(mtmp);
      tau = gsl_vector_alloc(QR->size1);
      gsl_linalg_QR_decomp(QR, tau);
      flagq = 1;
    }
    fsolve = &gsl_linalg_R_svx;
    break;
    */
#ifdef GSL_1_6_LATER
  case LINALG_LQ_LSVX:
    if (CLASS_OF(omatrix) != cgsl_matrix_LQ) {
      QR = make_matrix_clone(mtmp);
      tau = gsl_vector_alloc(QR->size1);
      gsl_linalg_LQ_decomp(QR, tau);
      flagq = 1;
    }
    fsolve = &gsl_linalg_LQ_Lsvx_T;
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  (*fsolve)(QR, x);
  if (flagq == 1) {
    gsl_matrix_free(QR);
    gsl_vector_free(tau);
  }
  return retval;
}

static VALUE rb_gsl_linalg_QR_Rsolve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQ_RLsolve(argc, argv, obj, LINALG_QR_RSOLVE);
}

static VALUE rb_gsl_linalg_QR_Rsvx(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQ_RLsvx(argc, argv, obj, LINALG_QR_RSVX);
}

static VALUE rb_gsl_linalg_R_solve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQ_RLsolve(argc, argv, obj, LINALG_R_SOLVE);
}

/* singleton */
static VALUE rb_gsl_linalg_QR_QRsolve(int argc, VALUE *argv, VALUE obj,
					  int flag)
{
  return rb_gsl_linalg_QRLQ_QRLQsolve(argc, argv, obj, LINALG_QR_DECOMP);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_LQ_Lsolve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQ_RLsolve(argc, argv, obj, LINALG_LQ_LSOLVE);
}

static VALUE rb_gsl_linalg_LQ_Lsvx(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQ_RLsvx(argc, argv, obj, LINALG_LQ_LSVX);
}

static VALUE rb_gsl_linalg_L_solve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQ_RLsolve(argc, argv, obj, LINALG_L_SOLVE);
}

/* singleton */
static VALUE rb_gsl_linalg_LQ_LQsolve(int argc, VALUE *argv, VALUE obj,
					  int flag)
{
  return rb_gsl_linalg_QRLQ_QRLQsolve(argc, argv, obj, LINALG_LQ_DECOMP);
}
#endif

static VALUE rb_gsl_linalg_QRLQ_update(VALUE obj, VALUE qq, VALUE rr, VALUE ww,
				     VALUE vv, int flag)
{
  gsl_matrix *Q = NULL, *R = NULL;
  gsl_vector *w = NULL, *v = NULL;
  int status;
  CHECK_MATRIX(qq); CHECK_MATRIX(rr);
  CHECK_VECTOR(ww); CHECK_VECTOR(vv);
  Data_Get_Struct(qq, gsl_matrix, Q);
  Data_Get_Struct(rr, gsl_matrix, R);
  Data_Get_Struct(ww, gsl_vector, w);
  Data_Get_Struct(vv, gsl_vector, v);
  switch (flag) {
  case LINALG_QR_DECOMP:
    status = gsl_linalg_QR_update(Q, R, w, v);
    break;
#ifdef GSL_1_6_LATER
  case LINALG_LQ_DECOMP:
    status = gsl_linalg_LQ_update(Q, R, w, v);
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  return INT2FIX(status);
}

/* singleton */
static VALUE rb_gsl_linalg_QR_update(VALUE obj, VALUE qq, VALUE rr, VALUE ww,
				     VALUE vv)
{
  return rb_gsl_linalg_QRLQ_update(obj, qq, rr, ww, vv, LINALG_QR_DECOMP);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_LQ_update(VALUE obj, VALUE qq, VALUE rr, VALUE ww,
				     VALUE vv)
{
  return rb_gsl_linalg_QRLQ_update(obj, qq, rr, ww, vv, LINALG_LQ_DECOMP);
}
#endif

/******/
enum {
  LINALG_QRPT,
  LINALG_PTLQ,
};

static VALUE rb_gsl_linalg_QRLQPT_decomp(int argc, VALUE *argv, VALUE obj, int flag)
{
  gsl_matrix *A = NULL, *QR = NULL;
  gsl_vector *tau = NULL, *norm = NULL;
  gsl_permutation *p = NULL;
  int signum;
  size_t size0;
  VALUE vtau, vp, vA, vQR;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
    vA = argv[0];
    break;
  default:
    vA = obj;
    break;
  }
  CHECK_MATRIX(vA);
  Data_Get_Struct(vA, gsl_matrix, A);
  QR = make_matrix_clone(A);
  size0 = GSL_MIN(A->size1, A->size2);
  tau = gsl_vector_alloc(size0);
  p = gsl_permutation_alloc(size0);
  norm = gsl_vector_alloc(size0);
  switch (flag) {
  case LINALG_QRPT:
    vQR = Data_Wrap_Struct(cgsl_matrix_QRPT, 0, gsl_matrix_free, QR);
    vtau = Data_Wrap_Struct(cgsl_vector_tau, 0, gsl_vector_free, tau);
    vp = Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p);
    gsl_linalg_QRPT_decomp(QR, tau, p, &signum, norm);
    break;
#ifdef GSL_1_6_LATER
  case LINALG_PTLQ:
    vQR = Data_Wrap_Struct(cgsl_matrix_PTLQ, 0, gsl_matrix_free, QR);
    vtau = Data_Wrap_Struct(cgsl_vector_tau, 0, gsl_vector_free, tau);
    vp = Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p);
    gsl_linalg_PTLQ_decomp(QR, tau, p, &signum, norm);
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  gsl_vector_free(norm);
  return rb_ary_new3(4, vQR, vtau, vp, INT2FIX(signum));
}

static VALUE rb_gsl_linalg_QRLQPT_decomp_bang(int argc, VALUE *argv, VALUE obj, int flag)
{
  gsl_matrix *A = NULL;
  gsl_vector *tau = NULL, *norm = NULL;
  gsl_permutation *p = NULL;
  int signum;
  size_t size0;
  VALUE vtau, vp, vA;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
    vA = argv[0];
    break;
  default:
    vA = obj;
    break;
  }
  CHECK_MATRIX(vA);
  Data_Get_Struct(vA, gsl_matrix, A);
  size0 = GSL_MIN(A->size1, A->size2);
  tau = gsl_vector_alloc(size0);
  p = gsl_permutation_alloc(size0);
  norm = gsl_vector_alloc(size0);
  switch (flag) {
  case LINALG_QRPT:
    RBASIC(vA)->klass = cgsl_matrix_QRPT;
    vtau = Data_Wrap_Struct(cgsl_vector_tau, 0, gsl_vector_free, tau);
    vp = Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p);
    gsl_linalg_QRPT_decomp(A, tau, p, &signum, norm);
    break;
#ifdef GSL_1_6_LATER
  case LINALG_PTLQ:
    RBASIC(vA)->klass = cgsl_matrix_PTLQ;
    vtau = Data_Wrap_Struct(cgsl_vector_tau, 0, gsl_vector_free, tau);
    vp = Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p);
    gsl_linalg_PTLQ_decomp(A, tau, p, &signum, norm);
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  gsl_vector_free(norm);
  return rb_ary_new3(3, vtau, vp, INT2FIX(signum));
}

static VALUE rb_gsl_linalg_QRPT_decomp(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_decomp(argc, argv, obj, LINALG_QRPT);
}

static VALUE rb_gsl_linalg_QRPT_decomp_bang(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_decomp_bang(argc, argv, obj, LINALG_QRPT);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_PTLQ_decomp(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_decomp(argc, argv, obj, LINALG_PTLQ);
}

static VALUE rb_gsl_linalg_PTLQ_decomp_bang(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_decomp_bang(argc, argv, obj, LINALG_PTLQ);
}
#endif

static VALUE rb_gsl_linalg_QRLQPT_decomp2(int argc, VALUE *argv, VALUE obj,int flag)
{
  gsl_matrix *A = NULL, *Q = NULL, *R = NULL;
  gsl_vector *tau = NULL, *norm = NULL;
  gsl_permutation *p = NULL;
  int signum;
  size_t size0;
  VALUE vtau, vp, vA, vQ, vR;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments");
    vA = argv[0];
    break;
  default:
    if (argc != 0) rb_raise(rb_eArgError, "wrong number of arguments");
    vA = obj;
    break;
  }
  CHECK_MATRIX(vA);
  Data_Get_Struct(vA, gsl_matrix, A);
  Q = gsl_matrix_alloc(A->size1, A->size2);
  R = gsl_matrix_alloc(A->size1, A->size2);
  size0 = GSL_MIN(A->size1, A->size2);
  tau = gsl_vector_alloc(size0);
  p = gsl_permutation_alloc(size0);
  norm = gsl_vector_alloc(size0);
  /*  vQ = Data_Wrap_Struct(cgsl_matrix_Q, 0, gsl_matrix_free, Q);
      vR = Data_Wrap_Struct(cgsl_matrix_R, 0, gsl_matrix_free, R);*/
  vtau = Data_Wrap_Struct(cgsl_vector_tau, 0, gsl_vector_free, tau);
  vp = Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p);
  switch (flag) {
  case LINALG_QRPT:
    vQ = Data_Wrap_Struct(cgsl_matrix_Q, 0, gsl_matrix_free, Q);
    vR = Data_Wrap_Struct(cgsl_matrix_R, 0, gsl_matrix_free, R);
    gsl_linalg_QRPT_decomp2(A, Q, R, tau, p, &signum, norm);
    break;
#ifdef GSL_1_6_LATER
  case LINALG_PTLQ:
    vR = Data_Wrap_Struct(cgsl_matrix_L, 0, gsl_matrix_free, R);
    vQ = Data_Wrap_Struct(cgsl_matrix_Q, 0, gsl_matrix_free, Q);
    gsl_linalg_PTLQ_decomp2(A, Q, R, tau, p, &signum, norm);
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
  }
  gsl_vector_free(norm);
  return rb_ary_new3(5, vQ, vR, vtau, vp, INT2FIX(signum));
}

static VALUE rb_gsl_linalg_QRPT_decomp2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_decomp2(argc, argv, obj, LINALG_QRPT);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_PTLQ_decomp2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_decomp2(argc, argv, obj, LINALG_PTLQ);
}
#endif

#ifdef GSL_1_6_LATER
int gsl_linalg_PTLQ_solve_T(const gsl_matrix * QR, const gsl_vector * tau,
			    const gsl_permutation * p, const gsl_vector * b, 
			    gsl_vector * x);
int gsl_linalg_PTLQ_svx_T(const gsl_matrix * LQ,
                         const gsl_vector * tau,
                         const gsl_permutation * p,
                         gsl_vector * x);
int gsl_linalg_PTLQ_LQsolve_T (const gsl_matrix * Q, const gsl_matrix * L,
                           const gsl_permutation * p,
                           const gsl_vector * b,
                           gsl_vector * x);
int gsl_linalg_PTLQ_Lsolve_T (const gsl_matrix * LQ,
                        const gsl_permutation * p,
                        const gsl_vector * b,
			  gsl_vector * x);
int gsl_linalg_PTLQ_Lsvx_T (const gsl_matrix * LQ,
                        const gsl_permutation * p,
                        gsl_vector * x);
#endif

static VALUE rb_gsl_linalg_QRLQPT_solve(int argc, VALUE *argv, VALUE obj, int flag)
{
  gsl_matrix *QR = NULL, *A = NULL;
  gsl_vector *tau = NULL, *b = NULL, *x = NULL, *norm = NULL;
  gsl_permutation *p = NULL;
  int signum, itmp, flagb = 0, flagq = 0;
  VALUE vtmp, klass;
  size_t size0;
  int (*fdecomp)(gsl_matrix*, gsl_vector*, gsl_permutation*, int *, gsl_vector*);
  int (*fsolve)(const gsl_matrix*, const gsl_vector*, const gsl_permutation*, 
		const gsl_vector*, gsl_vector *);
  switch (flag) {
  case LINALG_QRPT:
    klass = cgsl_matrix_QRPT;
    fdecomp = &gsl_linalg_QRPT_decomp;
    fsolve = &gsl_linalg_QRPT_solve;
    break;
#ifdef GSL_1_6_LATER
  case LINALG_PTLQ:
    klass = cgsl_matrix_PTLQ;
    fdecomp = &gsl_linalg_PTLQ_decomp;
    fsolve = &gsl_linalg_PTLQ_solve_T;
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments.");
    vtmp = argv[0];
    itmp = 1;
    break;
  default:
    vtmp = obj;
    itmp = 0;
    break;
  }
  CHECK_MATRIX(vtmp);
  if (CLASS_OF(vtmp) == klass) {
    if (argc-itmp != 3) rb_raise(rb_eArgError, 
				 "wrong number of arguments (%d for %d)", 
				 argc, 4-itmp);
    CHECK_VECTOR(argv[itmp]);
    if (CLASS_OF(argv[itmp]) != cgsl_vector_tau) 
      rb_raise(rb_eTypeError, "not a tau vector");
    CHECK_PERMUTATION(argv[itmp+1]);
    Data_Get_Struct(argv[itmp], gsl_vector, tau);
    Data_Get_Struct(argv[itmp+1], gsl_permutation, p);
    Data_Get_Struct(vtmp, gsl_matrix, QR);
    size0 = GSL_MIN(QR->size1, QR->size2);
    itmp += 2;
  } else {
    if (argc-itmp != 1) rb_raise(rb_eArgError, 
			    "wrong number of arguments (%d for %d)", argc, 2-itmp);
    Data_Get_Struct(vtmp, gsl_matrix, A);
    QR = make_matrix_clone(A);
    size0 = GSL_MIN(QR->size1, QR->size2);
    flagq = 1;
    p = gsl_permutation_alloc(size0);
    tau = gsl_vector_alloc(size0);
  }
  norm = gsl_vector_alloc(size0);
  if (TYPE(argv[itmp]) == T_ARRAY) {
    b = make_cvector_from_rarray(argv[itmp]);
    flagb = 1;
  } else {
    CHECK_VECTOR(argv[itmp]);
    Data_Get_Struct(argv[itmp], gsl_vector, b);
  }
  x = gsl_vector_alloc(b->size);
  if (flagq == 1) (*fdecomp)(QR, tau, p, &signum, norm);
  (*fsolve)(QR, tau, p, b, x);
  if (flagb == 1) gsl_vector_free(b);
  if (flagq == 1) {
    gsl_matrix_free(QR);
    gsl_permutation_free(p);
    gsl_vector_free(tau);
    gsl_vector_free(norm);
  }
  return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
}

static VALUE rb_gsl_linalg_QRPT_solve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_solve(argc, argv, obj, LINALG_QRPT);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_PTLQ_solve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_solve(argc, argv, obj, LINALG_PTLQ);
}
#endif

static VALUE rb_gsl_linalg_QRLQPT_svx(int argc, VALUE *argv, VALUE obj, int flag)
{
  gsl_matrix *QR = NULL, *A = NULL;
  gsl_vector *tau = NULL, *b = NULL, *norm = NULL;
  gsl_permutation *p = NULL;
  int signum, itmp, flagq = 0;
  VALUE vtmp, klass;
  size_t size0;
  int (*fdecomp)(gsl_matrix*, gsl_vector*, gsl_permutation*, int *, gsl_vector*);
  int (*fsvx)(const gsl_matrix*, const gsl_vector*, const gsl_permutation*, 
	      gsl_vector *);
  switch (flag) {
  case LINALG_QRPT:
    klass = cgsl_matrix_QRPT;
    fdecomp = &gsl_linalg_QRPT_decomp;
    fsvx = &gsl_linalg_QRPT_svx;
    break;
#ifdef GSL_1_6_LATER
  case LINALG_PTLQ:
    klass = cgsl_matrix_PTLQ;
    fdecomp = &gsl_linalg_PTLQ_decomp;
    fsvx = &gsl_linalg_PTLQ_svx_T;
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
    vtmp = argv[0];
    itmp = 1;
    break;
  default:
    vtmp = obj;
    itmp = 0;
    break;
  }
  CHECK_MATRIX(vtmp);
  if (CLASS_OF(vtmp) == klass) {
    if (argc-itmp != 3) rb_raise(rb_eArgError, 
				 "wrong number of arguments (%d for %d)", 
				 argc, 3+itmp);
    CHECK_VECTOR(argv[itmp]);
    if (CLASS_OF(argv[itmp]) != cgsl_vector_tau) 
      rb_raise(rb_eTypeError, "not a tau vector");
    CHECK_PERMUTATION(argv[itmp+1]);
    Data_Get_Struct(argv[itmp], gsl_vector, tau);
    Data_Get_Struct(argv[itmp+1], gsl_permutation, p);
    Data_Get_Struct(vtmp, gsl_matrix, QR);
    size0 = GSL_MIN(QR->size1, QR->size2);
    itmp += 2;
  } else {
    if (argc-itmp != 1) rb_raise(rb_eArgError, 
			    "wrong number of arguments (%d for %d)", argc, 2+itmp);
    Data_Get_Struct(vtmp, gsl_matrix, A);
    QR = make_matrix_clone(A);
    size0 = GSL_MIN(QR->size1, QR->size2);
    flagq = 1;
    p = gsl_permutation_alloc(size0);
    tau = gsl_vector_alloc(size0);
  }
  norm = gsl_vector_alloc(size0);
  CHECK_VECTOR(argv[itmp]);
  Data_Get_Struct(argv[itmp], gsl_vector, b);
  if (flagq == 1) (*fdecomp)(QR, tau, p, &signum, norm);
  (*fsvx)(QR, tau, p, b);
  if (flagq == 1) {
    gsl_matrix_free(QR);
    gsl_permutation_free(p);
    gsl_vector_free(tau);
    gsl_vector_free(norm);
  }
  return argv[itmp];
}

static VALUE rb_gsl_linalg_QRPT_svx(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_svx(argc, argv, obj, LINALG_QRPT);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_PTLQ_svx(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_svx(argc, argv, obj, LINALG_PTLQ);
}
#endif

/* singleton */
static VALUE rb_gsl_linalg_QRLQPT_QRLQsolve(VALUE obj, VALUE qq, VALUE rr, 
					    VALUE pp, VALUE bb, int flag)
{
  gsl_matrix *Q = NULL, *R = NULL;
  gsl_vector *b = NULL, *x = NULL;
  gsl_permutation *p = NULL;
  int flagb = 0;
  int (*fsolve)(const gsl_matrix*, const gsl_matrix*, const gsl_permutation*, 
		const gsl_vector*, gsl_vector*);
  switch (flag) {
  case LINALG_QRPT:
    if (CLASS_OF(qq) != cgsl_matrix_Q) rb_raise(rb_eTypeError, "not a Q matrix");
    if (CLASS_OF(rr) != cgsl_matrix_R) rb_raise(rb_eTypeError, "not a R matrix");
    fsolve = &gsl_linalg_QRPT_QRsolve;
    break;
#ifdef GSL_1_6_LATER
  case LINALG_PTLQ:
    if (CLASS_OF(qq) != cgsl_matrix_Q) rb_raise(rb_eTypeError, "not a Q matrix");
    if (CLASS_OF(rr) != cgsl_matrix_L) rb_raise(rb_eTypeError, "not a L matrix");
    fsolve = &gsl_linalg_PTLQ_LQsolve_T;
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  if (TYPE(bb) == T_ARRAY) {
    b = make_cvector_from_rarray(bb);
    flagb = 1;
  } else {
    CHECK_VECTOR(bb);
    Data_Get_Struct(bb, gsl_vector, b);
  }
  CHECK_PERMUTATION(pp);
  Data_Get_Struct(qq, gsl_matrix, Q);
  Data_Get_Struct(rr, gsl_matrix, R);
  Data_Get_Struct(pp, gsl_permutation, p);
  x = gsl_vector_alloc(b->size);
  (*fsolve)(Q, R, p, b, x);
  if (flagb == 1) gsl_vector_free(b);
  return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
}

static VALUE rb_gsl_linalg_QRPT_QRsolve(VALUE obj, VALUE qq, VALUE rr, 
					VALUE pp, VALUE bb)
{
  return rb_gsl_linalg_QRLQPT_QRLQsolve(obj, qq, rr, pp, bb, LINALG_QRPT);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_PTLQ_LQsolve(VALUE obj, VALUE qq, VALUE rr, 
					VALUE pp, VALUE bb)
{
  return rb_gsl_linalg_QRLQPT_QRLQsolve(obj, qq, rr, pp, bb, LINALG_PTLQ);
}
#endif

/* singleton */
static VALUE rb_gsl_linalg_QRLQPT_update(VALUE obj, VALUE qq, VALUE rr, 
				       VALUE pp, VALUE ww, VALUE vv, int flag)
{
  gsl_matrix *Q = NULL, *R = NULL;
  gsl_vector *w = NULL, *v = NULL;
  gsl_permutation *p = NULL;
  switch (flag) {
  case LINALG_QRPT:
    if (CLASS_OF(qq) != cgsl_matrix_Q) rb_raise(rb_eTypeError, "not a Q matrix");
    if (CLASS_OF(rr) != cgsl_matrix_R) rb_raise(rb_eTypeError, "not a R matrix");
    break;
#ifdef GSL_1_6_LATER
  case LINALG_PTLQ:
    if (CLASS_OF(qq) != cgsl_matrix_Q) rb_raise(rb_eTypeError, "not a Q matrix");
    if (CLASS_OF(rr) != cgsl_matrix_L) rb_raise(rb_eTypeError, "not a L matrix");
    break;
#endif
  }
  CHECK_PERMUTATION(pp);
  Data_Get_Struct(qq, gsl_matrix, Q);
  Data_Get_Struct(rr, gsl_matrix, R);
  Data_Get_Struct(pp, gsl_permutation, p);
  Data_Get_Struct(ww, gsl_vector, w);
  Data_Get_Struct(vv, gsl_vector, v);
  switch (flag) {
  case LINALG_QRPT:
    gsl_linalg_QRPT_update(Q, R, p, w, v);
    break;
#ifdef GSL_1_6_LATER
  case LINALG_PTLQ:
    gsl_linalg_PTLQ_update(Q, R, p, w, v);
    break;
#endif
  }
  return obj;
}

static VALUE rb_gsl_linalg_QRPT_update(VALUE obj, VALUE qq, VALUE rr, 
				       VALUE pp, VALUE ww, VALUE vv)
{
  return rb_gsl_linalg_QRLQPT_update(obj, qq, rr, pp, ww, vv, LINALG_QRPT);
}

#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_PTLQ_update(VALUE obj, VALUE qq, VALUE rr, 
				       VALUE pp, VALUE ww, VALUE vv)
{
  return rb_gsl_linalg_QRLQPT_update(obj, qq, rr, pp, ww, vv, LINALG_PTLQ);
}
#endif

static VALUE rb_gsl_linalg_QRLQPT_RLsolve(int argc, VALUE *argv, VALUE obj, int flag)
{
  gsl_matrix *QR = NULL;
  gsl_vector *b = NULL, *x = NULL;
  gsl_permutation *p = NULL;
  int itmp, flagb = 0;
  VALUE vtmp, klass;
  int (*fsolve)(const gsl_matrix*, const gsl_permutation*, const gsl_vector*, 
		gsl_vector*);
  switch (flag) {
  case LINALG_QRPT:
    klass = cgsl_matrix_QRPT;
    fsolve = &gsl_linalg_QRPT_Rsolve;
    break;
#ifdef GSL_1_6_LATER
  case LINALG_PTLQ:
    klass = cgsl_matrix_PTLQ;
    fsolve = &gsl_linalg_PTLQ_Lsolve_T;
    break;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
    vtmp = argv[0];
    itmp = 1;
    break;
  default:
    vtmp = obj;
    itmp = 0;
    break;
  }
  if (argc-itmp != 2)
    rb_raise(rb_eArgError, "wrong number of argument (%d for %d)", argc, 2+itmp);
  CHECK_MATRIX(vtmp);
  if (CLASS_OF(vtmp) != klass) {
    rb_raise(rb_eArgError, "not a QR matrix");
  }
  CHECK_PERMUTATION(argv[itmp]);
  Data_Get_Struct(argv[itmp], gsl_permutation, p);
  Data_Get_Struct(vtmp, gsl_matrix, QR);
  itmp++;
  if (TYPE(argv[itmp]) == T_ARRAY) {
    b = make_cvector_from_rarray(argv[itmp]);
    flagb = 1;
  } else {
    CHECK_VECTOR(argv[itmp]);
    Data_Get_Struct(argv[itmp], gsl_vector, b);
  }
  x = gsl_vector_alloc(b->size);
  (*fsolve)(QR, p, b, x);
  if (flagb == 1) gsl_vector_free(b);
  return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
}

static VALUE rb_gsl_linalg_QRPT_Rsolve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_RLsolve(argc, argv, obj, LINALG_QRPT);
}
#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_PTLQ_Lsolve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_RLsolve(argc, argv, obj, LINALG_PTLQ);
}
#endif

static VALUE rb_gsl_linalg_QRLQPT_RLsvx(int argc, VALUE *argv, VALUE obj, int flag)
{
  gsl_matrix *QR = NULL;
  gsl_vector *b = NULL;
  gsl_permutation *p = NULL;
  int itmp;
  VALUE vtmp, klass;
  int (*fsvx)(const gsl_matrix*, const gsl_permutation*, gsl_vector*);
  switch (flag) {
  case LINALG_QRPT:
    klass = cgsl_matrix_QRPT;
    fsvx = &gsl_linalg_QRPT_Rsvx;
    break;
#ifdef GSL_1_6_LATER
  case LINALG_PTLQ:
    klass = cgsl_matrix_PTLQ;
    fsvx = &gsl_linalg_PTLQ_Lsvx_T;
#endif
  default:
    rb_raise(rb_eRuntimeError, "unknown operation");
    break;
  }
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
    vtmp = argv[0];
    itmp = 1;
    break;
  default:
    vtmp = obj;
    itmp = 0;
    break;
  }
  if (argc-itmp != 2)
    rb_raise(rb_eArgError, "wrong number of argument (%d for %d)", argc, 2+itmp);
  CHECK_MATRIX(vtmp);
  if (CLASS_OF(vtmp) != klass) {
    rb_raise(rb_eArgError, "not a QR matrix");
  }
  CHECK_PERMUTATION(argv[itmp]);
  Data_Get_Struct(argv[itmp], gsl_permutation, p);
  Data_Get_Struct(vtmp, gsl_matrix, QR);
  itmp++;
  if (TYPE(argv[itmp]) == T_ARRAY) {
    b = make_cvector_from_rarray(argv[itmp]);
  } else {
    CHECK_VECTOR(argv[itmp]);
    Data_Get_Struct(argv[itmp], gsl_vector, b);
  }
  (*fsvx)(QR, p, b);
  return argv[itmp];
}

static VALUE rb_gsl_linalg_QRPT_Rsvx(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_RLsvx(argc, argv, obj, LINALG_QRPT);
}
#ifdef GSL_1_6_LATER
static VALUE rb_gsl_linalg_PTLQ_Lsvx(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_linalg_QRLQPT_RLsvx(argc, argv, obj, LINALG_PTLQ);
}
#endif

/*******/
#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_SV_decomp_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *A;
  gsl_matrix_view uv, vv;
  gsl_vector_view sv;
  gsl_vector *work;
  VALUE u, s, v;
  int shape[2];
  GetNArray(argv[0], A);
  shape[0] = A->shape[0];
  shape[1] = A->shape[0];
  u = na_make_object(NA_DFLOAT, 2, A->shape, CLASS_OF(argv[0]));
  v = na_make_object(NA_DFLOAT, 2, shape, CLASS_OF(argv[0]));
  s = na_make_object(NA_DFLOAT, 1, &(shape[0]), cNVector);
  uv = gsl_matrix_view_array(NA_PTR_TYPE(u,double*), A->shape[1], A->shape[0]);
  vv = gsl_matrix_view_array(NA_PTR_TYPE(v,double*), shape[1], shape[0]);
  sv = gsl_vector_view_array(NA_PTR_TYPE(s,double*), shape[0]);
  work = gsl_vector_alloc(shape[0]);
  memcpy(NA_PTR_TYPE(u,double*), (double*)A->ptr, sizeof(double)*A->total);
  gsl_linalg_SV_decomp(&uv.matrix, &vv.matrix, &sv.vector, work);
  gsl_vector_free(work);
  return rb_ary_new3(3, u, v, s);
}

static VALUE rb_gsl_linalg_SV_decomp_jacobi_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *A;
  gsl_matrix_view uv, vv;
  gsl_vector_view sv;
  VALUE u, s, v;
  int shape[2];
  GetNArray(argv[0], A);
  shape[0] = A->shape[0];
  shape[1] = A->shape[0];
  u = na_make_object(NA_DFLOAT, 2, A->shape, CLASS_OF(argv[0]));
  v = na_make_object(NA_DFLOAT, 2, shape, CLASS_OF(argv[0]));
  s = na_make_object(NA_DFLOAT, 1, &(shape[0]), cNVector);
  uv = gsl_matrix_view_array(NA_PTR_TYPE(u,double*), A->shape[1], A->shape[0]);
  vv = gsl_matrix_view_array(NA_PTR_TYPE(v,double*), shape[1], shape[0]);
  sv = gsl_vector_view_array(NA_PTR_TYPE(s,double*), shape[0]);
  memcpy(NA_PTR_TYPE(u,double*), (double*)A->ptr, sizeof(double)*A->total);
  gsl_linalg_SV_decomp_jacobi(&uv.matrix, &vv.matrix, &sv.vector);
  return rb_ary_new3(3, u, v, s);
}

static VALUE rb_gsl_linalg_SV_solve_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *A;
  gsl_matrix_view uv, vv;
  gsl_vector_view sv, bv, xv;
  VALUE x;
  if (argc != 4)
    rb_raise(rb_eArgError, "Usage: SV.solve(u, v, s, b)");
  GetNArray(argv[0], A);
  uv = gsl_matrix_view_array(NA_PTR_TYPE(argv[0],double*), A->shape[1], A->shape[0]);
  vv = gsl_matrix_view_array(NA_PTR_TYPE(argv[1],double*), A->shape[0], A->shape[0]);
  sv = gsl_vector_view_array(NA_PTR_TYPE(argv[2],double*), A->shape[0]);
  bv = gsl_vector_view_array(NA_PTR_TYPE(argv[3],double*), A->shape[0]);
  x = na_make_object(NA_DFLOAT, 1, &(A->shape[0]), CLASS_OF(argv[3]));
  xv = gsl_vector_view_array(NA_PTR_TYPE(x,double*), A->shape[0]);
  gsl_linalg_SV_solve(&uv.matrix, &vv.matrix, &sv.vector, &bv.vector, &xv.vector);
  return x;
}

#endif

static VALUE rb_gsl_linalg_SV_decomp(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *V = NULL, *U = NULL;
  gsl_vector *w = NULL, *S = NULL;
  int flag = 1;
  VALUE vs, vv, vu;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    switch (argc) {
    case 2:
      CHECK_VECTOR(argv[1]);
      Data_Get_Struct(argv[1], gsl_vector, w);
      flag = 0;
      /* no break, do next */
    case 1:
#ifdef HAVE_NARRAY_H
      if (NA_IsNArray(argv[0]))
	return rb_gsl_linalg_SV_decomp_narray(argc, argv, obj);
#endif
      CHECK_MATRIX(argv[0]);
      Data_Get_Struct(argv[0], gsl_matrix, A);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 0:
      /* do nothing */
      break;
    case 1:
      CHECK_VECTOR(argv[0]);
      Data_Get_Struct(argv[0], gsl_vector, w);
      flag = 0;
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
      break;
    }
    Data_Get_Struct(obj, gsl_matrix, A);
    break;
  }
  U = make_matrix_clone(A);
  S = gsl_vector_alloc(A->size2);   /* see manual p 123 */
  V = gsl_matrix_alloc(A->size2, A->size2);
  if (flag == 1) w = gsl_vector_alloc(A->size2);
  gsl_linalg_SV_decomp(U, V, S, w);
  if (flag == 1) gsl_vector_free(w);
  vu = Data_Wrap_Struct(cgsl_matrix_U, 0, gsl_matrix_free, U);
  vv = Data_Wrap_Struct(cgsl_matrix_V, 0, gsl_matrix_free, V);
  vs = Data_Wrap_Struct(cgsl_vector_S, 0, gsl_vector_free, S);
  return rb_ary_new3(3, vu, vv, vs);
}

static VALUE rb_gsl_linalg_SV_decomp_mod(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *V = NULL, *U = NULL, *X = NULL;
  gsl_vector *w = NULL, *S = NULL;
  VALUE vs, vv, vu;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError,
			    "wrong number of argument (%d for 1)", argc);
    CHECK_MATRIX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix, A);
    break;
  }
  U = make_matrix_clone(A);
  S = gsl_vector_alloc(A->size2);   /* see manual p 123 */
  V = gsl_matrix_alloc(A->size2, A->size2);
  X = gsl_matrix_alloc(A->size2, A->size2);
  w = gsl_vector_alloc(A->size2);
  gsl_linalg_SV_decomp_mod(U, X, V, S, w);
  gsl_vector_free(w);
  gsl_matrix_free(X);
  vu = Data_Wrap_Struct(cgsl_matrix_U, 0, gsl_matrix_free, U);
  vv = Data_Wrap_Struct(cgsl_matrix_V, 0, gsl_matrix_free, V);
  vs = Data_Wrap_Struct(cgsl_vector_S, 0, gsl_vector_free, S);
  return rb_ary_new3(3, vu, vv, vs);
}

static VALUE rb_gsl_linalg_SV_decomp_jacobi(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *V = NULL, *U = NULL;
  gsl_vector *S = NULL;
  VALUE vs, vv, vu;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError,
			    "wrong number of argument (%d for 1)", argc);
#ifdef HAVE_NARRAY_H
      if (NA_IsNArray(argv[0]))
	return rb_gsl_linalg_SV_decomp_jacobi_narray(argc, argv, obj);
#endif
    CHECK_MATRIX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix, A);
    break;
  }
  U = make_matrix_clone(A);
  S = gsl_vector_alloc(A->size2);   /* see manual p 123 */
  V = gsl_matrix_alloc(A->size2, A->size2);
  gsl_linalg_SV_decomp_jacobi(U, V, S);
  vu = Data_Wrap_Struct(cgsl_matrix_U, 0, gsl_matrix_free, U);
  vv = Data_Wrap_Struct(cgsl_matrix_V, 0, gsl_matrix_free, V);
  vs = Data_Wrap_Struct(cgsl_vector_S, 0, gsl_vector_free, S);
  return rb_ary_new3(3, vu, vv, vs);
}

static VALUE rb_gsl_linalg_SV_solve(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *U = NULL, *V = NULL;
  gsl_vector *S = NULL, *b = NULL, *x = NULL;
  int flagb = 0, flagv = 0;
  
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments.");
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0]))
      return rb_gsl_linalg_SV_solve_narray(argc, argv, obj);
#endif

    CHECK_MATRIX(argv[0]);
    if (CLASS_OF(argv[0]) == cgsl_matrix_U) {
      if (argc != 4) rb_raise(rb_eArgError, 
			      "wrong number of arguments (%d for 4)", argc);
      Data_Get_Struct(argv[0], gsl_matrix, U);
      CHECK_MATRIX(argv[1]);
      if (CLASS_OF(argv[1]) != cgsl_matrix_V) 
	rb_raise(rb_eTypeError, "not a V matrix");
      Data_Get_Struct(argv[1], gsl_matrix, V);
      CHECK_VECTOR(argv[2]);
      if (CLASS_OF(argv[2]) != cgsl_vector_S) 
	rb_raise(rb_eTypeError, "not a S vector");
      Data_Get_Struct(argv[2], gsl_vector, S);
      if (TYPE(argv[3]) == T_ARRAY) {
	b = make_cvector_from_rarray(argv[3]);
	flagb = 1;
      } else {
	CHECK_VECTOR(argv[3]);
	Data_Get_Struct(argv[3], gsl_vector, b);
      }
    } else {
      if (argc != 2) rb_raise(rb_eArgError, 
			      "wrong number of arguments (%d for 2)", argc);
      Data_Get_Struct(argv[0], gsl_matrix, A);
      U = make_matrix_clone(A);
      if (TYPE(argv[1]) == T_ARRAY) {
	b = make_cvector_from_rarray(argv[1]);
	flagb = 1;
      } else {
	CHECK_VECTOR(argv[1]);
	Data_Get_Struct(argv[1], gsl_vector, b);
      }
      S = gsl_vector_alloc(A->size2);   /* see manual p 123 */
      V = gsl_matrix_alloc(A->size2, A->size2);
      gsl_linalg_SV_decomp_jacobi(U, V, S);
      flagv = 1;
    }
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, 
			    "wrong number of arguments (%d for 1)", argc);
    Data_Get_Struct(obj, gsl_matrix, A);
    U = make_matrix_clone(A);
    if (TYPE(argv[0]) == T_ARRAY) {
      b = make_cvector_from_rarray(argv[0]);
      flagb = 1;
    } else {
      CHECK_VECTOR(argv[0]);
      Data_Get_Struct(argv[0], gsl_vector, b);
    }
    S = gsl_vector_alloc(A->size2);   /* see manual p 123 */
    V = gsl_matrix_alloc(A->size2, A->size2);
    gsl_linalg_SV_decomp_jacobi(U, V, S);
    flagv = 1;
    break;
  }
  //  x = gsl_vector_alloc(b->size);  
  // Bug report #25842
  x = gsl_vector_alloc(S->size);
  gsl_linalg_SV_solve(U, V, S, b, x);
  if (flagv == 1) {
    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_vector_free(S);
  }
  if (flagb == 1) gsl_vector_free(b);
  return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
}

/*****/

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_cholesky_decomp_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na;
  VALUE chol;
  gsl_matrix_view mv;
  GetNArray(argv[0], na);
  chol = na_make_object(NA_DFLOAT, 2, na->shape, CLASS_OF(argv[0]));
  memcpy(NA_PTR_TYPE(chol,double*), (double*)na->ptr, sizeof(double)*na->total);
  mv = gsl_matrix_view_array(NA_PTR_TYPE(chol,double*), na->shape[1], na->shape[0]);
  gsl_linalg_cholesky_decomp(&mv.matrix);
  return chol;
}

static VALUE rb_gsl_linalg_cholesky_solve_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *nm, *nb;
  VALUE x;
  gsl_matrix_view mv;
  gsl_vector_view bv, xv;
  switch (argc) {
  case 2:
    GetNArray(argv[0], nm);
    GetNArray(argv[1], nb);
    x = na_make_object(NA_DFLOAT, 1, nb->shape, CLASS_OF(argv[1]));
    break;
  case 3:
    GetNArray(argv[0], nm);
    GetNArray(argv[1], nb);
    x = argv[2];
    break;
  default:
    rb_raise(rb_eArgError, 
	     "Usage: Cholesky.solve(chol, b) or Cholesky.solve(chol, b, x)");
    break;
  }
  mv = gsl_matrix_view_array((double*)nm->ptr, nm->shape[1], nm->shape[0]);
  bv = gsl_vector_view_array((double*)nb->ptr, nb->shape[0]);
  xv = gsl_vector_view_array(NA_PTR_TYPE(x,double*), nb->shape[0]);
  gsl_linalg_cholesky_solve(&mv.matrix, &bv.vector, &xv.vector);
  return x;
}

static VALUE rb_gsl_linalg_cholesky_svx_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *nm, *nb;
  gsl_matrix_view mv;
  gsl_vector_view bv;
  GetNArray(argv[0], nm);  GetNArray(argv[1], nb);
  mv = gsl_matrix_view_array((double*)nm->ptr, nm->shape[1], nm->shape[0]);
  bv = gsl_vector_view_array((double*)nb->ptr, nb->shape[0]);
  gsl_linalg_cholesky_svx(&mv.matrix, &bv.vector);
  return argv[1];
}

#endif

static VALUE rb_gsl_linalg_cholesky_decomp(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *Atmp = NULL;
  switch(TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0]))
      return rb_gsl_linalg_cholesky_decomp_narray(argc, argv, obj);
#endif
    CHECK_MATRIX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix, Atmp);
    break;
  default:
    CHECK_MATRIX(obj);
    Data_Get_Struct(obj, gsl_matrix, Atmp);
    break;
  }
  A = make_matrix_clone(Atmp);
  gsl_linalg_cholesky_decomp(A);
  return Data_Wrap_Struct(cgsl_matrix_C, 0, gsl_matrix_free, A);
}

static VALUE rb_gsl_linalg_cholesky_solve(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *Atmp = NULL;
  gsl_vector *b = NULL, *x = NULL;
  int flagb = 0, flaga = 0;
  VALUE vA, vb;
  switch(TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of argument (%d for 2)",
			    argc);
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0]))
      return rb_gsl_linalg_cholesky_solve_narray(argc, argv, obj);
#endif
    vA = argv[0];
    vb = argv[1];
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    vA = obj;
    vb = argv[0];
    break;
  }
  CHECK_MATRIX(vA);
  Data_Get_Struct(vA, gsl_matrix, Atmp);
  if (TYPE(vb) == T_ARRAY) {
    b = make_cvector_from_rarray(vb);
    flagb = 1;
  } else {
    CHECK_VECTOR(vb);
    Data_Get_Struct(vb, gsl_vector, b);
  }
  if (CLASS_OF(vA) == cgsl_matrix_C) {
    A = Atmp;
  } else {
    A = make_matrix_clone(Atmp);
    flaga = 1;
    gsl_linalg_cholesky_decomp(A);
  }
  x = gsl_vector_alloc(b->size);
  gsl_linalg_cholesky_solve(A, b, x);
  if (flaga == 1) gsl_matrix_free(A);
  if (flagb == 1) gsl_vector_free(b);
  return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
}


static VALUE rb_gsl_linalg_cholesky_svx(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *Atmp = NULL;
  gsl_vector *b = NULL;
  int flaga = 0;
  VALUE vA, vb;
  switch(TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of argument (%d for 2)",
			    argc);
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0]))
      return rb_gsl_linalg_cholesky_svx_narray(argc, argv, obj);
#endif
    vA = argv[0];
    vb = argv[1];
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    vA = obj;
    vb = argv[0];
    break;
  }
  CHECK_MATRIX(vA);
  Data_Get_Struct(vA, gsl_matrix, Atmp);
  CHECK_VECTOR(vb);
  Data_Get_Struct(vb, gsl_vector, b);
  if (CLASS_OF(vA) == cgsl_matrix_C) {
    A = Atmp;
  } else {
    A = make_matrix_clone(Atmp);
    flaga = 1;
    gsl_linalg_cholesky_decomp(A);
  }
  gsl_linalg_cholesky_svx(A, b);
  if (flaga == 1) gsl_matrix_free(A);
  return vb;
}

static VALUE rb_gsl_linalg_symmtd_decomp(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *Atmp = NULL;
  gsl_vector *tau = NULL;
  VALUE vQ, vtau;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    CHECK_MATRIX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix, Atmp);
    break;
  default:
    CHECK_MATRIX(obj);
    Data_Get_Struct(obj, gsl_matrix, Atmp);
    break;
  }
  A = make_matrix_clone(Atmp);
  tau = gsl_vector_alloc(A->size1);
  gsl_linalg_symmtd_decomp(A, tau);
  vQ = Data_Wrap_Struct(cgsl_matrix_Q, 0, gsl_matrix_free, A);
  vtau = Data_Wrap_Struct(cgsl_vector_tau, 0, gsl_vector_free, tau);
  return rb_ary_new3(2, vQ, vtau);
}


static VALUE rb_gsl_linalg_symmtd_decomp2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *tau = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    CHECK_MATRIX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    break;
  default:
    CHECK_MATRIX(obj);
    Data_Get_Struct(obj, gsl_matrix, A);
    break;
  }
  tau = gsl_vector_alloc(A->size1);
  gsl_linalg_symmtd_decomp(A, tau);
  return Data_Wrap_Struct(cgsl_vector_tau, 0, gsl_vector_free, tau);
}

static VALUE rb_gsl_linalg_symmtd_unpack(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *Q = NULL;
  gsl_vector *tau = NULL, *d = NULL, *sd = NULL;
  VALUE vq, vd, vsd;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of argument (%d for 2)",
			    argc);
    CHECK_MATRIX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    Data_Get_Struct(argv[1], gsl_vector, tau);
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    CHECK_MATRIX(obj);
    Data_Get_Struct(obj, gsl_matrix, A);
    Data_Get_Struct(argv[0], gsl_vector, tau);
    break;
  }
  Q = gsl_matrix_alloc(A->size1, A->size2);
  d = gsl_vector_alloc(tau->size);
  sd = gsl_vector_alloc(tau->size);
  gsl_linalg_symmtd_unpack(A, tau, Q, d, sd);

  vq = Data_Wrap_Struct(cgsl_matrix_Q, 0, gsl_matrix_free, Q);
  vd = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, d);
  vsd = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, sd);

  return rb_ary_new3(3, vq, vd, vsd);
}

static VALUE rb_gsl_linalg_symmtd_unpack_T(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *d = NULL, *sd = NULL;
  VALUE vd, vsd;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 2)",
			    argc);
    CHECK_MATRIX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix, A);
    break;
  }
  d = gsl_vector_alloc(A->size1);
  sd = gsl_vector_alloc(A->size1);
  gsl_linalg_symmtd_unpack_T(A, d, sd);

  vd = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, d);
  vsd = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, sd);

  return rb_ary_new3(2, vd, vsd);
}

/*****/

static VALUE rb_gsl_linalg_hermtd_decomp(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL, *Atmp = NULL;
  gsl_vector_complex *tau = NULL;
  VALUE vQ, vtau;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    CHECK_MATRIX_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, Atmp);
    break;
  default:
    CHECK_MATRIX_COMPLEX(obj);
    Data_Get_Struct(obj, gsl_matrix_complex, Atmp);
    break;
  }
  A = make_matrix_complex_clone(Atmp);
  tau = gsl_vector_complex_alloc(A->size1);
  gsl_linalg_hermtd_decomp(A, tau);
  vQ = Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, A);
  vtau = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, tau);
  return rb_ary_new3(2, vQ, vtau);
}

static VALUE rb_gsl_linalg_hermtd_decomp2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector_complex *tau = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    CHECK_MATRIX_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, A);
    break;
  default:
    CHECK_MATRIX_COMPLEX(obj);
    Data_Get_Struct(obj, gsl_matrix_complex, A);
    break;
  }
  tau = gsl_vector_complex_alloc(A->size1);
  gsl_linalg_hermtd_decomp(A, tau);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, tau);
}

static VALUE rb_gsl_linalg_hermtd_unpack(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL, *Q = NULL;
  gsl_vector_complex *tau = NULL;
  gsl_vector *d = NULL, *sd = NULL;
  VALUE vq, vd, vsd;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    CHECK_MATRIX_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, A);
    Data_Get_Struct(argv[1], gsl_vector_complex, tau);
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    CHECK_MATRIX_COMPLEX(obj);
    Data_Get_Struct(obj, gsl_matrix_complex, A);
    Data_Get_Struct(argv[0], gsl_vector_complex, tau);
    break;
  }
  Q = gsl_matrix_complex_alloc(A->size1, A->size2);
  d = gsl_vector_alloc(tau->size);
  sd = gsl_vector_alloc(tau->size);
  gsl_linalg_hermtd_unpack(A, tau, Q, d, sd);

  vq = Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, Q);
  vd = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, d);
  vsd = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, sd);

  return rb_ary_new3(3, vq, vd, vsd);
}

static VALUE rb_gsl_linalg_hermtd_unpack_T(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL;
  gsl_vector *d = NULL, *sd = NULL;
  VALUE vd, vsd;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    CHECK_MATRIX_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, A);
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix_complex, A);
    break;
  }
  d = gsl_vector_alloc(A->size1);
  sd = gsl_vector_alloc(A->size1);
  gsl_linalg_hermtd_unpack_T(A, d, sd);

  vd = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, d);
  vsd = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, sd);

  return rb_ary_new3(2, vd, vsd);
}

/******/

static VALUE rb_gsl_linalg_bidiag_decomp(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *Atmp = NULL;
  gsl_vector *tau_U = NULL, *tau_V = NULL;
  size_t size0;
  // local variable "status" was defined and set, but never used
  //int status;
  VALUE vu, vv, vA;

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
    Data_Get_Struct(argv[0], gsl_matrix, Atmp);
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix, Atmp);
    break;
  }
  A = make_matrix_clone(Atmp);
  size0 = GSL_MIN(A->size1, A->size2);
  tau_U = gsl_vector_alloc(size0);
  tau_V = gsl_vector_alloc(size0-1);
  /*status =*/ gsl_linalg_bidiag_decomp(A, tau_U, tau_V);
  vA = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, A);
  vu = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, tau_U);
  vv = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, tau_V);
  return rb_ary_new3(3, vA, vu, vv);
}

static VALUE rb_gsl_linalg_bidiag_decomp2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *tau_U = NULL, *tau_V = NULL;
  size_t size0;
  VALUE vu, vv;

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix, A);
    break;
  }
  size0 = GSL_MIN(A->size1, A->size2);
  tau_U = gsl_vector_alloc(size0);
  tau_V = gsl_vector_alloc(size0-1);
  gsl_linalg_bidiag_decomp(A, tau_U, tau_V);
  vu = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, tau_U);
  vv = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, tau_V);
  return rb_ary_new3(2, vu, vv);
}

static VALUE rb_gsl_linalg_bidiag_unpack(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *U = NULL, *V = NULL;
  gsl_vector *tau_U = NULL, *tau_V = NULL, *d = NULL, *s = NULL;
  size_t size0;
  VALUE vu, vv, vd, vs;

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 3) rb_raise(rb_eArgError, "wrong number of arguments (%d for 3)",
			    argc);
    CHECK_MATRIX(argv[0]);
    CHECK_VECTOR(argv[1]);
    CHECK_VECTOR(argv[2]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    Data_Get_Struct(argv[1], gsl_vector, tau_U);
    Data_Get_Struct(argv[2], gsl_vector, tau_V);
    break;
  default:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
			    argc);
    CHECK_MATRIX(obj);
    CHECK_VECTOR(argv[0]);
    CHECK_VECTOR(argv[1]);
    Data_Get_Struct(obj, gsl_matrix, A);
    Data_Get_Struct(argv[0], gsl_vector, tau_U);
    Data_Get_Struct(argv[1], gsl_vector, tau_V);
    break;
  }
  size0 = GSL_MIN(A->size1, A->size2);
  U = gsl_matrix_alloc(A->size1, A->size2);
  V = gsl_matrix_alloc(size0, size0);

  d = gsl_vector_alloc(size0);
  s = gsl_vector_alloc(size0-1);
  gsl_linalg_bidiag_unpack(A, tau_U, U, tau_V, V, d, s);
  vu = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, U);
  vv = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, V);
  vd = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, d);
  vs = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, s);
  return rb_ary_new3(4, vu, vv, vd, vs);
}

static VALUE rb_gsl_linalg_bidiag_unpack2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *V = NULL;
  gsl_vector *tau_V = NULL, *tau_U = NULL;
  VALUE vv;

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
       if (argc != 3) rb_raise(rb_eArgError, "wrong number of arguments (%d for 3)",
			    argc);
    CHECK_MATRIX(argv[0]);
    CHECK_VECTOR(argv[1]);
    CHECK_VECTOR(argv[2]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    Data_Get_Struct(argv[1], gsl_vector, tau_U);
    Data_Get_Struct(argv[2], gsl_vector, tau_V);
    break;
  default:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
			    argc);
    CHECK_MATRIX(obj);
    CHECK_VECTOR(argv[0]);
    CHECK_VECTOR(argv[1]);
    Data_Get_Struct(obj, gsl_matrix, A);
    Data_Get_Struct(argv[0], gsl_vector, tau_U);
    Data_Get_Struct(argv[1], gsl_vector, tau_V);
    break;
  } 
  V = gsl_matrix_alloc(A->size2, A->size2);
  gsl_linalg_bidiag_unpack2(A, tau_U, tau_V, V);
  vv = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, V);
  return vv;
}

static VALUE rb_gsl_linalg_bidiag_unpack_B(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *d = NULL, *s = NULL;
  size_t size0;
  VALUE vd, vs;

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 3)",
			    argc);
    CHECK_MATRIX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    break;
  default:
    CHECK_MATRIX(obj);
    Data_Get_Struct(obj, gsl_matrix, A);
    break;
  } 
  size0 = GSL_MIN(A->size1, A->size2);
  d = gsl_vector_alloc(size0);
  s = gsl_vector_alloc(size0);
  gsl_linalg_bidiag_unpack_B(A, d, s);
  vd = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, d);
  vs = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, s);
  return rb_ary_new3(2, vd, vs);
}

/* Householder Transformations 11.Jul.2004 */
static VALUE rb_gsl_linalg_householder_transform(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments.");
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, v);
    break;
  default:
    Data_Get_Struct(obj, gsl_vector, v);
    break;
  }
  return rb_float_new(gsl_linalg_householder_transform(v));
}

/* singleton */
static VALUE rb_gsl_linalg_householder_hm(VALUE obj, VALUE t, VALUE vv, VALUE aa)
{
  gsl_vector *v = NULL;
  double tau;
  gsl_matrix *A = NULL;
  CHECK_VECTOR(vv);
  CHECK_MATRIX(aa);
  tau = NUM2DBL(t);
  Data_Get_Struct(vv, gsl_vector, v);
  Data_Get_Struct(aa, gsl_matrix, A);
  gsl_linalg_householder_hm(tau, v, A);
  return aa;
}

static VALUE rb_gsl_linalg_householder_mh(VALUE obj, VALUE t, VALUE vv, VALUE aa)
{
  gsl_vector *v = NULL;
  double tau;
  gsl_matrix *A = NULL;
  CHECK_VECTOR(vv);
  CHECK_MATRIX(aa);
  tau = NUM2DBL(t);
  Data_Get_Struct(vv, gsl_vector, v);
  Data_Get_Struct(aa, gsl_matrix, A);
  gsl_linalg_householder_mh(tau, v, A);
  return aa;
}

static VALUE rb_gsl_linalg_householder_hv(VALUE obj, VALUE t, VALUE vv, VALUE ww)
{
  gsl_vector *v = NULL, *w = NULL;
  double tau;
  CHECK_VECTOR(vv);
  CHECK_VECTOR(ww);
  tau = NUM2DBL(t);
  Data_Get_Struct(vv, gsl_vector, v);
  Data_Get_Struct(ww, gsl_vector, w);
  gsl_linalg_householder_hv(tau, v, w);
  return ww;
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_HH_solve_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na;
  gsl_vector_view bv, xv;
  VALUE x;
  gsl_matrix *mtmp;
  GetNArray(argv[0], na);
  bv = gsl_vector_view_array(NA_PTR_TYPE(argv[1],double*), na->shape[1]);
  x = na_make_object(NA_DFLOAT, 1, &na->shape[1], CLASS_OF(argv[1]));
  xv = gsl_vector_view_array(NA_PTR_TYPE(x,double*), na->shape[1]);
  mtmp = gsl_matrix_alloc(na->shape[1], na->shape[0]);
  memcpy(mtmp->data, (double*)na->ptr, sizeof(double)*na->total);
  gsl_linalg_HH_solve(mtmp, &bv.vector, &xv.vector);
  gsl_matrix_free(mtmp);
  return x;
}
static VALUE rb_gsl_linalg_HH_svx_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na;
  gsl_matrix *mtmp;
  gsl_vector_view bv;
  GetNArray(argv[0], na);
  bv = gsl_vector_view_array(NA_PTR_TYPE(argv[1],double*), na->shape[1]);
  mtmp = gsl_matrix_alloc(na->shape[1], na->shape[0]);
  memcpy(mtmp->data, (double*)na->ptr, sizeof(double)*na->total);
  gsl_linalg_HH_svx(mtmp, &bv.vector);
  gsl_matrix_free(mtmp);
  return argv[1];
}
#endif

/* 17.Apr.2004 */
static VALUE rb_gsl_linalg_HH_solve(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *Atmp = NULL;
  gsl_vector *b = NULL, *x = NULL;
  int flagb = 0;
  VALUE vA, vb;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of argument (%d for 2)",
			    argc);
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0]))
      return rb_gsl_linalg_HH_solve_narray(argc, argv, obj);
#endif
    vA = argv[0];
    vb = argv[1];
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    vA = obj;
    vb = argv[0];
    break;
  }
  CHECK_MATRIX(vA);
  Data_Get_Struct(vA, gsl_matrix, Atmp);
  if (TYPE(vb) == T_ARRAY) {
    b = make_cvector_from_rarray(vb);
    flagb = 1;
  } else {
    CHECK_VECTOR(vb);
    Data_Get_Struct(vb, gsl_vector, b);
  }
  A = make_matrix_clone(Atmp);
  x = gsl_vector_alloc(b->size);
  gsl_linalg_HH_solve(A, b, x);
  gsl_matrix_free(A);
  if (flagb == 1) gsl_vector_free(b);
  return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
}

static VALUE rb_gsl_linalg_HH_solve_bang(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *b = NULL, *x = NULL;
  int flagb = 0;
  VALUE vA, vb;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of argument (%d for 2)",
			    argc);
    vA = argv[0];
    vb = argv[1];
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    vA = obj;
    vb = argv[0];
    break;
  }
  CHECK_MATRIX(vA);
  Data_Get_Struct(vA, gsl_matrix, A);
  if (TYPE(vb) == T_ARRAY) {
    b = make_cvector_from_rarray(vb);
    flagb = 1;
  } else {
    CHECK_VECTOR(vb);
    Data_Get_Struct(vb, gsl_vector, b);
  }
  x = gsl_vector_alloc(b->size);
  gsl_linalg_HH_solve(A, b, x);
  if (flagb == 1) gsl_vector_free(b);
  return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
}

static VALUE rb_gsl_linalg_HH_svx(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *Atmp = NULL;
  gsl_vector *b = NULL;
  VALUE vA, vb;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of argument (%d for 2)",
			    argc);
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0]))
      return rb_gsl_linalg_HH_svx_narray(argc, argv, obj);
#endif
    vA = argv[0];
    vb = argv[1];
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of argument (%d for 1)",
			    argc);
    vA = obj;
    vb = argv[0];
    break;
  }
  CHECK_MATRIX(vA);
  Data_Get_Struct(vA, gsl_matrix, Atmp);
  CHECK_VECTOR(vb);
  Data_Get_Struct(vb, gsl_vector, b);
  A = make_matrix_clone(Atmp);
  gsl_linalg_HH_svx(A, b);
  gsl_matrix_free(A);
  return vb;
}

static VALUE rb_gsl_linalg_solve_symm_tridiag(VALUE obj, VALUE dd, VALUE ee, VALUE bb)
{
  gsl_vector *b = NULL, *x = NULL, *d = NULL, *e = NULL;

  Data_Get_Struct(dd, gsl_vector, d);
  Data_Get_Struct(ee, gsl_vector, e);
  Data_Get_Struct(bb, gsl_vector, b);
  x = gsl_vector_alloc(b->size);

  gsl_linalg_solve_symm_tridiag(d, e, b, x);
  return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
}

#ifdef GSL_1_2_LATER
static VALUE rb_gsl_linalg_solve_tridiag(VALUE obj, VALUE dd, VALUE ee, VALUE ff,
					 VALUE bb)
{
  gsl_vector *b = NULL, *x = NULL, *d = NULL, *e = NULL, *f = NULL;

  Data_Get_Struct(dd, gsl_vector, d);
  Data_Get_Struct(ee, gsl_vector, e);
  Data_Get_Struct(ff, gsl_vector, f);
  Data_Get_Struct(bb, gsl_vector, b);
  x = gsl_vector_alloc(b->size);

  gsl_linalg_solve_tridiag(d, e, f, b, x);
  return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
}

static VALUE rb_gsl_linalg_solve_symm_cyc_tridiag(VALUE obj, VALUE dd, VALUE ee, VALUE bb)
{
  gsl_vector *b = NULL, *x = NULL, *d = NULL, *e = NULL;

  Data_Get_Struct(dd, gsl_vector, d);
  Data_Get_Struct(ee, gsl_vector, e);
  Data_Get_Struct(bb, gsl_vector, b);
  x = gsl_vector_alloc(b->size);

  gsl_linalg_solve_symm_cyc_tridiag(d, e, b, x);
  return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
}

static VALUE rb_gsl_linalg_solve_cyc_tridiag(VALUE obj, VALUE dd, VALUE ee, 
					     VALUE ff, VALUE bb)
{
  gsl_vector *b = NULL, *x = NULL, *d = NULL, *e = NULL, *f = NULL;
  Data_Get_Struct(dd, gsl_vector, d);
  Data_Get_Struct(ee, gsl_vector, e);
  Data_Get_Struct(ff, gsl_vector, f);
  Data_Get_Struct(bb, gsl_vector, b);
  x = gsl_vector_alloc(b->size);
  gsl_linalg_solve_cyc_tridiag(d, e, f, b, x);
  return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, x);
}
#endif

static void rb_gsl_linalg_balance_columns_init(int argc, VALUE *argv, VALUE obj,
						VALUE *mat, VALUE *vec,
						gsl_matrix **M, gsl_vector **V)
{
  gsl_matrix *A = NULL;
  gsl_vector *D = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    switch (argc) {
    case 2:
      CHECK_MATRIX(argv[0]); CHECK_VECTOR(argv[1]);
      Data_Get_Struct(argv[0], gsl_matrix, A);
      Data_Get_Struct(argv[1], gsl_vector, D);
      *vec = argv[1];
      break;
    case 1:
      CHECK_MATRIX(argv[0]);
      Data_Get_Struct(argv[0], gsl_matrix, A);
      D = gsl_vector_alloc(A->size2);
      *vec = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, D);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
      break;
    }
    *mat = argv[0];
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix, A);
    switch (argc) {
    case 1:
      CHECK_VECTOR(argv[0]);
      Data_Get_Struct(argv[0], gsl_vector, D);
      *vec = argv[0];
      break;
    case 0:
      D = gsl_vector_alloc(A->size2);
      *vec = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, D);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
      break;
    }
    *mat = obj;
    break;
  }
  *M = A;
  *V = D;
}

static VALUE rb_gsl_linalg_balance_columns_bang(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL;
  gsl_vector *D = NULL;
  VALUE mat, vec;
  // local variable "status" was defined and set, but never used
  //int status;
  rb_gsl_linalg_balance_columns_init(argc, argv, obj, &mat, &vec, &A, &D);
  /*status =*/ gsl_linalg_balance_columns(A, D);
  return rb_ary_new3(2, mat, vec);
}

static VALUE rb_gsl_linalg_balance_columns(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *A = NULL, *Anew;
  gsl_vector *D = NULL;
  VALUE mat, vec;
  // local variable "status" was defined and set, but never used
  //int status;
  rb_gsl_linalg_balance_columns_init(argc, argv, obj, &mat, &vec, &A, &D);
  Anew = make_matrix_clone(A);
  mat = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Anew);
  /*status =*/ gsl_linalg_balance_columns(Anew, D);
  return rb_ary_new3(2, mat, vec);
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_linalg_QR_decomp_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na;
  gsl_matrix_view mv;
  gsl_vector_view vv;
  int shapem[2], shapev[1];
  VALUE qr, tau;
  if (argc < 1) rb_raise(rb_eArgError, "too few arguments.");
  GetNArray(argv[0], na);
  shapem[0] = na->shape[1];
  shapem[1] = na->shape[1];
  shapev[0] = shapem[0];
  qr = na_make_object(NA_DFLOAT, 2, shapem, CLASS_OF(argv[0]));
  tau = na_make_object(NA_DFLOAT, 1, shapev, cNVector);
  memcpy(NA_PTR_TYPE(qr,double*),na->ptr,sizeof(double)*shapem[0]*shapem[1]);
  mv = gsl_matrix_view_array(NA_PTR_TYPE(qr,double*), shapem[0], shapem[1]);
  vv = gsl_vector_view_array(NA_PTR_TYPE(tau,double*), shapev[0]);
  gsl_linalg_QR_decomp(&mv.matrix, &vv.vector);
  return rb_ary_new3(2, qr, tau);
}

static VALUE rb_gsl_linalg_QR_unpack_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *m, *tau;
  gsl_matrix_view mv, mq, mr;
  gsl_vector_view vv;
  int shape[2];
  VALUE q, r;
  if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
			  argc);
  GetNArray(argv[0], m);
  GetNArray(argv[1], tau);
  mv = gsl_matrix_view_array((double*)m->ptr, m->shape[1], m->shape[0]);
  vv = gsl_vector_view_array((double*)tau->ptr, tau->shape[0]);
  shape[0] = m->shape[1];
  shape[1] = m->shape[1];
  q = na_make_object(NA_DFLOAT, 2, shape, CLASS_OF(argv[0]));
  shape[0] = m->shape[1];
  shape[1] = m->shape[0];
  r = na_make_object(NA_DFLOAT, 2, shape, CLASS_OF(argv[0]));
  mq = gsl_matrix_view_array(NA_PTR_TYPE(q,double*), m->shape[1], m->shape[1]);
  mr = gsl_matrix_view_array(NA_PTR_TYPE(r,double*), m->shape[1], m->shape[0]);
  //  printf("OK 4 %d %d\n", mq.matrix.size1, mr.matrix.size2);
  gsl_linalg_QR_unpack(&mv.matrix, &vv.vector, &mq.matrix, &mr.matrix);
  //  printf("OK 5\n");
  return rb_ary_new3(2, q, r);
}

static VALUE rb_gsl_linalg_QR_solve_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *qr, *tau, *b;
  VALUE x;
  gsl_matrix_view mv;
  gsl_vector_view tv, bv, xv;
  if (argc != 3) rb_raise(rb_eArgError, "Usage: QR.solve(qr, tau, b)");
  GetNArray(argv[0], qr);
  GetNArray(argv[1], tau);
  GetNArray(argv[2], b);
  x = na_make_object(NA_DFLOAT, 1, b->shape, CLASS_OF(argv[2]));
  mv = gsl_matrix_view_array((double*)qr->ptr, qr->shape[1], qr->shape[0]);
  tv = gsl_vector_view_array((double*)tau->ptr, tau->shape[0]);
  bv = gsl_vector_view_array((double*)b->ptr, b->shape[0]);
  xv = gsl_vector_view_array(NA_PTR_TYPE(x,double*), b->shape[0]);
  gsl_linalg_QR_solve(&mv.matrix, &tv.vector, &bv.vector, &xv.vector);
  return x;
}
static VALUE rb_gsl_linalg_QR_svx_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *qr, *tau, *b;
  gsl_matrix_view mv;
  gsl_vector_view tv, bv;
  if (argc != 3) rb_raise(rb_eArgError, "Usage: QR.solve(qr, tau, b)");
  GetNArray(argv[0], qr);
  GetNArray(argv[1], tau);
  GetNArray(argv[2], b);
  mv = gsl_matrix_view_array((double*)qr->ptr, qr->shape[1], qr->shape[0]);
  tv = gsl_vector_view_array((double*)tau->ptr, tau->shape[0]);
  bv = gsl_vector_view_array((double*)b->ptr, b->shape[0]);
  gsl_linalg_QR_svx(&mv.matrix, &tv.vector, &bv.vector);
  return argv[2];
}

#endif

#ifdef GSL_1_9_LATER
static VALUE rb_gsl_linalg_hessenberg_decomp(VALUE module, VALUE AA)
{
  gsl_matrix *A = NULL, *Atmp = NULL;
  gsl_vector *tau = NULL;
  VALUE vH, vtau;
  CHECK_MATRIX(AA);
  Data_Get_Struct(AA, gsl_matrix, Atmp);
  A = make_matrix_clone(Atmp);
  tau = gsl_vector_alloc(A->size1);
  gsl_linalg_hessenberg_decomp(A, tau);
  vH = Data_Wrap_Struct(cgsl_matrix_Q, 0, gsl_matrix_free, A);
  vtau = Data_Wrap_Struct(cgsl_vector_tau, 0, gsl_vector_free, tau);
  return rb_ary_new3(2, vH, vtau);
}

static VALUE rb_gsl_linalg_hessenberg_unpack(VALUE module, VALUE HH, VALUE tt)
{
  gsl_matrix *H = NULL, *U = NULL;
  gsl_vector *tau = NULL;
  CHECK_MATRIX(HH);
  CHECK_VECTOR(tt);
  Data_Get_Struct(HH, gsl_matrix, H);
  Data_Get_Struct(tt, gsl_vector, tau);
  U = gsl_matrix_alloc(H->size1, H->size2);
  gsl_linalg_hessenberg_unpack(H, tau, U);

  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, U);
}

static VALUE rb_gsl_linalg_hessenberg_unpack_accum(int argc, VALUE *argv, VALUE module)
{
  gsl_matrix *H = NULL, *V = NULL;
  gsl_vector *tau = NULL;
  size_t i;
  VALUE val = Qnil;
  switch (argc) {
  case 2:
    /* nothing to do */
    break;
  case 3:
    CHECK_MATRIX(argv[2]);
    Data_Get_Struct(argv[2], gsl_matrix, V);
    val = argv[2];
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 2 or 3)", argc);
  }
  CHECK_MATRIX(argv[0]);
  CHECK_VECTOR(argv[1]);
  Data_Get_Struct(argv[0], gsl_matrix, H);
  Data_Get_Struct(argv[1], gsl_vector, tau);
  if (argc == 2) {
    V = gsl_matrix_alloc(H->size1, H->size2);
    val = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, V);
    for (i = 0; i < V->size1; i++) gsl_matrix_set(V, i, i, 1.0);
  }
  gsl_linalg_hessenberg_unpack_accum(H, tau, V);
  return val;
}
static VALUE rb_gsl_linalg_hessenberg_set_zero(VALUE module, VALUE HH)
{
  gsl_matrix *H;
  CHECK_MATRIX(HH);
  Data_Get_Struct(HH, gsl_matrix, H);
  return INT2FIX(gsl_linalg_hessenberg_set_zero(H));
  /*  gsl_linalg_hessenberg_set_zero(H);
      return INT2FIX(0);*/
}
static VALUE rb_gsl_linalg_hesstri_decomp(int argc, VALUE *argv, VALUE module)
{
  gsl_matrix *A = NULL, *B = NULL, *Anew, *Bnew;
  gsl_matrix *U = NULL, *V = NULL;
  gsl_vector *work = NULL;
  VALUE vH, vR, vU = Qnil, vV = Qnil, ary;
  int flag = 0;
  switch (argc) {
  case 2:
    flag = 1;
    break;
  case 3:
    CHECK_VECTOR(argv[2]);
    Data_Get_Struct(argv[2], gsl_vector, work); 
    break;
  case 4:
    CHECK_MATRIX(argv[2]);
    CHECK_MATRIX(argv[3]);
    Data_Get_Struct(argv[2], gsl_matrix, U);
    Data_Get_Struct(argv[3], gsl_matrix, V);
    flag = 1;
    break;  
  case 5:
    CHECK_MATRIX(argv[2]);
    CHECK_MATRIX(argv[3]);
    CHECK_VECTOR(argv[4]);
    Data_Get_Struct(argv[2], gsl_matrix, U);
    Data_Get_Struct(argv[3], gsl_matrix, V);
    Data_Get_Struct(argv[4], gsl_vector, work);
    vU = argv[2];
    vV = argv[3];
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 2-55)", argc);
  }
  CHECK_MATRIX(argv[0]);
  CHECK_MATRIX(argv[1]);      
  Data_Get_Struct(argv[0], gsl_matrix, A);
  Data_Get_Struct(argv[1], gsl_matrix, B);
  Anew = make_matrix_clone(A);
  Bnew = make_matrix_clone(B);  
  if (flag == 1) work = gsl_vector_alloc(A->size1);
  gsl_linalg_hesstri_decomp(Anew, Bnew, U, V, work);
  if (flag == 1) gsl_vector_free(work);
  vH = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Anew);
  vR = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Bnew);  
  if (argc == 2 || argc == 3) {    
    ary = rb_ary_new3(2, vH, vR);
  } else {
    ary = rb_ary_new3(4, vH, vR, vU, vV);
  }
  return ary;
}
static VALUE rb_gsl_linalg_hesstri_decomp_bang(int argc, VALUE *argv, VALUE module)
{
  gsl_matrix *A = NULL, *B = NULL;
  gsl_matrix *U = NULL, *V = NULL;
  gsl_vector *work = NULL;
  VALUE vH, vR, vU = Qnil, vV = Qnil, ary;
  int flag = 0;
  switch (argc) {
  case 2:
    flag = 1;
    break;
  case 3:
    CHECK_VECTOR(argv[2]);
    Data_Get_Struct(argv[2], gsl_vector, work); 
    break;
  case 4:
    CHECK_MATRIX(argv[2]);
    CHECK_MATRIX(argv[3]);
    Data_Get_Struct(argv[2], gsl_matrix, U);
    Data_Get_Struct(argv[3], gsl_matrix, V);
    flag = 1;
    break;  
  case 5:
    CHECK_MATRIX(argv[2]);
    CHECK_MATRIX(argv[3]);
    CHECK_VECTOR(argv[4]);
    Data_Get_Struct(argv[2], gsl_matrix, U);
    Data_Get_Struct(argv[3], gsl_matrix, V);
    Data_Get_Struct(argv[4], gsl_vector, work);
    vU = argv[2];
    vV = argv[3];
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 2-55)", argc);
  }
  CHECK_MATRIX(argv[0]);
  CHECK_MATRIX(argv[1]);      
  Data_Get_Struct(argv[0], gsl_matrix, A);
  Data_Get_Struct(argv[1], gsl_matrix, B);
  if (flag == 1) work = gsl_vector_alloc(A->size1);
  gsl_linalg_hesstri_decomp(A, B, U, V, work);
  if (flag == 1) gsl_vector_free(work);
  vH = argv[0];
  vR = argv[1];
  if (argc == 2 || argc == 3) {    
    ary = rb_ary_new3(2, vH, vR);
  } else {
    ary = rb_ary_new3(4, vH, vR, vU, vV);
  }
  return ary;
}

static VALUE rb_gsl_linalg_balance_matrix(int argc, VALUE *argv, VALUE module)
{
  gsl_matrix *A, *Anew;
  gsl_vector *D;
  VALUE vA, vD;
  switch (argc) {
  case 1:
    CHECK_MATRIX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    Anew = make_matrix_clone(A);
    D = gsl_vector_alloc(A->size1);
    vD = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, D);
    break;
  case 2:
    CHECK_MATRIX(argv[0]);
    CHECK_VECTOR(argv[1]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    Data_Get_Struct(argv[1], gsl_vector, D);
    Anew = make_matrix_clone(A);  
    vD = argv[1];
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 1 or 2)", argc);
  }
  gsl_linalg_balance_matrix(Anew, D);
  vA = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Anew);
  return rb_ary_new3(2, vA, vD);

}
static VALUE rb_gsl_linalg_balance_matrix2(int argc, VALUE *argv, VALUE module)
{
  gsl_matrix *A;
  gsl_vector *D;
  switch (argc) {
  case 1:
    CHECK_MATRIX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    D = gsl_vector_alloc(A->size1);
    gsl_linalg_balance_matrix(A, D);
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, D);
    break;
  case 2:
    CHECK_MATRIX(argv[0]);
    CHECK_VECTOR(argv[1]);
    Data_Get_Struct(argv[0], gsl_matrix, A);
    Data_Get_Struct(argv[1], gsl_vector, D);
    return INT2FIX(gsl_linalg_balance_matrix(A, D));
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 1 or 2)", argc);
  }
  return Qtrue;
}
#endif

void Init_gsl_linalg_complex(VALUE module);
void Init_gsl_linalg(VALUE module)
{
  VALUE mgsl_linalg;
  VALUE mgsl_linalg_LU;
  VALUE mgsl_linalg_QR;
  VALUE mgsl_linalg_QRPT;
  VALUE mgsl_linalg_LQ;
  VALUE mgsl_linalg_PTLQ;
  VALUE mgsl_linalg_SV;
  VALUE mgsl_linalg_cholesky;
  VALUE mgsl_linalg_symmtd;
  VALUE mgsl_linalg_hermtd;
  VALUE mgsl_linalg_bidiag;
  VALUE mgsl_linalg_tridiag;
  VALUE mgsl_linalg_HH;
  VALUE mgsl_linalg_Householder;
#ifdef GSL_1_9_LATER
  VALUE mhessen;
#endif

  mgsl_linalg = rb_define_module_under(module, "Linalg");
  mgsl_linalg_LU = rb_define_module_under(mgsl_linalg, "LU");
  cgsl_matrix_LU = rb_define_class_under(mgsl_linalg_LU, "LUMatrix", cgsl_matrix);
  mgsl_linalg_QR = rb_define_module_under(mgsl_linalg, "QR");
  mgsl_linalg_QRPT = rb_define_module_under(mgsl_linalg, "QRPT");
  cgsl_matrix_QR = rb_define_class_under(mgsl_linalg, "QRMatrix", cgsl_matrix);
  cgsl_matrix_QRPT = rb_define_class_under(mgsl_linalg, "QRPTMatrix", cgsl_matrix);
  cgsl_vector_tau = rb_define_class_under(mgsl_linalg, "TauVector", cgsl_vector);
  cgsl_matrix_Q = rb_define_class_under(mgsl_linalg, "QMatrix", cgsl_matrix);
  cgsl_matrix_R = rb_define_class_under(mgsl_linalg, "RMatrix", cgsl_matrix);

  mgsl_linalg_LQ = rb_define_module_under(mgsl_linalg, "LQ");
  mgsl_linalg_PTLQ = rb_define_module_under(mgsl_linalg, "PTLQ");
  cgsl_matrix_LQ = rb_define_class_under(mgsl_linalg, "LQMatrix", cgsl_matrix);
  cgsl_matrix_PTLQ = rb_define_class_under(mgsl_linalg, "PTLQMatrix", cgsl_matrix);
  cgsl_matrix_L = rb_define_class_under(mgsl_linalg, "LMatrix", cgsl_matrix);

  /*****/
  mgsl_linalg_SV = rb_define_module_under(mgsl_linalg, "SV");
  cgsl_matrix_SV = rb_define_class_under(mgsl_linalg_SV, "SVMatrix", cgsl_matrix);
  cgsl_matrix_U = rb_define_class_under(mgsl_linalg_SV, "UMatrix", cgsl_matrix);
  cgsl_matrix_V = rb_define_class_under(mgsl_linalg_SV, "VMatrix", cgsl_matrix);
  cgsl_vector_S = rb_define_class_under(mgsl_linalg_SV, "SingularValues", cgsl_vector);

  /*****/
  mgsl_linalg_cholesky = rb_define_module_under(mgsl_linalg, "Cholesky");
  cgsl_matrix_C = rb_define_class_under(mgsl_linalg_cholesky, "CholeskyMatrix", cgsl_matrix);
  mgsl_linalg_symmtd = rb_define_module_under(mgsl_linalg, "Symmtd");

  mgsl_linalg_hermtd = rb_define_module_under(mgsl_linalg, "Hermtd");
  mgsl_linalg_bidiag = rb_define_module_under(mgsl_linalg, "Bidiag");
  mgsl_linalg_tridiag = rb_define_module_under(mgsl_linalg, "Tridiag");

  mgsl_linalg_HH = rb_define_module_under(mgsl_linalg, "HH");
  mgsl_linalg_Householder = rb_define_module_under(mgsl_linalg, "Householder");

  /*****/
  rb_define_module_function(mgsl_linalg, "LU_decomp!", rb_gsl_linalg_LU_decomp_bang, -1);
  rb_define_module_function(mgsl_linalg_LU, "decomp!", rb_gsl_linalg_LU_decomp_bang, -1);
  rb_define_module_function(mgsl_linalg, "LU_decomp", rb_gsl_linalg_LU_decomp, -1);
  rb_define_module_function(mgsl_linalg_LU, "decomp", rb_gsl_linalg_LU_decomp, -1);
  rb_define_method(cgsl_matrix, "LU_decomp!", rb_gsl_linalg_LU_decomp_bang, -1);
  rb_define_method(cgsl_matrix, "LU_decomp", rb_gsl_linalg_LU_decomp, -1);

  rb_define_module_function(mgsl_linalg, "LU_solve", rb_gsl_linalg_LU_solve, -1);
  rb_define_module_function(mgsl_linalg_LU, "solve", rb_gsl_linalg_LU_solve, -1);
  rb_define_method(cgsl_matrix, "LU_solve", rb_gsl_linalg_LU_solve, -1);
  rb_define_method(cgsl_matrix_LU, "solve", rb_gsl_linalg_LU_solve, -1);

  rb_define_module_function(mgsl_linalg, "LU_svx", rb_gsl_linalg_LU_svx, -1);
  rb_define_module_function(mgsl_linalg_LU, "svx", rb_gsl_linalg_LU_svx, -1);
  rb_define_method(cgsl_matrix, "LU_svx", rb_gsl_linalg_LU_svx, -1);
  rb_define_method(cgsl_matrix_LU, "svx", rb_gsl_linalg_LU_svx, -1);

  rb_define_module_function(mgsl_linalg, "LU_invert", rb_gsl_linalg_LU_invert, -1);
  rb_define_module_function(mgsl_linalg_LU, "invert", rb_gsl_linalg_LU_invert, -1);
  rb_define_module_function(mgsl_linalg_LU, "inv", rb_gsl_linalg_LU_invert, -1);
  rb_define_module_function(mgsl_linalg_LU, "refine", rb_gsl_linalg_LU_refine, 5);

  rb_define_method(cgsl_matrix, "invert", rb_gsl_linalg_LU_invert, -1);
  rb_define_alias(cgsl_matrix, "LU_invert", "invert");
  rb_define_alias(cgsl_matrix, "inv", "invert");

  rb_define_module_function(mgsl_linalg, "LU_det", rb_gsl_linalg_LU_det, -1);
  rb_define_module_function(mgsl_linalg_LU, "det", rb_gsl_linalg_LU_det, -1);
  rb_define_method(cgsl_matrix, "LU_det", rb_gsl_linalg_LU_det, -1);
  rb_define_alias(cgsl_matrix, "det", "LU_det");

  rb_define_module_function(mgsl_linalg, "LU_lndet", rb_gsl_linalg_LU_lndet, -1);
  rb_define_module_function(mgsl_linalg_LU, "lndet", rb_gsl_linalg_LU_lndet, -1);
  rb_define_method(cgsl_matrix, "LU_lndet", rb_gsl_linalg_LU_lndet, -1);
  rb_define_alias(cgsl_matrix, "lndet", "LU_lndet");

  rb_define_module_function(mgsl_linalg, "LU_sgndet", rb_gsl_linalg_LU_sgndet, -1);
  rb_define_module_function(mgsl_linalg_LU, "sgndet", rb_gsl_linalg_LU_sgndet, -1);
  rb_define_method(cgsl_matrix, "LU_sgndet", rb_gsl_linalg_LU_sgndet, -1);
  rb_define_alias(cgsl_matrix, "sgndet", "LU_sgndet");

  /*****/
  rb_define_module_function(mgsl_linalg, "QR_decomp", rb_gsl_linalg_QR_decomp, -1);
  rb_define_module_function(mgsl_linalg_QR, "decomp", rb_gsl_linalg_QR_decomp, -1);
  rb_define_method(cgsl_matrix, "QR_decomp", rb_gsl_linalg_QR_decomp, -1);
  rb_define_module_function(mgsl_linalg, "QR_decomp!", rb_gsl_linalg_QR_decomp_bang, -1);
  rb_define_module_function(mgsl_linalg_QR, "decomp!", rb_gsl_linalg_QR_decomp_bang, -1);
  rb_define_method(cgsl_matrix, "QR_decomp!", rb_gsl_linalg_QR_decomp_bang, -1);

  rb_define_module_function(mgsl_linalg, "QR_solve", rb_gsl_linalg_QR_solve, -1);
  rb_define_module_function(mgsl_linalg_QR, "solve", rb_gsl_linalg_QR_solve, -1);
  rb_define_module_function(mgsl_linalg, "QR_svx", rb_gsl_linalg_QR_svx, -1);
  rb_define_module_function(mgsl_linalg_QR, "svx", rb_gsl_linalg_QR_svx, -1);
  rb_define_method(cgsl_matrix, "QR_solve", rb_gsl_linalg_QR_solve, -1);
  rb_define_method(cgsl_matrix_QR, "solve", rb_gsl_linalg_QR_solve, -1);
  rb_define_method(cgsl_matrix, "QR_svx", rb_gsl_linalg_QR_svx, -1);
  rb_define_method(cgsl_matrix_QR, "svx", rb_gsl_linalg_QR_svx, -1);

  rb_define_module_function(mgsl_linalg_QR, "lssolve", rb_gsl_linalg_QR_lssolve, -1);
  rb_define_method(cgsl_matrix, "QR_lssolve", rb_gsl_linalg_QR_lssolve, -1);
  rb_define_method(cgsl_matrix_QR, "lssolve", rb_gsl_linalg_QR_lssolve, -1);

  rb_define_module_function(mgsl_linalg_QR, "QTvec", rb_gsl_linalg_QR_QTvec, -1);
  rb_define_method(cgsl_matrix_QR, "QTvec", rb_gsl_linalg_QR_QTvec, -1);
 rb_define_module_function(mgsl_linalg_QR, "Qvec", rb_gsl_linalg_QR_Qvec, -1);
  rb_define_method(cgsl_matrix_QR, "Qvec", rb_gsl_linalg_QR_Qvec, -1);

  rb_define_module_function(mgsl_linalg_QR, "Rsolve", rb_gsl_linalg_QR_Rsolve, -1);
  rb_define_method(cgsl_matrix, "QR_Rsolve", rb_gsl_linalg_QR_Rsolve, -1);
  rb_define_method(cgsl_matrix_QR, "Rsolve", rb_gsl_linalg_QR_Rsolve, -1);

  rb_define_module_function(mgsl_linalg_QR, "Rsvx", rb_gsl_linalg_QR_Rsvx, -1);
  rb_define_method(cgsl_matrix_QR, "Rsvx", rb_gsl_linalg_QR_Rsvx, 1);

  rb_define_module_function(mgsl_linalg_QR, "unpack", rb_gsl_linalg_QR_unpack, -1);
  rb_define_method(cgsl_matrix_QR, "unpack", rb_gsl_linalg_QR_unpack, -1);

  rb_define_module_function(mgsl_linalg_QR, "QRsolve", rb_gsl_linalg_QR_QRsolve, -1);
  rb_define_module_function(mgsl_linalg_QR, "update", rb_gsl_linalg_QR_update, 4);

  rb_define_method(mgsl_linalg, "R_solve", rb_gsl_linalg_R_solve, -1);
  rb_define_method(cgsl_matrix_R, "solve", rb_gsl_linalg_R_solve, -1);
  /*
  rb_define_method(cgsl_matrix_R, "svx", rb_gsl_linalg_R_svx, -1);
  */
  rb_define_module_function(mgsl_linalg_QRPT, "decomp", rb_gsl_linalg_QRPT_decomp, -1);
  rb_define_method(cgsl_matrix, "QRPT_decomp", rb_gsl_linalg_QRPT_decomp, -1);
  rb_define_module_function(mgsl_linalg_QRPT, "decomp!", rb_gsl_linalg_QRPT_decomp_bang, -1);
  rb_define_method(cgsl_matrix, "QRPT_decomp!", rb_gsl_linalg_QRPT_decomp_bang, -1);

 rb_define_module_function(mgsl_linalg_QRPT, "decomp2", rb_gsl_linalg_QRPT_decomp2, -1);
  rb_define_method(cgsl_matrix, "QRPT_decomp2", rb_gsl_linalg_QRPT_decomp2, -1);

  rb_define_module_function(mgsl_linalg_QRPT, "solve", rb_gsl_linalg_QRPT_solve, -1);
  rb_define_method(cgsl_matrix, "QRPT_solve", rb_gsl_linalg_QRPT_solve, -1);
  rb_define_method(cgsl_matrix_QRPT, "solve", rb_gsl_linalg_QRPT_solve, -1);

  rb_define_module_function(mgsl_linalg_QRPT, "svx", rb_gsl_linalg_QRPT_svx, -1);
  rb_define_method(cgsl_matrix, "QRPT_svx", rb_gsl_linalg_QRPT_svx, -1);
  rb_define_method(cgsl_matrix_QRPT, "svx", rb_gsl_linalg_QRPT_svx, -1);

  rb_define_module_function(mgsl_linalg_QRPT, "QRsolve", rb_gsl_linalg_QRPT_QRsolve, 4);
  rb_define_module_function(mgsl_linalg_QRPT, "update", rb_gsl_linalg_QRPT_update, 5);

  rb_define_module_function(mgsl_linalg_QRPT, "Rsolve", rb_gsl_linalg_QRPT_Rsolve, -1);
  rb_define_method(cgsl_matrix_QRPT, "Rsolve", rb_gsl_linalg_QRPT_Rsolve, -1);
  rb_define_module_function(mgsl_linalg_QRPT, "Rsvx", rb_gsl_linalg_QRPT_Rsvx, -1);
  rb_define_method(cgsl_matrix_QRPT, "Rsvx", rb_gsl_linalg_QRPT_Rsvx, -1);

  /*****/
  rb_define_module_function(mgsl_linalg_SV, "decomp", rb_gsl_linalg_SV_decomp, -1);
  rb_define_method(cgsl_matrix, "SV_decomp", rb_gsl_linalg_SV_decomp, -1);
  rb_define_alias(cgsl_matrix, "SVD", "SV_decomp");
  rb_define_alias(cgsl_matrix, "svd", "SV_decomp");
  rb_define_module_function(mgsl_linalg_SV, "decomp_mod", rb_gsl_linalg_SV_decomp_mod, -1);
  rb_define_method(cgsl_matrix, "SV_decomp_mod", rb_gsl_linalg_SV_decomp_mod, -1);
  rb_define_module_function(mgsl_linalg_SV, "decomp_jacobi", rb_gsl_linalg_SV_decomp_jacobi, -1);
  rb_define_method(cgsl_matrix, "SV_decomp_jacobi", rb_gsl_linalg_SV_decomp_jacobi, -1);

  rb_define_module_function(mgsl_linalg_SV, "solve", rb_gsl_linalg_SV_solve, -1);

  rb_define_method(cgsl_matrix, "SV_solve", rb_gsl_linalg_SV_solve, -1);

  /*****/
  rb_define_module_function(mgsl_linalg_cholesky, "decomp", rb_gsl_linalg_cholesky_decomp, -1);
  rb_define_method(cgsl_matrix, "cholesky_decomp", rb_gsl_linalg_cholesky_decomp, -1);

  rb_define_module_function(mgsl_linalg_cholesky, "solve", rb_gsl_linalg_cholesky_solve, -1);
  rb_define_method(cgsl_matrix, "cholesky_solve", rb_gsl_linalg_cholesky_solve, -1);
  rb_define_method(cgsl_matrix_C, "solve", rb_gsl_linalg_cholesky_solve, -1);
  rb_define_module_function(mgsl_linalg_cholesky, "svx", rb_gsl_linalg_cholesky_svx, -1);
  rb_define_method(cgsl_matrix, "cholesky_svx", rb_gsl_linalg_cholesky_svx, -1);
  rb_define_method(cgsl_matrix_C, "svx", rb_gsl_linalg_cholesky_svx, -1);

  /*****/

  rb_define_module_function(mgsl_linalg_symmtd, "decomp", rb_gsl_linalg_symmtd_decomp, -1);
  rb_define_method(cgsl_matrix, "symmtd_decomp", rb_gsl_linalg_symmtd_decomp, -1);
  rb_define_module_function(mgsl_linalg_symmtd, "decomp!", rb_gsl_linalg_symmtd_decomp2, -1);
  rb_define_method(cgsl_matrix, "symmtd_decomp!", rb_gsl_linalg_symmtd_decomp2, -1);

  rb_define_method(cgsl_matrix, "symmtd_unpack", rb_gsl_linalg_symmtd_unpack, -1);
  rb_define_method(cgsl_matrix, "symmtd_unpack_T", rb_gsl_linalg_symmtd_unpack_T, -1);

  rb_define_module_function(mgsl_linalg_symmtd, "unpack", rb_gsl_linalg_symmtd_unpack, -1);
  rb_define_module_function(mgsl_linalg_symmtd, "unpack_T", rb_gsl_linalg_symmtd_unpack_T, -1);
  /*****/
  rb_define_module_function(mgsl_linalg_hermtd, "decomp", rb_gsl_linalg_hermtd_decomp, -1);
  rb_define_method(cgsl_matrix, "hermtd_decomp", rb_gsl_linalg_hermtd_decomp, -1);
   rb_define_module_function(mgsl_linalg_hermtd, "decomp!", rb_gsl_linalg_hermtd_decomp2, -1);
  rb_define_method(cgsl_matrix, "hermtd_decomp!", rb_gsl_linalg_hermtd_decomp2, -1);
 
  rb_define_method(cgsl_matrix_complex, "hermtd_unpack", rb_gsl_linalg_hermtd_unpack, -1);
  rb_define_module_function(mgsl_linalg_hermtd, "unpack", rb_gsl_linalg_hermtd_unpack, -1);
  rb_define_method(cgsl_matrix_complex, "hermtd_unpack_T", rb_gsl_linalg_hermtd_unpack_T, -1);
  rb_define_module_function(mgsl_linalg_hermtd, "unpack_T", rb_gsl_linalg_hermtd_unpack_T, -1);

  /*****/
  rb_define_method(cgsl_matrix, "bidiag_decomp", rb_gsl_linalg_bidiag_decomp, -1);
  rb_define_method(cgsl_matrix, "bidiag_decomp!", rb_gsl_linalg_bidiag_decomp2, -1);

  rb_define_module_function(mgsl_linalg, "bidiag_decomp", rb_gsl_linalg_bidiag_decomp, -1);
  rb_define_module_function(mgsl_linalg, "bidiag_decomp!", rb_gsl_linalg_bidiag_decomp2, -1);
  rb_define_module_function(mgsl_linalg_bidiag, "decomp", rb_gsl_linalg_bidiag_decomp, -1);
  rb_define_module_function(mgsl_linalg_bidiag, "decomp!", rb_gsl_linalg_bidiag_decomp2, -1);

  rb_define_method(cgsl_matrix, "bidiag_unpack", rb_gsl_linalg_bidiag_unpack, -1);
  rb_define_method(cgsl_matrix, "bidiag_unpack2", rb_gsl_linalg_bidiag_unpack2, -1);
  rb_define_module_function(mgsl_linalg, "bidiag_unpack", rb_gsl_linalg_bidiag_unpack, -1);
  rb_define_module_function(mgsl_linalg, "bidiag_unpack2", rb_gsl_linalg_bidiag_unpack2, -1);
  rb_define_module_function(mgsl_linalg_bidiag, "unpack", rb_gsl_linalg_bidiag_unpack, -1);
  rb_define_module_function(mgsl_linalg_bidiag, "unpack2", rb_gsl_linalg_bidiag_unpack2, -1);

  rb_define_method(cgsl_matrix, "bidiag_unpack_B", rb_gsl_linalg_bidiag_unpack_B, -1);
  rb_define_module_function(mgsl_linalg, "bidiag_unpack_B", rb_gsl_linalg_bidiag_unpack_B, -1);
  rb_define_module_function(mgsl_linalg_bidiag, "unpack_B", rb_gsl_linalg_bidiag_unpack_B, -1);
  /*****/
  rb_define_module_function(mgsl_linalg, "householder_transform", 
			     rb_gsl_linalg_householder_transform, -1);
  rb_define_module_function(mgsl_linalg_Householder, "transform", 
			     rb_gsl_linalg_householder_transform, -1);
  rb_define_module_function(mgsl_linalg_HH, "transform", 
			     rb_gsl_linalg_householder_transform, -1);
  rb_define_method(cgsl_vector, "householder_transform", 
		   rb_gsl_linalg_householder_transform, -1);

  rb_define_module_function(mgsl_linalg, "householder_hm", 
			     rb_gsl_linalg_householder_hm, 3);
  rb_define_module_function(mgsl_linalg_Householder, "hm", 
			     rb_gsl_linalg_householder_hm, 3);
  rb_define_module_function(mgsl_linalg_HH, "hm", 
			     rb_gsl_linalg_householder_hm, 3);

  rb_define_module_function(mgsl_linalg, "householder_mh", 
			     rb_gsl_linalg_householder_mh, 3);
  rb_define_module_function(mgsl_linalg_Householder, "mh", 
			     rb_gsl_linalg_householder_mh, 3);
  rb_define_module_function(mgsl_linalg_HH, "mh", 
			     rb_gsl_linalg_householder_mh, 3);

  rb_define_module_function(mgsl_linalg, "householder_hv", 
			     rb_gsl_linalg_householder_hv, 3);
  rb_define_module_function(mgsl_linalg_Householder, "hv", 
			     rb_gsl_linalg_householder_hv, 3);
  rb_define_module_function(mgsl_linalg_HH, "hv", 
			     rb_gsl_linalg_householder_hv, 3);

  rb_define_module_function(mgsl_linalg_HH, "solve", rb_gsl_linalg_HH_solve, -1);
  rb_define_module_function(mgsl_linalg_HH, "solve!", rb_gsl_linalg_HH_solve_bang, -1);
  rb_define_method(cgsl_matrix, "HH_solve", rb_gsl_linalg_HH_solve, -1);
  rb_define_method(cgsl_matrix, "HH_solve!", rb_gsl_linalg_HH_solve_bang, -1);
  rb_define_module_function(mgsl_linalg_HH, "svx", rb_gsl_linalg_HH_svx, -1);
  rb_define_method(cgsl_matrix, "HH_svx", rb_gsl_linalg_HH_svx, -1);

  /*****/

  rb_define_module_function(mgsl_linalg, "solve_symm_tridiag", rb_gsl_linalg_solve_symm_tridiag, 3);

  rb_define_module_function(mgsl_linalg_tridiag, "solve_symm", rb_gsl_linalg_solve_symm_tridiag, 3);

#ifdef GSL_1_2_LATER
  rb_define_module_function(mgsl_linalg, "solve_tridiag", rb_gsl_linalg_solve_tridiag, 4);
  rb_define_module_function(mgsl_linalg_tridiag, "solve", rb_gsl_linalg_solve_tridiag, 4);
  rb_define_module_function(mgsl_linalg, "solve_symm_cyc_tridiag", rb_gsl_linalg_solve_symm_cyc_tridiag, 3);
  rb_define_module_function(mgsl_linalg, "solve_cyc_tridiag", rb_gsl_linalg_solve_cyc_tridiag, 4);
  rb_define_module_function(mgsl_linalg_tridiag, "solve_symm_cyc", rb_gsl_linalg_solve_symm_cyc_tridiag, 3);
  rb_define_module_function(mgsl_linalg_tridiag, "solve_cyc", rb_gsl_linalg_solve_cyc_tridiag, 4);
#endif

  /*****/
  rb_define_module_function(mgsl_linalg, "balance_columns!", 
			     rb_gsl_linalg_balance_columns_bang, -1);
  rb_define_method(cgsl_matrix, "balance_columns!", 
			     rb_gsl_linalg_balance_columns_bang, -1);
  rb_define_module_function(mgsl_linalg, "balance_columns", 
			     rb_gsl_linalg_balance_columns, -1);
  rb_define_method(cgsl_matrix, "balance_columns", 
			     rb_gsl_linalg_balance_columns, -1);
  rb_define_alias(cgsl_matrix, "balance", "balance_columns");
  rb_define_alias(cgsl_matrix, "balanc", "balance_columns");
  /*****/

  Init_gsl_linalg_complex(mgsl_linalg);			     

  /** GSL-1.6 **/
#ifdef GSL_1_6_LATER
  rb_define_module_function(mgsl_linalg, "LQ_decomp", rb_gsl_linalg_LQ_decomp, -1);
  rb_define_module_function(mgsl_linalg_LQ, "decomp", rb_gsl_linalg_LQ_decomp, -1);
  rb_define_method(cgsl_matrix, "LQ_decomp", rb_gsl_linalg_LQ_decomp, -1);
  rb_define_module_function(mgsl_linalg, "LQ_decomp!", rb_gsl_linalg_LQ_decomp_bang, -1);
  rb_define_module_function(mgsl_linalg_LQ, "decomp!", rb_gsl_linalg_LQ_decomp_bang, -1);
  rb_define_method(cgsl_matrix, "LQ_decomp!", rb_gsl_linalg_LQ_decomp_bang, -1);

  rb_define_module_function(mgsl_linalg, "LQ_solve_T", rb_gsl_linalg_LQ_solve, -1);
  rb_define_module_function(mgsl_linalg_LQ, "solve_T", rb_gsl_linalg_LQ_solve, -1);
  rb_define_module_function(mgsl_linalg, "LQ_svx_T", rb_gsl_linalg_LQ_svx, -1);
  rb_define_module_function(mgsl_linalg_LQ, "svx_T", rb_gsl_linalg_LQ_svx, -1);
  rb_define_method(cgsl_matrix, "LQ_solve_T", rb_gsl_linalg_LQ_solve, -1);
  rb_define_method(cgsl_matrix_LQ, "solve_T", rb_gsl_linalg_LQ_solve, -1);
  rb_define_method(cgsl_matrix, "LQ_svx_T", rb_gsl_linalg_LQ_svx, -1);
  rb_define_method(cgsl_matrix_LQ, "svx_T", rb_gsl_linalg_LQ_svx, -1);

  rb_define_module_function(mgsl_linalg_LQ, "lssolve_T", rb_gsl_linalg_LQ_lssolve, -1);
  rb_define_method(cgsl_matrix, "LQ_lssolve_T", rb_gsl_linalg_LQ_lssolve, -1);
  rb_define_method(cgsl_matrix_LQ, "lssolve_T", rb_gsl_linalg_LQ_lssolve, -1);

  rb_define_module_function(mgsl_linalg_LQ, "vecQT", rb_gsl_linalg_LQ_vecQT, -1);
  rb_define_method(cgsl_matrix_LQ, "vecQT", rb_gsl_linalg_LQ_vecQT, -1);
  rb_define_module_function(mgsl_linalg_LQ, "vecQ", rb_gsl_linalg_LQ_vecQ, -1);
  rb_define_method(cgsl_matrix_LQ, "vecQ", rb_gsl_linalg_LQ_vecQ, -1);

  rb_define_module_function(mgsl_linalg_LQ, "Lsolve_T", rb_gsl_linalg_LQ_Lsolve, -1);
  rb_define_method(cgsl_matrix, "LQ_Lsolve_T", rb_gsl_linalg_LQ_Lsolve, -1);
  rb_define_method(cgsl_matrix_LQ, "Lsolve_T", rb_gsl_linalg_LQ_Lsolve, -1);

  rb_define_module_function(mgsl_linalg_LQ, "Lsvx_T", rb_gsl_linalg_LQ_Lsvx, -1);
  rb_define_method(cgsl_matrix_LQ, "Lsvx_T", rb_gsl_linalg_LQ_Lsvx, 1);

  rb_define_module_function(mgsl_linalg_LQ, "unpack", rb_gsl_linalg_LQ_unpack, -1);
  rb_define_method(cgsl_matrix_LQ, "unpack", rb_gsl_linalg_LQ_unpack, -1);

  rb_define_module_function(mgsl_linalg_LQ, "LQsolve_T", rb_gsl_linalg_LQ_LQsolve, -1);
  rb_define_module_function(mgsl_linalg_LQ, "update", rb_gsl_linalg_LQ_update, 4);

  rb_define_method(mgsl_linalg, "L_solve_T", rb_gsl_linalg_L_solve, -1);
  rb_define_method(cgsl_matrix_L, "solve_T", rb_gsl_linalg_L_solve, -1);
  /*
  rb_define_method(cgsl_matrix_R, "svx", rb_gsl_linalg_R_svx, -1);
  */
  rb_define_module_function(mgsl_linalg_PTLQ, "decomp", rb_gsl_linalg_PTLQ_decomp, -1);
  rb_define_method(cgsl_matrix, "PTLQ_decomp", rb_gsl_linalg_PTLQ_decomp, -1);
  rb_define_module_function(mgsl_linalg_PTLQ, "decomp!", rb_gsl_linalg_PTLQ_decomp_bang, -1);
  rb_define_method(cgsl_matrix, "PTLQ_decomp!", rb_gsl_linalg_PTLQ_decomp_bang, -1);

 rb_define_module_function(mgsl_linalg_PTLQ, "decomp2", rb_gsl_linalg_PTLQ_decomp2, -1);
  rb_define_method(cgsl_matrix, "PTLQ_decomp2", rb_gsl_linalg_PTLQ_decomp2, -1);

  rb_define_module_function(mgsl_linalg_PTLQ, "solve_T", rb_gsl_linalg_PTLQ_solve, -1);
  rb_define_method(cgsl_matrix, "PTLQ_solve_T", rb_gsl_linalg_PTLQ_solve, -1);
  rb_define_method(cgsl_matrix_PTLQ, "solve_T", rb_gsl_linalg_PTLQ_solve, -1);

  rb_define_module_function(mgsl_linalg_PTLQ, "svx_T", rb_gsl_linalg_PTLQ_svx, -1);
  rb_define_method(cgsl_matrix, "PTLQ_svx_T", rb_gsl_linalg_PTLQ_svx, -1);
  rb_define_method(cgsl_matrix_PTLQ, "svx_T", rb_gsl_linalg_PTLQ_svx, -1);

  rb_define_module_function(mgsl_linalg_PTLQ, "LQsolve_T", rb_gsl_linalg_PTLQ_LQsolve, 4);
  rb_define_module_function(mgsl_linalg_PTLQ, "update", rb_gsl_linalg_PTLQ_update, 5);

  rb_define_module_function(mgsl_linalg_PTLQ, "Lsolve_T", rb_gsl_linalg_PTLQ_Lsolve, -1);
  rb_define_method(cgsl_matrix_PTLQ, "Lsolve_T", rb_gsl_linalg_PTLQ_Lsolve, -1);
  rb_define_module_function(mgsl_linalg_PTLQ, "Lsvx_T", rb_gsl_linalg_PTLQ_Lsvx, -1);
  rb_define_method(cgsl_matrix_PTLQ, "Lsvx_T", rb_gsl_linalg_PTLQ_Lsvx, -1);

#endif

#ifdef GSL_1_9_LATER
  mhessen = rb_define_module_under(mgsl_linalg, "Hessenberg");
  rb_define_module_function(mhessen, "decomp", rb_gsl_linalg_hessenberg_decomp, 1);
  rb_define_module_function(mgsl_linalg, "heesenberg_decomp", rb_gsl_linalg_hessenberg_decomp, 1);  
  rb_define_module_function(mhessen, "unpack", rb_gsl_linalg_hessenberg_unpack, 2);  
  rb_define_module_function(mgsl_linalg, "hessenberg_unpack", rb_gsl_linalg_hessenberg_unpack, 2);    
  rb_define_module_function(mhessen, "unpack_accum", rb_gsl_linalg_hessenberg_unpack_accum, -1);    
  rb_define_module_function(mgsl_linalg, "hessenberg_unpack_accum", rb_gsl_linalg_hessenberg_unpack_accum, -1);      
  rb_define_module_function(mhessen, "set_zero", rb_gsl_linalg_hessenberg_set_zero, 1);    
  rb_define_module_function(mgsl_linalg, "hessenberg_set_zero", rb_gsl_linalg_hessenberg_set_zero, 1);      
  
  rb_define_module_function(mgsl_linalg, "hesstri_decomp", rb_gsl_linalg_hesstri_decomp, -1);
  rb_define_module_function(mgsl_linalg, "hesstri_decomp!", rb_gsl_linalg_hesstri_decomp_bang, -1);
  
  rb_define_module_function(mgsl_linalg, "balance_matrix", rb_gsl_linalg_balance_matrix, -1);  
  rb_define_module_function(mgsl_linalg, "balance_matrix!", rb_gsl_linalg_balance_matrix2, -1);    
  rb_define_module_function(mgsl_linalg, "balance", rb_gsl_linalg_balance_matrix, -1);    
  rb_define_module_function(mgsl_linalg, "balance!", rb_gsl_linalg_balance_matrix2, -1);      
#endif

}
