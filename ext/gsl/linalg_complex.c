/*
  linalg_complex.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include <gsl/gsl_math.h>
#include "rb_gsl_array.h"
#include "rb_gsl_common.h"
#include "rb_gsl_linalg.h"

EXTERN VALUE mgsl_linalg;
EXTERN VALUE cgsl_complex;

static VALUE cgsl_matrix_complex_LU;
#ifdef GSL_1_10_LATER
static VALUE cgsl_matrix_complex_C;
#endif

VALUE rb_gsl_linalg_complex_LU_decomp(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  gsl_permutation *p = NULL;
  int signum, itmp;
  size_t size;
  VALUE obj2;

  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
    CHECK_MATRIX_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, m);
    itmp = 1;
    break;
  default:
    CHECK_MATRIX_COMPLEX(obj);
    Data_Get_Struct(obj, gsl_matrix_complex, m);
    itmp = 0;
  }
  size = m->size1;
  switch (argc-itmp) {
  case 0:
    p = gsl_permutation_alloc(size);
    gsl_linalg_complex_LU_decomp(m, p, &signum);
    if (itmp == 1) RBGSL_SET_CLASS(argv[0], cgsl_matrix_complex_LU);
    else RBGSL_SET_CLASS(obj, cgsl_matrix_complex_LU);
    obj2 = Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p);
    return rb_ary_new3(2, obj2, INT2FIX(signum));
    break;
  case 1:  /* when a permutation object is given */
    CHECK_PERMUTATION(argv[itmp]);
    Data_Get_Struct(argv[itmp], gsl_permutation, p);
    gsl_linalg_complex_LU_decomp(m, p, &signum);
    if (itmp == 1) RBGSL_SET_CLASS(argv[0], cgsl_matrix_complex_LU);
    else RBGSL_SET_CLASS(obj, cgsl_matrix_complex_LU);
    return INT2FIX(signum);
    break;
  default:
    rb_raise(rb_eArgError, "Usage: LU_decomp!() or LU_decomp!(permutation)");
  }
}

VALUE rb_gsl_linalg_complex_LU_decomp2(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL, *mnew = NULL;
  gsl_permutation *p = NULL;
  int signum, itmp;
  size_t size;
  VALUE objm, obj2;

  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
    CHECK_MATRIX_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, m);
    itmp = 1;
    break;
  default:
    CHECK_MATRIX_COMPLEX(obj);
    Data_Get_Struct(obj, gsl_matrix_complex, m);
    itmp = 0;
  }
  size = m->size1;
  mnew = gsl_matrix_complex_alloc(m->size1, m->size2);
  gsl_matrix_complex_memcpy(mnew, m);
  objm = Data_Wrap_Struct(cgsl_matrix_complex_LU, 0, gsl_matrix_complex_free, mnew);
  switch (argc-itmp) {
  case 0:
    p = gsl_permutation_alloc(size);
    gsl_linalg_complex_LU_decomp(mnew, p, &signum);
    obj2 = Data_Wrap_Struct(cgsl_permutation, 0, gsl_permutation_free, p);
    return rb_ary_new3(3, objm ,obj2, INT2FIX(signum));
    break;
  case 1:  /* when a permutation object is given */
    CHECK_PERMUTATION(argv[itmp]);
    Data_Get_Struct(argv[itmp], gsl_permutation, p);
    gsl_linalg_complex_LU_decomp(m, p, &signum);
    return rb_ary_new3(3, objm , argv[itmp], INT2FIX(signum));
    break;
  default:
    rb_raise(rb_eArgError, "Usage: LU_decomp!() or LU_decomp!(permutation)");
  }
}

static VALUE rb_gsl_linalg_complex_LU_solve(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL, *mtmp = NULL;
  gsl_permutation *p = NULL;
  gsl_vector_complex *b = NULL, *x = NULL;
  int flagm = 0, flagx = 0, itmp, signum;
  
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 2 || argc > 4) 
      rb_raise(rb_eArgError, "Usage: solve(m, b), solve(m, b, x), solve(lu, p, b), solve(lu, p, b, x)");

    CHECK_MATRIX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, m);
    if (CLASS_OF(argv[0]) != cgsl_matrix_complex_LU) {
      flagm = 1;
      mtmp = gsl_matrix_complex_alloc(m->size1, m->size2);
      gsl_matrix_complex_memcpy(mtmp, m);
    } else {
      mtmp = m;
    }
    itmp = 1;
    break;
  default:
    if (argc < 1 || argc > 3) 
      rb_raise(rb_eArgError, "Usage: LU_solve(b), LU_solve(p, b), LU_solve(b, x), solve(p, b, x)");
    Data_Get_Struct(obj, gsl_matrix_complex, m);
    if (CLASS_OF(obj) != cgsl_matrix_complex_LU) {
      flagm = 1;
      mtmp = gsl_matrix_complex_alloc(m->size1, m->size2);
      gsl_matrix_complex_memcpy(mtmp, m);
    } else {
      mtmp = m;
    }
    itmp = 0;
  }
  if (flagm == 1) {
    if (itmp != argc-1) rb_raise(rb_eArgError, "Usage: m.LU_solve(b)");
    Data_Get_Struct(argv[itmp], gsl_vector_complex, b);
    x = gsl_vector_complex_alloc(b->size);
    p = gsl_permutation_alloc(b->size);
    gsl_linalg_complex_LU_decomp(mtmp, p, &signum);
  } else {
    Data_Get_Struct(argv[itmp], gsl_permutation, p);
    itmp++;
    Data_Get_Struct(argv[itmp], gsl_vector_complex, b);
    itmp++;
    if (itmp == argc-1) {
      Data_Get_Struct(argv[itmp], gsl_vector_complex, x);
      flagx = 1;
    } else {
      x = gsl_vector_complex_alloc(m->size1);
    }
  }
  gsl_linalg_complex_LU_solve(mtmp, p, b, x);
  if (flagm == 1) {
    gsl_matrix_complex_free(mtmp);
    gsl_permutation_free(p);
  }
  if (flagx == 0) return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, x);
  else return argv[argc-1];
}


static VALUE rb_gsl_linalg_complex_LU_svx(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL, *mtmp = NULL;
  gsl_permutation *p = NULL;
  gsl_vector_complex *x = NULL;
  int flagm = 0, itmp, signum;
  
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    CHECK_MATRIX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, m);
    if (CLASS_OF(argv[0]) != cgsl_matrix_complex_LU) {
      flagm = 1;
      mtmp = gsl_matrix_complex_alloc(m->size1, m->size2);
      gsl_matrix_complex_memcpy(mtmp, m);
    } else {
      mtmp = m;
    }
    itmp = 1;
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix_complex, m);
    if (CLASS_OF(obj) != cgsl_matrix_complex_LU) {
      flagm = 1;
      mtmp = gsl_matrix_complex_alloc(m->size1, m->size2);
      gsl_matrix_complex_memcpy(mtmp, m);
    } else {
      mtmp = m;
    }
    itmp = 0;
  }
  if (flagm == 1) {
    if (itmp != argc-1) rb_raise(rb_eArgError, "Usage: m.LU_solve(b)");
    Data_Get_Struct(argv[itmp], gsl_vector_complex, x);
    p = gsl_permutation_alloc(x->size);
    gsl_linalg_complex_LU_decomp(mtmp, p, &signum);
  } else {
    Data_Get_Struct(argv[itmp], gsl_permutation, p);
    itmp++;
    Data_Get_Struct(argv[itmp], gsl_vector_complex, x);
    itmp++;
  }
  gsl_linalg_complex_LU_svx(mtmp, p, x);
  if (flagm == 1) {
    gsl_matrix_complex_free(mtmp);
    gsl_permutation_free(p);
  }
  return argv[argc-1];
}

static VALUE rb_gsl_linalg_complex_LU_refine(VALUE obj, VALUE vm,
					     VALUE lu, VALUE pp, VALUE bb,
					     VALUE xx)
{
  gsl_matrix_complex *m = NULL, *mlu = NULL;
  gsl_permutation *p = NULL;
  gsl_vector_complex *b = NULL, *x = NULL, *r = NULL;
  int flagb = 0;
  VALUE vr;

  if (CLASS_OF(obj) != cgsl_matrix_complex_LU)
    rb_raise(rb_eRuntimeError, "Decompose first!");
  CHECK_MATRIX_COMPLEX(vm);
  CHECK_MATRIX_COMPLEX(lu);
  CHECK_PERMUTATION(pp);
  CHECK_VECTOR_COMPLEX(xx);
  Data_Get_Struct(vm, gsl_matrix_complex, m);
  Data_Get_Struct(lu, gsl_matrix_complex, mlu);
  Data_Get_Struct(pp, gsl_permutation, p);
  CHECK_VECTOR_COMPLEX(bb);
  Data_Get_Struct(bb, gsl_vector_complex, b);
  Data_Get_Struct(xx, gsl_vector_complex, x);
  r = gsl_vector_complex_alloc(m->size1);
  gsl_linalg_complex_LU_refine(m, mlu, p, b, x, r);
  vr = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, r);
  if (flagb == 1) gsl_vector_complex_free(b);
  return rb_ary_new3(2, xx, vr);
}

static VALUE rb_gsl_linalg_complex_LU_invert(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL, *mtmp = NULL, *inverse = NULL;
  gsl_permutation *p = NULL;
  int flagm = 0, signum, itmp;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    CHECK_MATRIX_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, m);
    if (CLASS_OF(argv[0]) != cgsl_matrix_complex_LU) {
      mtmp = gsl_matrix_complex_alloc(m->size1, m->size2);
      gsl_matrix_complex_memcpy(mtmp, m);
      flagm = 1;
    } else {
      mtmp = m;
    }
    itmp = 1;
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix_complex, m);
    if (CLASS_OF(obj) != cgsl_matrix_complex_LU) {
      mtmp = gsl_matrix_complex_alloc(m->size1, m->size2);
      gsl_matrix_complex_memcpy(mtmp, m);
      flagm = 1;
    } else {
      mtmp = m;
    }
    itmp = 0;
  }

  if (flagm == 1) {
    p = gsl_permutation_alloc(m->size1);
    gsl_linalg_complex_LU_decomp(mtmp, p, &signum);
  } else {
    Data_Get_Struct(argv[itmp], gsl_permutation, p);
  }
  inverse = gsl_matrix_complex_alloc(m->size1, m->size2);
  gsl_linalg_complex_LU_invert(mtmp, p, inverse);
  if (flagm == 1) {
    gsl_matrix_complex_free(mtmp);
    gsl_permutation_free(p);
  }
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, inverse);
}
#ifdef GSL_1_1_1_LATER
static VALUE rb_gsl_linalg_complex_LU_det(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL, *mtmp = NULL;
  gsl_permutation *p = NULL;
  gsl_complex *z = NULL;
  VALUE vz;
  int flagm = 0, signum, itmp;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    CHECK_MATRIX_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, m);
    if (CLASS_OF(argv[0]) != cgsl_matrix_complex_LU) {
      mtmp = gsl_matrix_complex_alloc(m->size1, m->size2);
      gsl_matrix_complex_memcpy(mtmp, m);
      flagm = 1;
    } else {
      mtmp = m;
    }
    itmp = 1;
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix_complex, m);
    if (CLASS_OF(obj) != cgsl_matrix_complex_LU) {
      mtmp = gsl_matrix_complex_alloc(m->size1, m->size2);
      gsl_matrix_complex_memcpy(mtmp, m);
      flagm = 1;
    } else {
      mtmp = m;
    }
    itmp = 0;
  }
  if (flagm == 1) {
    p = gsl_permutation_alloc(m->size1);
    gsl_linalg_complex_LU_decomp(mtmp, p, &signum);
  } else {
    if (itmp != argc-1) rb_raise(rb_eArgError, "signum not given");
    signum = NUM2DBL(argv[itmp]);
  }
  vz = Data_Make_Struct(cgsl_complex, gsl_complex, 0, free, z);
  *z = gsl_linalg_complex_LU_det(mtmp, signum);
  if (flagm == 1) {
    gsl_matrix_complex_free(mtmp);
    gsl_permutation_free(p);
  }
  return vz;
}

static VALUE rb_gsl_linalg_complex_LU_lndet(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL, *mtmp = NULL;
  gsl_permutation *p = NULL;
  double lndet;
  int flagm = 0, signum;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    CHECK_MATRIX_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, m);
    if (CLASS_OF(argv[0]) != cgsl_matrix_complex_LU) {
      mtmp = gsl_matrix_complex_alloc(m->size1, m->size2);
      gsl_matrix_complex_memcpy(mtmp, m);
      flagm = 1;
    } else {
      mtmp = m;
    }
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix_complex, m);
    if (CLASS_OF(obj) != cgsl_matrix_complex_LU) {
      mtmp = gsl_matrix_complex_alloc(m->size1, m->size2);
      gsl_matrix_complex_memcpy(mtmp, m);
      flagm = 1;
    } else {
      mtmp = m;
    }
  }
  if (flagm == 1) {
    p = gsl_permutation_alloc(m->size1);
    gsl_linalg_complex_LU_decomp(mtmp, p, &signum);
  } 
  lndet = gsl_linalg_complex_LU_lndet(mtmp);
  if (flagm == 1) {
    gsl_matrix_complex_free(mtmp);
    gsl_permutation_free(p);
  }
  return rb_float_new(lndet);
}

static VALUE rb_gsl_linalg_complex_LU_sgndet(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *m = NULL, *mtmp = NULL;
  gsl_permutation *p = NULL;
  gsl_complex *z = NULL;
  VALUE vz;
  int flagm = 0, signum, itmp;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    CHECK_MATRIX_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, m);
    if (CLASS_OF(argv[0]) != cgsl_matrix_complex_LU) {
      mtmp = gsl_matrix_complex_alloc(m->size1, m->size2);
      gsl_matrix_complex_memcpy(mtmp, m);
      flagm = 1;
    } else {
      mtmp = m;
    }
    itmp = 1;
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix_complex, m);
    if (CLASS_OF(obj) != cgsl_matrix_complex_LU) {
      mtmp = gsl_matrix_complex_alloc(m->size1, m->size2);
      gsl_matrix_complex_memcpy(mtmp, m);
      flagm = 1;
    } else {
      mtmp = m;
    }
    itmp = 0;
  }
  if (flagm == 1) {
    p = gsl_permutation_alloc(m->size1);
    gsl_linalg_complex_LU_decomp(mtmp, p, &signum);
  } else {
    if (itmp != argc-1) rb_raise(rb_eArgError, "signum not given");
    signum = NUM2DBL(argv[itmp]);
  }
  vz = Data_Make_Struct(cgsl_complex, gsl_complex, 0, free, z);
  *z = gsl_linalg_complex_LU_sgndet(mtmp, signum);
  if (flagm == 1) {
    gsl_matrix_complex_free(mtmp);
    gsl_permutation_free(p);
  }
  return vz;
}
#endif

#ifdef GSL_1_10_LATER

static VALUE rb_gsl_linalg_cholesky_decomp(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL, *Atmp = NULL;
  switch(TYPE(obj)) {
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
  gsl_linalg_complex_cholesky_decomp(A);
  return Data_Wrap_Struct(cgsl_matrix_complex_C, 0, gsl_matrix_complex_free, A);
}

static VALUE rb_gsl_linalg_cholesky_solve(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL, *Atmp = NULL;
  gsl_vector_complex *b = NULL, *x = NULL;
  int flagb = 0, flaga = 0;
  VALUE vA, vb;
  switch(TYPE(obj)) {
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
  CHECK_MATRIX_COMPLEX(vA);
  Data_Get_Struct(vA, gsl_matrix_complex, Atmp);
  CHECK_VECTOR_COMPLEX(vb);
  Data_Get_Struct(vb, gsl_vector_complex, b);

  if (CLASS_OF(vA) == cgsl_matrix_complex_C) {
    A = Atmp;
  } else {
    A = make_matrix_complex_clone(Atmp);
    flaga = 1;
    gsl_linalg_complex_cholesky_decomp(A);
  }
  x = gsl_vector_complex_alloc(b->size);
  gsl_linalg_complex_cholesky_solve(A, b, x);
  if (flaga == 1) gsl_matrix_complex_free(A);
  if (flagb == 1) gsl_vector_complex_free(b);
  return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, x);
}


static VALUE rb_gsl_linalg_cholesky_svx(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *A = NULL, *Atmp = NULL;
  gsl_vector_complex *b = NULL;
  int flaga = 0;
  VALUE vA, vb;
  switch(TYPE(obj)) {
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
  CHECK_MATRIX_COMPLEX(vA);
  Data_Get_Struct(vA, gsl_matrix_complex, Atmp);
  CHECK_VECTOR_COMPLEX(vb);
  Data_Get_Struct(vb, gsl_vector_complex, b);
  if (CLASS_OF(vA) == cgsl_matrix_complex_C) {
    A = Atmp;
  } else {
    A = make_matrix_complex_clone(Atmp);
    flaga = 1;
    gsl_linalg_complex_cholesky_decomp(A);
  }
  gsl_linalg_complex_cholesky_svx(A, b);
  if (flaga == 1) gsl_matrix_complex_free(A);
  return vb;
}


static VALUE rb_gsl_linalg_complex_householder_transform(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *v = NULL;
  gsl_complex *z;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments.");
    CHECK_VECTOR_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector_complex, v);
    break;
  default:
    Data_Get_Struct(obj, gsl_vector_complex, v);
    break;
  }
  z = (gsl_complex*) malloc(sizeof(gsl_complex));
  *z = gsl_linalg_complex_householder_transform(v);
  return Data_Wrap_Struct(cgsl_complex, 0, free, z);
}

/* singleton */
static VALUE rb_gsl_linalg_complex_householder_hm(VALUE obj, VALUE t, VALUE vv, VALUE aa)
{
  gsl_vector_complex *v = NULL;
  gsl_complex *tau;
  gsl_matrix_complex *A = NULL;
  CHECK_COMPLEX(t);  
  CHECK_VECTOR_COMPLEX(vv);
  CHECK_MATRIX_COMPLEX(aa);
  Data_Get_Struct(t, gsl_complex, tau);
  Data_Get_Struct(vv, gsl_vector_complex, v);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  gsl_linalg_complex_householder_hm(*tau, v, A);
  return aa;
}

static VALUE rb_gsl_linalg_complex_householder_mh(VALUE obj, VALUE t, VALUE vv, VALUE aa)
{
  gsl_vector_complex *v = NULL;
  gsl_complex *tau;
  gsl_matrix_complex *A = NULL;
  CHECK_COMPLEX(t);  
  CHECK_VECTOR_COMPLEX(vv);
  CHECK_MATRIX_COMPLEX(aa);
  Data_Get_Struct(t, gsl_complex, tau);
  Data_Get_Struct(vv, gsl_vector_complex, v);
  Data_Get_Struct(aa, gsl_matrix_complex, A);
  gsl_linalg_complex_householder_hm(*tau, v, A);
  return aa;
}

static VALUE rb_gsl_linalg_complex_householder_hv(VALUE obj, VALUE t, VALUE vv, VALUE ww)
{
  gsl_vector_complex *v = NULL, *w = NULL;
  gsl_complex *tau;
  CHECK_COMPLEX(t);    
  CHECK_VECTOR_COMPLEX(vv);
  CHECK_VECTOR_COMPLEX(ww);
  Data_Get_Struct(t, gsl_complex, tau);
  Data_Get_Struct(vv, gsl_vector_complex, v);
  Data_Get_Struct(ww, gsl_vector_complex, w);
  gsl_linalg_complex_householder_hv(*tau, v, w);
  return ww;
}
#endif

void Init_gsl_linalg_complex(VALUE module)
{
  VALUE mgsl_linalg_complex;
  VALUE mgsl_linalg_complex_LU;
#ifdef GSL_1_10_LATER
  VALUE mgsl_linalg_complex_chol, mgsl_linalg_complex_Householder;
#endif

  mgsl_linalg_complex = rb_define_module_under(module, "Complex");
  mgsl_linalg_complex_LU = rb_define_module_under(mgsl_linalg_complex, "LU");

  cgsl_matrix_complex_LU = rb_define_class_under(mgsl_linalg_complex_LU,
						 "LUMatrix", cgsl_matrix_complex);

  rb_define_singleton_method(mgsl_linalg_complex, "LU_decomp!",
			     rb_gsl_linalg_complex_LU_decomp, -1);
  rb_define_singleton_method(mgsl_linalg_complex_LU, "decomp!",
			     rb_gsl_linalg_complex_LU_decomp, -1);
  rb_define_method(cgsl_matrix_complex, "LU_decomp!", 
		   rb_gsl_linalg_complex_LU_decomp, -1);
  rb_define_alias(cgsl_matrix_complex, "decomp!", "LU_decomp!");

  rb_define_singleton_method(mgsl_linalg_complex, "LU_decomp",
			     rb_gsl_linalg_complex_LU_decomp2, -1);
  rb_define_singleton_method(mgsl_linalg_complex_LU, "decomp",
			     rb_gsl_linalg_complex_LU_decomp2, -1);
  rb_define_method(cgsl_matrix_complex, "LU_decomp", 
		   rb_gsl_linalg_complex_LU_decomp2, -1);
  rb_define_alias(cgsl_matrix_complex, "decomp", "LU_decomp");

  rb_define_singleton_method(mgsl_linalg_complex, "LU_solve",
			     rb_gsl_linalg_complex_LU_solve, -1);
  rb_define_singleton_method(mgsl_linalg_complex_LU, "solve",
			     rb_gsl_linalg_complex_LU_solve, -1);
  rb_define_method(cgsl_matrix_complex, "LU_solve", 
		   rb_gsl_linalg_complex_LU_solve, -1);
  rb_define_method(cgsl_matrix_complex_LU, "solve", 
		   rb_gsl_linalg_complex_LU_solve, -1);

  rb_define_singleton_method(mgsl_linalg_complex, "LU_svx",
			     rb_gsl_linalg_complex_LU_svx, -1);
  rb_define_singleton_method(mgsl_linalg_complex_LU, "svx",
			     rb_gsl_linalg_complex_LU_svx, -1);
  rb_define_method(cgsl_matrix_complex, "LU_svx", 
		   rb_gsl_linalg_complex_LU_svx, -1);
  rb_define_method(cgsl_matrix_complex_LU, "svx", 
		   rb_gsl_linalg_complex_LU_svx, -1);

  rb_define_singleton_method(mgsl_linalg_complex, "LU_refine",
			     rb_gsl_linalg_complex_LU_refine, 5);
  rb_define_singleton_method(mgsl_linalg_complex_LU, "refine",
			     rb_gsl_linalg_complex_LU_refine, 5);

  rb_define_singleton_method(mgsl_linalg_complex, "LU_invert",
			     rb_gsl_linalg_complex_LU_invert, -1);
  rb_define_singleton_method(mgsl_linalg_complex_LU, "invert",
			     rb_gsl_linalg_complex_LU_invert, -1);
  rb_define_method(cgsl_matrix_complex, "LU_invert", 
		   rb_gsl_linalg_complex_LU_invert, -1);
  rb_define_alias(cgsl_matrix_complex, "invert", "LU_invert");
  rb_define_alias(cgsl_matrix_complex, "inv", "LU_invert");
  rb_define_method(cgsl_matrix_complex_LU, "invert", 
		   rb_gsl_linalg_complex_LU_invert, -1);

#ifdef GSL_1_1_1_LATER
  rb_define_singleton_method(mgsl_linalg_complex, "LU_det", 
			     rb_gsl_linalg_complex_LU_det, -1);
  rb_define_singleton_method(mgsl_linalg_complex_LU, "det", 
			     rb_gsl_linalg_complex_LU_det, -1);
  rb_define_method(cgsl_matrix_complex, "LU_det", rb_gsl_linalg_complex_LU_det, -1);
  rb_define_alias(cgsl_matrix_complex, "det", "LU_det");
  rb_define_method(cgsl_matrix_complex_LU, "det", rb_gsl_linalg_complex_LU_det, -1);

  rb_define_singleton_method(mgsl_linalg_complex, "LU_lndet", 
			     rb_gsl_linalg_complex_LU_lndet, -1);
  rb_define_singleton_method(mgsl_linalg_complex_LU, "lndet", 
			     rb_gsl_linalg_complex_LU_lndet, -1);
  rb_define_method(cgsl_matrix_complex, "LU_lndet", rb_gsl_linalg_complex_LU_lndet, -1);
  rb_define_alias(cgsl_matrix_complex, "lndet", "LU_lndet");
  rb_define_method(cgsl_matrix_complex_LU, "LU_lndet", rb_gsl_linalg_complex_LU_lndet, -1);

  rb_define_singleton_method(mgsl_linalg_complex, "LU_sgndet", 
			     rb_gsl_linalg_complex_LU_sgndet, -1);
  rb_define_singleton_method(mgsl_linalg_complex_LU, "sgndet", 
			     rb_gsl_linalg_complex_LU_sgndet, -1);
  rb_define_method(cgsl_matrix_complex, "LU_sgndet", rb_gsl_linalg_complex_LU_sgndet, -1);
  rb_define_alias(cgsl_matrix_complex, "sgndet", "LU_sgndet");
  rb_define_method(cgsl_matrix_complex_LU, "LU_sgndet", rb_gsl_linalg_complex_LU_sgndet, -1);
#endif

#ifdef GSL_1_10_LATER
  mgsl_linalg_complex_chol = rb_define_module_under(mgsl_linalg_complex, "Cholesky");
  cgsl_matrix_complex_C = rb_define_class_under(mgsl_linalg_complex_chol, "CholeskyMatrix", cgsl_matrix_complex);    
  rb_define_singleton_method(mgsl_linalg_complex_chol, "decomp", rb_gsl_linalg_cholesky_decomp, -1);
  rb_define_method(cgsl_matrix_complex, "cholesky_decomp", rb_gsl_linalg_cholesky_decomp, -1);  
  rb_define_singleton_method(mgsl_linalg_complex_chol, "solve", rb_gsl_linalg_cholesky_solve, -1);
  rb_define_method(cgsl_matrix_complex, "cholesky_solve", rb_gsl_linalg_cholesky_solve, -1);
  rb_define_method(cgsl_matrix_complex_C, "solve", rb_gsl_linalg_cholesky_solve, -1);
  rb_define_singleton_method(mgsl_linalg_complex_chol, "svx", rb_gsl_linalg_cholesky_svx, -1);
  rb_define_method(cgsl_matrix_complex, "cholesky_svx", rb_gsl_linalg_cholesky_svx, -1);
  rb_define_method(cgsl_matrix_complex_C, "svx", rb_gsl_linalg_cholesky_svx, -1);

  mgsl_linalg_complex_Householder = rb_define_module_under(mgsl_linalg_complex, "Householder");
 rb_define_singleton_method(mgsl_linalg_complex, "householder_transform", 
			     rb_gsl_linalg_complex_householder_transform, -1);
  rb_define_singleton_method(mgsl_linalg_complex_Householder, "transform", 
			     rb_gsl_linalg_complex_householder_transform, -1);
  rb_define_method(cgsl_vector_complex, "householder_transform", 
		   rb_gsl_linalg_complex_householder_transform, -1);

  rb_define_singleton_method(mgsl_linalg_complex, "householder_hm", 
			     rb_gsl_linalg_complex_householder_hm, 3);
  rb_define_singleton_method(mgsl_linalg_complex_Householder, "hm", 
			     rb_gsl_linalg_complex_householder_hm, 3);

  rb_define_singleton_method(mgsl_linalg_complex, "householder_mh", 
			     rb_gsl_linalg_complex_householder_mh, 3);
  rb_define_singleton_method(mgsl_linalg_complex_Householder, "mh", 
			     rb_gsl_linalg_complex_householder_mh, 3);

  rb_define_singleton_method(mgsl_linalg_complex, "householder_hv", 
			     rb_gsl_linalg_complex_householder_hv, 3);
  rb_define_singleton_method(mgsl_linalg_complex_Householder, "hv", 
			     rb_gsl_linalg_complex_householder_hv, 3);
#endif			     
}
