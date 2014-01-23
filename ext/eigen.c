/*
  eigen.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl.h"
#include "rb_gsl_array.h"

#include "rb_gsl_eigen.h"
#include "rb_gsl_complex.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

static VALUE cgsl_eigen_symm_workspace;
static VALUE cgsl_eigen_symmv_workspace;
static VALUE cgsl_eigen_herm_workspace;
static VALUE cgsl_eigen_hermv_workspace;
static VALUE cgsl_eigen_values;
static VALUE cgsl_eigen_vectors;
static VALUE cgsl_eigen_vector;
static VALUE cgsl_eigen_vector_complex;
static VALUE cgsl_eigen_herm_vectors;

#ifdef HAVE_EIGEN_FRANCIS
static VALUE cgsl_eigen_francis_workspace;

#endif
#ifdef GSL_1_9_LATER
static VALUE cgsl_eigen_nonsymm_workspace;
static VALUE cgsl_eigen_nonsymmv_workspace;
#endif

#ifdef GSL_1_10_LATER
static VALUE cgensymm, mgensymm;
static VALUE cgensymmv, mgensymmv;
static VALUE cgenherm, mgenherm;
static VALUE cgenhermv, mgenhermv;

static VALUE mgen, mgenv;
static VALUE cgenw, cgenvw;
#endif

static VALUE rb_gsl_eigen_symm_alloc(VALUE klass, VALUE nn)
{
  gsl_eigen_symm_workspace *w = NULL;
  CHECK_FIXNUM(nn);
  w = gsl_eigen_symm_alloc(FIX2INT(nn));
  return Data_Wrap_Struct(cgsl_eigen_symm_workspace, 0, 
			  gsl_eigen_symm_free, w);
}

static VALUE rb_gsl_eigen_symmv_alloc(VALUE klass, VALUE nn)
{
  gsl_eigen_symmv_workspace *w = NULL;
  CHECK_FIXNUM(nn);
  w = gsl_eigen_symmv_alloc(FIX2INT(nn));
  return Data_Wrap_Struct(cgsl_eigen_symmv_workspace, 0, gsl_eigen_symmv_free, w);
}

static VALUE rb_gsl_eigen_herm_alloc(VALUE klass, VALUE nn)
{
  gsl_eigen_herm_workspace *w = NULL;
  CHECK_FIXNUM(nn);
  w = gsl_eigen_herm_alloc(FIX2INT(nn));
  return Data_Wrap_Struct(cgsl_eigen_herm_workspace, 0, 
			  gsl_eigen_herm_free, w);
}

static VALUE rb_gsl_eigen_hermv_alloc(VALUE klass, VALUE nn)
{
  gsl_eigen_hermv_workspace *w = NULL;
  CHECK_FIXNUM(nn);
  w = gsl_eigen_hermv_alloc(FIX2INT(nn));
  return Data_Wrap_Struct(cgsl_eigen_hermv_workspace, 0, 
			  gsl_eigen_hermv_free, w);
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_eigen_symm_narray(int argc, VALUE *argv, VALUE obj);
#endif

static VALUE rb_gsl_eigen_symm(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *Atmp = NULL, *A = NULL;
  gsl_eigen_symm_workspace *w = NULL;
  gsl_vector *v = NULL;
  int flagw = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    switch (argc) {
    case 2:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0])) return rb_gsl_eigen_symm_narray(argc, argv, obj);
#endif
      CHECK_MATRIX(argv[0]);
      Data_Get_Struct(argv[0], gsl_matrix, Atmp);
      if (CLASS_OF(argv[1]) != cgsl_eigen_symm_workspace)
	rb_raise(rb_eTypeError,
		 "argv[1]: wrong argument type %s (Eigen::Symm::Workspace expected)", 
		 rb_class2name(CLASS_OF(argv[1])));
      Data_Get_Struct(argv[1], gsl_eigen_symm_workspace, w);
      break;
    case 1:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0])) return rb_gsl_eigen_symm_narray(argc, argv, obj);
#endif
      CHECK_MATRIX(argv[0]);
      Data_Get_Struct(argv[0], gsl_matrix, Atmp);
      w = gsl_eigen_symm_alloc(Atmp->size1);
      flagw = 1;
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
      break;
    }
    break;
  default:
    CHECK_MATRIX(obj);
    Data_Get_Struct(obj, gsl_matrix, Atmp);
    switch (argc) {
    case 1:
      if (CLASS_OF(argv[0]) != cgsl_eigen_symm_workspace)
	rb_raise(rb_eTypeError, 
		 "argv[0]:  wrong argument type %s (Eigen::Symm::Workspace expected", 
		 rb_class2name(CLASS_OF(argv[0])));

      Data_Get_Struct(argv[0], gsl_eigen_symm_workspace, w);
      break;
    case 0:
      w = gsl_eigen_symm_alloc(Atmp->size1);
      flagw = 1;
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    }
  }
  A = make_matrix_clone(Atmp);
  v = gsl_vector_alloc(A->size1);
  gsl_eigen_symm(A, v, w);
  /*  gsl_sort_vector(v);*/
  gsl_matrix_free(A);
  if (flagw == 1) gsl_eigen_symm_free(w);
  return Data_Wrap_Struct(cgsl_eigen_values, 0, gsl_vector_free, v);
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_eigen_symm_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na;
  VALUE nary;
  gsl_matrix *A = NULL;
  gsl_eigen_symm_workspace *w = NULL;
  gsl_vector_view vv;
  int shape[1];
  int flagw = 0;
  switch (argc) {
  case 2:
    if (!NA_IsNArray(argv[0])) 
      rb_raise(rb_eTypeError, "wrong argument type %s (NArray expected)", 
	       rb_class2name(CLASS_OF(argv[0])));
    GetNArray(argv[0], na);
    if (na->rank < 2) rb_raise(rb_eRuntimeError, "rank >= 2 required");
    if (na->shape[0] != na->shape[1])
      rb_raise(rb_eRuntimeError, "square matrix required");
    A = gsl_matrix_alloc(na->shape[1], na->shape[0]);
    memcpy(A->data, (double*) na->ptr, sizeof(double)*A->size1*A->size2);
    if (CLASS_OF(argv[1]) != cgsl_eigen_symm_workspace)
      rb_raise(rb_eTypeError, 
	       "argv[1]:  wrong argument type %s (Eigen::Symm::Workspace expected", 
	       rb_class2name(CLASS_OF(argv[1])));
    Data_Get_Struct(argv[1], gsl_eigen_symm_workspace, w);
    flagw = 0;
    break;
  case 1:
    if (!NA_IsNArray(argv[0])) 
      rb_raise(rb_eTypeError, "wrong argument type %s (NArray expected)", 
	       rb_class2name(CLASS_OF(argv[0])));
    GetNArray(argv[0], na);
    if (na->rank < 2) rb_raise(rb_eRuntimeError, "rank >= 2 required");
    if (na->shape[0] != na->shape[1])
      rb_raise(rb_eRuntimeError, "square matrix required");
    A = gsl_matrix_alloc(na->shape[1], na->shape[0]);
    memcpy(A->data, (double*) na->ptr, sizeof(double)*A->size1*A->size2);
    w = gsl_eigen_symm_alloc(A->size1);
    flagw = 1;
    break;
  default:
    rb_raise(rb_eArgError, "matrix not given");
    break;
  }
  shape[0] = A->size1;
  nary = na_make_object(NA_DFLOAT, 1, shape, cNVector);
  vv = gsl_vector_view_array(NA_PTR_TYPE(nary,double*), A->size1);
  gsl_eigen_symm(A, &vv.vector, w);
  /*  gsl_sort_vector(v);*/
  gsl_matrix_free(A);
  if (flagw == 1) gsl_eigen_symm_free(w);
  return nary;
}
#endif

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_eigen_symmv_narray(int argc, VALUE *argv, VALUE obj);
#endif

static VALUE rb_gsl_eigen_symmv(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *Atmp = NULL, *A = NULL, *em = NULL;
  gsl_eigen_symmv_workspace *w = NULL;
  gsl_vector *v = NULL;
  int flagw = 0;
  VALUE vval, vvec;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    switch (argc) {
    case 2:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0])) return rb_gsl_eigen_symmv_narray(argc, argv, obj);
#endif
      CHECK_MATRIX(argv[0]);
      Data_Get_Struct(argv[0], gsl_matrix, Atmp);
      if (CLASS_OF(argv[1]) != cgsl_eigen_symmv_workspace)
	rb_raise(rb_eTypeError, 
		 "argv[1]: wrong argument type %s (Eigen::Symmv::Workspace expected)",
		 rb_class2name(CLASS_OF(argv[1])));
      Data_Get_Struct(argv[1], gsl_eigen_symmv_workspace, w);
      break;
    case 1:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv[0])) return rb_gsl_eigen_symmv_narray(argc, argv, obj);
#endif
      CHECK_MATRIX(argv[0]);
      Data_Get_Struct(argv[0], gsl_matrix, Atmp);
      w = gsl_eigen_symmv_alloc(Atmp->size1);
      flagw = 1;
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    }
    break;
  default:
    CHECK_MATRIX(obj);
    Data_Get_Struct(obj, gsl_matrix, Atmp);
    switch (argc) {
    case 1:
      if (CLASS_OF(argv[0]) != cgsl_eigen_symmv_workspace)
	rb_raise(rb_eTypeError, 
		 "argv[0]: wrong argument type %s (Eigen::Symmv::Workspace expected)", 
		 rb_class2name(CLASS_OF(argv[0])));

      Data_Get_Struct(argv[0], gsl_eigen_symmv_workspace, w);
      break;
    case 0:
      w = gsl_eigen_symmv_alloc(Atmp->size1);
      flagw = 1;
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
      break;
    }
  }
  A = make_matrix_clone(Atmp);
  em = gsl_matrix_alloc(A->size1, A->size2);
  v = gsl_vector_alloc(A->size1);
  gsl_eigen_symmv(A, v, em, w);
  /*  gsl_eigen_symmv_sort(v, em, GSL_EIGEN_SORT_VAL_DESC);*/
  gsl_matrix_free(A);
  if (flagw == 1) gsl_eigen_symmv_free(w);
  vval = Data_Wrap_Struct(cgsl_eigen_values, 0, gsl_vector_free, v);
  vvec = Data_Wrap_Struct(cgsl_eigen_vectors, 0, gsl_matrix_free, em);
  return rb_ary_new3(2, vval, vvec);
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_eigen_symmv_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na;
  VALUE eval, evec;
  gsl_matrix *A = NULL;
  gsl_eigen_symmv_workspace *w = NULL;
  gsl_matrix_view mv;
  gsl_vector_view vv;
  int shape1[1], shape2[2];
  int flagw = 0;
  switch (argc) {
  case 2:
    if (!NA_IsNArray(argv[0])) 
      rb_raise(rb_eTypeError, "wrong argument type %s (NArray expected)", 
	       rb_class2name(CLASS_OF(argv[0])));
    GetNArray(argv[0], na);
    if (na->rank < 2) rb_raise(rb_eRuntimeError, "rank >= 2 required");
    if (na->shape[0] != na->shape[1])
      rb_raise(rb_eRuntimeError, "square matrix required");
    A = gsl_matrix_alloc(na->shape[1], na->shape[0]);
    memcpy(A->data, (double*) na->ptr, sizeof(double)*A->size1*A->size2);
    if (CLASS_OF(argv[1]) != cgsl_eigen_symmv_workspace)
      rb_raise(rb_eTypeError, 
	       "argv[1]:  wrong argument type %s (Eigen::Symm::Workspace expected", 
	       rb_class2name(CLASS_OF(argv[1])));
    Data_Get_Struct(argv[1], gsl_eigen_symmv_workspace, w);
    flagw = 0;
    break;
  case 1:
    if (!NA_IsNArray(argv[0])) 
      rb_raise(rb_eTypeError, "wrong argument type %s (NArray expected)", 
	       rb_class2name(CLASS_OF(argv[0])));
    GetNArray(argv[0], na);
    if (na->rank < 2) rb_raise(rb_eRuntimeError, "rank >= 2 required");
    if (na->shape[0] != na->shape[1])
      rb_raise(rb_eRuntimeError, "square matrix required");
    A = gsl_matrix_alloc(na->shape[1], na->shape[0]);
    memcpy(A->data, (double*) na->ptr, sizeof(double)*A->size1*A->size2);
    w = gsl_eigen_symmv_alloc(A->size1);
    flagw = 1;
    break;
  default:
    rb_raise(rb_eArgError, "matrix not given");
    break;
  }
  shape1[0] = A->size1;
  shape2[0] = A->size1;
  shape2[1] = A->size1;
  eval = na_make_object(NA_DFLOAT, 1, shape1, cNVector);
  evec = na_make_object(NA_DFLOAT, 2, shape2, CLASS_OF(argv[0]));
  vv = gsl_vector_view_array(NA_PTR_TYPE(eval,double*), A->size1);
  mv = gsl_matrix_view_array(NA_PTR_TYPE(evec,double*), A->size1, A->size2);
  gsl_eigen_symmv(A, &vv.vector, &mv.matrix, w);
  /*  gsl_sort_vector(v);*/
  gsl_matrix_free(A);
  if (flagw == 1) gsl_eigen_symmv_free(w);
  return rb_ary_new3(2, eval, evec);
}
#endif

static VALUE rb_gsl_eigen_herm(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *Atmp = NULL, *A = NULL;
  gsl_eigen_herm_workspace *w = NULL;
  gsl_vector *v = NULL;
  int flagw = 0;

  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    switch (argc) {
    case 2:
      CHECK_MATRIX_COMPLEX(argv[0]);
      Data_Get_Struct(argv[0], gsl_matrix_complex, Atmp);
      if (CLASS_OF(argv[1]) != cgsl_eigen_herm_workspace)
	rb_raise(rb_eTypeError, 
		 "argv[1]: wrong argument type %s (Eigen::Herm::Workspace expected)", 
		 rb_class2name(CLASS_OF(argv[1])));
      Data_Get_Struct(argv[1], gsl_eigen_herm_workspace, w);
      break;
    case 1:
      CHECK_MATRIX_COMPLEX(argv[0]);
      Data_Get_Struct(argv[0], gsl_matrix_complex, Atmp);
      w = gsl_eigen_herm_alloc(Atmp->size1);
      flagw = 1;
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    }
    break;
  default:
    CHECK_MATRIX_COMPLEX(obj);
    Data_Get_Struct(obj, gsl_matrix_complex, Atmp);
    switch (argc) {
    case 1:
      if (CLASS_OF(argv[0]) != cgsl_eigen_herm_workspace)
	rb_raise(rb_eTypeError, 
		 "argv[0]: wrong argument type %s (Eigen::Herm::Workspace expected)", 
		 rb_class2name(CLASS_OF(argv[0])));

      Data_Get_Struct(argv[0], gsl_eigen_herm_workspace, w);
      break;
    case 0:
      w = gsl_eigen_herm_alloc(Atmp->size1);
      flagw = 1;
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    }
  }
  A = make_matrix_complex_clone(Atmp);
  v = gsl_vector_alloc(A->size1);
  gsl_eigen_herm(A, v, w);
  /*  gsl_sort_vector(v);*/
  gsl_matrix_complex_free(A);
  if (flagw == 1) gsl_eigen_herm_free(w);
  return Data_Wrap_Struct(cgsl_eigen_values, 0, gsl_vector_free, v);
}

static VALUE rb_gsl_eigen_hermv(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix_complex *Atmp = NULL, *A = NULL, *em = NULL;
  gsl_eigen_hermv_workspace *w = NULL;
  gsl_vector *v = NULL;
  int flagw = 0;
  VALUE vval, vvec;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    switch (argc) {
    case 2:
      CHECK_MATRIX_COMPLEX(argv[0]);
      Data_Get_Struct(argv[0], gsl_matrix_complex, Atmp);
      if (CLASS_OF(argv[1]) != cgsl_eigen_hermv_workspace)
	rb_raise(rb_eTypeError, 
		 "argv[1]: wrong argument type %s (Eigen::Hermv::Workspace expected)", 
		 rb_class2name(CLASS_OF(argv[1])));
      Data_Get_Struct(argv[1], gsl_eigen_hermv_workspace, w);
      break;
    case 1:
      CHECK_MATRIX_COMPLEX(argv[0]);
      Data_Get_Struct(argv[0], gsl_matrix_complex, Atmp);
      w = gsl_eigen_hermv_alloc(Atmp->size1);
      flagw = 1;
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    }
    break;
  default:
    CHECK_MATRIX_COMPLEX(obj);
    Data_Get_Struct(obj, gsl_matrix_complex, Atmp);
    switch (argc) {
    case 1:
      if (CLASS_OF(argv[0]) != cgsl_eigen_hermv_workspace)
	rb_raise(rb_eTypeError, 
		 "argv[0]: wrong argument type %s (Eigen::Hermv::Workspace expected)", 
		 rb_class2name(CLASS_OF(argv[0])));

      Data_Get_Struct(argv[0], gsl_eigen_hermv_workspace, w);
      break;
    case 0:
      w = gsl_eigen_hermv_alloc(Atmp->size1);
      flagw = 1;
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    }
  }
  A = make_matrix_complex_clone(Atmp);
  em = gsl_matrix_complex_alloc(A->size1, A->size2);
  v = gsl_vector_alloc(A->size1);
  gsl_eigen_hermv(A, v, em, w);
  /*  gsl_eigen_hermv_sort(v, em, GSL_EIGEN_SORT_VAL_DESC);*/
  gsl_matrix_complex_free(A);
  if (flagw == 1) gsl_eigen_hermv_free(w);
  vval = Data_Wrap_Struct(cgsl_eigen_values, 0, gsl_vector_free, v);
  vvec = Data_Wrap_Struct(cgsl_eigen_herm_vectors, 0, gsl_matrix_complex_free, em);
  return rb_ary_new3(2, vval, vvec);
}

static VALUE rb_gsl_eigen_vectors_unpack(VALUE obj)
{
  gsl_matrix *m = NULL;
  gsl_vector *v = NULL;
  size_t i, j;
  double val;
  VALUE ary, tmp;
  Data_Get_Struct(obj, gsl_matrix, m);
  ary = rb_ary_new2(m->size1);
  for (i = 0; i < m->size1; i++) {
    v = gsl_vector_alloc(m->size2);
    for (j = 0; j < m->size2; j++) {
      val = gsl_matrix_get(m, j, i);
      gsl_vector_set(v, j, val);
    }
    tmp = Data_Wrap_Struct(cgsl_eigen_vector, 0, gsl_vector_free, v);
    rb_ary_store(ary, i, tmp);
  }
  return ary;
}

static VALUE rb_gsl_eigen_vectors_complex_unpack(VALUE obj)
{
  gsl_matrix_complex *m = NULL;
  gsl_vector_complex *v = NULL;
  size_t i, j;
  gsl_complex z;
  VALUE ary, tmp;
  Data_Get_Struct(obj, gsl_matrix_complex, m);
  ary = rb_ary_new2(m->size1);
  for (i = 0; i < m->size1; i++) {
    v = gsl_vector_complex_alloc(m->size2);
    for (j = 0; j < m->size2; j++) {
      z= gsl_matrix_complex_get(m, j, i);
      gsl_vector_complex_set(v, j, z);
    }
    tmp = Data_Wrap_Struct(cgsl_eigen_vector_complex, 0, gsl_vector_complex_free, v);
    rb_ary_store(ary, i, tmp);
  }
  return ary;
}

static void rb_gsl_eigen_define_const(VALUE topmodule, VALUE module)
{
  rb_define_const(topmodule, "EIGEN_SORT_VAL_ASC", INT2FIX(GSL_EIGEN_SORT_VAL_ASC));
  rb_define_const(topmodule, "EIGEN_SORT_VAL_DESC", INT2FIX(GSL_EIGEN_SORT_VAL_DESC));
  rb_define_const(topmodule, "EIGEN_SORT_ABS_ASC", INT2FIX(GSL_EIGEN_SORT_ABS_ASC));
  rb_define_const(topmodule, "EIGEN_SORT_ABS_DESC", INT2FIX(GSL_EIGEN_SORT_ABS_DESC));
    
  rb_define_const(module, "SORT_VAL_ASC", INT2FIX(GSL_EIGEN_SORT_VAL_ASC));
  rb_define_const(module, "SORT_VAL_DESC", INT2FIX(GSL_EIGEN_SORT_VAL_DESC));
  rb_define_const(module, "SORT_ABS_ASC", INT2FIX(GSL_EIGEN_SORT_ABS_ASC));
  rb_define_const(module, "SORT_ABS_DESC", INT2FIX(GSL_EIGEN_SORT_ABS_DESC));

  rb_define_const(module, "VAL_ASC", INT2FIX(GSL_EIGEN_SORT_VAL_ASC));
  rb_define_const(module, "VAL_DESC", INT2FIX(GSL_EIGEN_SORT_VAL_DESC));
  rb_define_const(module, "ABS_ASC", INT2FIX(GSL_EIGEN_SORT_ABS_ASC));
  rb_define_const(module, "ABS_DESC", INT2FIX(GSL_EIGEN_SORT_ABS_DESC));
}

static VALUE rb_gsl_eigen_real_sort(int argc, VALUE *argv, VALUE obj,
  int (*sortfunc)(gsl_vector*, gsl_matrix*, gsl_eigen_sort_t))
{
  gsl_vector *v = NULL;
  gsl_matrix *m = NULL;
  gsl_eigen_sort_t type = GSL_EIGEN_SORT_VAL_DESC;
  switch (argc) {
  case 3:
    CHECK_FIXNUM(argv[2]);
    type = FIX2INT(argv[2]);
    /* no break, do next */
  case 2:
    if (argv[0] == Qnil) {
       v = NULL;
    } else {
      CHECK_VECTOR(argv[0]);
      Data_Get_Struct(argv[0], gsl_vector, v);
    }
    if (argv[1] == Qnil) {
      m = NULL;
    } else {
      CHECK_MATRIX(argv[1]);
      Data_Get_Struct(argv[1], gsl_matrix, m);
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
  }
  return INT2FIX((*sortfunc)(v, m, type));
}


static VALUE rb_gsl_eigen_complex_sort(int argc, VALUE *argv, VALUE obj,
  int (*sortfunc)(gsl_vector*, gsl_matrix_complex*, gsl_eigen_sort_t))
{
  gsl_vector *v = NULL;
  gsl_matrix_complex *m = NULL;
  gsl_eigen_sort_t type = GSL_EIGEN_SORT_VAL_DESC;

  switch (argc) {
  case 3:
    CHECK_FIXNUM(argv[2]);
    type = FIX2INT(argv[2]);
    /* no break, do next */
  case 2:
    if (argv[0] == Qnil) {
      v = NULL;
    } else {
      CHECK_VECTOR(argv[0]);
      Data_Get_Struct(argv[0], gsl_vector, v);
    }
    if (argv[1] == Qnil) {
      m = NULL;
    } else {
      CHECK_MATRIX_COMPLEX(argv[1]);
      Data_Get_Struct(argv[1], gsl_matrix_complex, m);        
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
  }
  return INT2FIX((*sortfunc)(v, m, type));
}
static VALUE rb_gsl_eigen_symmv_sort(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_eigen_real_sort(argc, argv, obj, gsl_eigen_symmv_sort);
}
static VALUE rb_gsl_eigen_hermv_sort(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_eigen_complex_sort(argc, argv, obj, gsl_eigen_hermv_sort);  
}


#ifdef HAVE_EIGEN_FRANCIS
static VALUE rb_gsl_eigen_francis_alloc(VALUE klass)
{
  gsl_eigen_francis_workspace *w = NULL;
  w = gsl_eigen_francis_alloc();
  return Data_Wrap_Struct(klass, 0, gsl_eigen_francis_free, w);
}

static VALUE rb_gsl_eigen_francis_T(int argc, VALUE *argv, VALUE obj)
{
  gsl_eigen_francis_workspace *w = NULL;
  int istart = 0;
  if (CLASS_OF(obj) == cgsl_eigen_francis_workspace) {
    Data_Get_Struct(obj, gsl_eigen_francis_workspace, w);
    istart = 0;
  } else {
    if (argc != 1) rb_raise(rb_eArgError, "too few arguments (%d for 1)\n", argc);
    Data_Get_Struct(argv[0], gsl_eigen_francis_workspace, w);
    istart = 1;
  }
  gsl_eigen_francis_T(FIX2INT(argv[istart]), w);

  return Qtrue;
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_eigen_francis_narray(int argc, VALUE *argv, VALUE obj);
#endif

static VALUE rb_gsl_eigen_francis(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m, *mtmp;
  gsl_vector_complex *v;
  gsl_eigen_francis_workspace *w;
  int vflag = 0, wflag = 0;
  int istart = 0;
  VALUE *argv2;

#ifdef HAVE_NARRAY_H
  if (NA_IsNArray(obj)) return rb_gsl_eigen_francis_narray(argc, argv, obj);
  if (argc >= 1 && NA_IsNArray(argv[0])) 
    return rb_gsl_eigen_francis_narray(argc, argv, obj);
#endif

  if (MATRIX_P(obj)) {
    Data_Get_Struct(obj, gsl_matrix, m);
    argv2 = argv;
    istart = 0;
  } else {
    if (argc < 1) rb_raise(rb_eArgError, "Wrong number of arguments.\n");
    Data_Get_Struct(argv[0], gsl_matrix, m);
    istart = 1;
    argv2 = argv + 1;
  }

  switch (argc-istart) {
  case 0:
    v = gsl_vector_complex_alloc(m->size1);
    w = gsl_eigen_francis_alloc();
    vflag = 1;
    wflag = 1;
    break;
  case 1:
    if (CLASS_OF(argv2[0]) == cgsl_vector_complex) {
      Data_Get_Struct(argv2[0], gsl_vector_complex, v);
      w = gsl_eigen_francis_alloc();
      wflag = 1;      
    } else if (CLASS_OF(argv2[0]) == cgsl_eigen_francis_workspace) {
      v = gsl_vector_complex_alloc(m->size1);
      vflag = 1;
      Data_Get_Struct(argv2[0], gsl_eigen_francis_workspace, w);
    } else {
      rb_raise(rb_eArgError, "Wrong argument type.\n");
    }
    break;
  case 2:
    CHECK_VECTOR_COMPLEX(argv2[0]);
    if (CLASS_OF(argv2[1]) != cgsl_eigen_francis_workspace) {
      rb_raise(rb_eArgError, "argv[1] must be a GSL::Eigen::francis::Workspace.\n");
    }
    Data_Get_Struct(argv2[0], gsl_vector_complex, v);
    Data_Get_Struct(argv2[1], gsl_eigen_francis_workspace, w);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-2).\n", argc);
  }
  mtmp = make_matrix_clone(m);
  gsl_eigen_francis(mtmp, v, w);
  gsl_matrix_free(mtmp);
  if (wflag == 1) gsl_eigen_francis_free(w);
  if (vflag == 1)
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, v);
  else
    return argv2[0];
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_eigen_francis_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na;
  VALUE nary;
  gsl_matrix *A = NULL;
  gsl_eigen_francis_workspace *w = NULL;
  gsl_vector_complex_view vv;
  int shape[1];
  int flagw = 0;
  switch (argc) {
  case 2:
    if (!NA_IsNArray(argv[0])) 
      rb_raise(rb_eTypeError, "wrong argument type %s (NArray expected)", 
	       rb_class2name(CLASS_OF(argv[0])));
    GetNArray(argv[0], na);
    if (na->rank < 2) rb_raise(rb_eRuntimeError, "rank >= 2 required");
    if (na->shape[0] != na->shape[1])
      rb_raise(rb_eRuntimeError, "square matrix required");
    A = gsl_matrix_alloc(na->shape[1], na->shape[0]);
    memcpy(A->data, (double*) na->ptr, sizeof(double)*A->size1*A->size2);
    if (CLASS_OF(argv[1]) != cgsl_eigen_francis_workspace)
      rb_raise(rb_eTypeError, 
	       "argv[1]:  wrong argument type %s (Eigen::Symm::Workspace expected", 
	       rb_class2name(CLASS_OF(argv[1])));
    Data_Get_Struct(argv[1], gsl_eigen_francis_workspace, w);
    flagw = 0;
    break;
  case 1:
    if (!NA_IsNArray(argv[0])) 
      rb_raise(rb_eTypeError, "wrong argument type %s (NArray expected)", 
	       rb_class2name(CLASS_OF(argv[0])));
    GetNArray(argv[0], na);
    if (na->rank < 2) rb_raise(rb_eRuntimeError, "rank >= 2 required");
    if (na->shape[0] != na->shape[1])
      rb_raise(rb_eRuntimeError, "square matrix required");
    A = gsl_matrix_alloc(na->shape[1], na->shape[0]);
    memcpy(A->data, (double*) na->ptr, sizeof(double)*A->size1*A->size2);
    w = gsl_eigen_francis_alloc();
    flagw = 1;
    break;
  default:
    rb_raise(rb_eArgError, "matrix not given");
    break;
  }
  shape[0] = A->size1;
  nary = na_make_object(NA_DCOMPLEX, 1, shape, cNVector);
  vv = gsl_vector_complex_view_array(NA_PTR_TYPE(nary,double*), A->size1);
  gsl_eigen_francis(A, &vv.vector, w);
  /*  gsl_sort_vector(v);*/
  gsl_matrix_free(A);
  if (flagw == 1) gsl_eigen_francis_free(w);
  return nary;
}
#endif

static VALUE rb_gsl_eigen_francis_Z(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m, *mtmp, *Z;
  gsl_vector_complex *v;
  gsl_eigen_francis_workspace *w;
  int vflag = 0, wflag = 0;
  int istart = 0;
  VALUE *argv2, vv, ZZ;

  if (MATRIX_P(obj)) {
    Data_Get_Struct(obj, gsl_matrix, m);
    argv2 = argv;
    istart = 0;
  } else {
    if (argc < 1) rb_raise(rb_eArgError, "Wrong number of arguments.\n");
    Data_Get_Struct(argv[0], gsl_matrix, m);
    istart = 1;
    argv2 = argv + 1;
  }

  switch (argc-istart) {
  case 0:
    v = gsl_vector_complex_alloc(m->size1);
    Z = gsl_matrix_alloc(m->size1, m->size2);
    w = gsl_eigen_francis_alloc();
    vflag = 1;
    wflag = 1;
    break;
  case 1:
    if (CLASS_OF(argv2[0]) == cgsl_eigen_francis_workspace) {
      v = gsl_vector_complex_alloc(m->size1);
      Z = gsl_matrix_alloc(m->size1, m->size2);
      vflag = 1;
      Data_Get_Struct(argv2[0], gsl_eigen_francis_workspace, w);
    } else {
      rb_raise(rb_eArgError, "Wrong argument type.\n");
    }
    break;
  case 3:
    CHECK_VECTOR_COMPLEX(argv2[0]);
    CHECK_MATRIX(argv2[1]);
    if (CLASS_OF(argv2[2]) != cgsl_eigen_francis_workspace) {
      rb_raise(rb_eArgError, "argv[1] must be a GSL::Eigen::francis::Workspace.\n");
    }
    Data_Get_Struct(argv2[0], gsl_vector_complex, v);
    Data_Get_Struct(argv2[1], gsl_matrix, Z);
    Data_Get_Struct(argv2[2], gsl_eigen_francis_workspace, w);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-2).\n", argc);
  }
  mtmp = make_matrix_clone(m);
  gsl_eigen_francis_Z(mtmp, v, Z, w);
  gsl_matrix_free(mtmp);

  if (wflag == 1) gsl_eigen_francis_free(w);
  if (vflag == 1) {
    vv = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, v);
    ZZ = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Z);
  } else {
    vv = argv2[0];
    ZZ = argv2[1];
  }
  return rb_ary_new3(2, vv, ZZ);
}
#endif

#ifdef GSL_1_9_LATER
static VALUE rb_gsl_eigen_nonsymm_alloc(VALUE klass, VALUE nn)
{
  size_t n;
  gsl_eigen_nonsymm_workspace *w = NULL;
  n = (size_t) FIX2UINT(nn);
  w = gsl_eigen_nonsymm_alloc(n);
  return Data_Wrap_Struct(cgsl_eigen_nonsymm_workspace, 0, gsl_eigen_nonsymm_free, w);
}

static VALUE rb_gsl_eigen_nonsymm_params(int argc, VALUE *argv, VALUE obj)
{
  gsl_eigen_nonsymm_workspace *w = NULL;
  int istart = 0;
  if (CLASS_OF(obj) == cgsl_eigen_nonsymm_workspace) {
    Data_Get_Struct(obj, gsl_eigen_nonsymm_workspace, w);
    istart = 0;
  } else {
    if (argc != 3) rb_raise(rb_eArgError, "too few arguments (%d for 3)\n", argc);
    Data_Get_Struct(argv[2], gsl_eigen_nonsymm_workspace, w);
    istart = 1;
  }
  switch (argc - istart) {
  case 2:
    gsl_eigen_nonsymm_params(FIX2INT(argv[0]), FIX2INT(argv[1]), w);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments.\n");
  }
  return Qtrue;
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_eigen_nonsymm_narray(int argc, VALUE *argv, VALUE obj);
#endif

static VALUE rb_gsl_eigen_nonsymm(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m;
  gsl_vector_complex *v;
  gsl_eigen_nonsymm_workspace *w;
  int vflag = 0, wflag = 0;
  int istart = 0;
  VALUE *argv2;

#ifdef HAVE_NARRAY_H
  if (NA_IsNArray(obj)) return rb_gsl_eigen_nonsymm_narray(argc, argv, obj);
  if (argc >= 1 && NA_IsNArray(argv[0])) 
    return rb_gsl_eigen_nonsymm_narray(argc, argv, obj);
#endif

  if (MATRIX_P(obj)) {
    Data_Get_Struct(obj, gsl_matrix, m);
    argv2 = argv;
    istart = 0;
  } else {
    if (argc < 1) rb_raise(rb_eArgError, "Wrong number of arguments.\n");
    Data_Get_Struct(argv[0], gsl_matrix, m);
    istart = 1;
    argv2 = argv + 1;
  }

  switch (argc-istart) {
  case 0:
    v = gsl_vector_complex_alloc(m->size1);
    w = gsl_eigen_nonsymm_alloc(m->size1);
    vflag = 1;
    wflag = 1;
    break;
  case 1:
    if (CLASS_OF(argv2[0]) == cgsl_vector_complex) {
      Data_Get_Struct(argv2[0], gsl_vector_complex, v);
      w = gsl_eigen_nonsymm_alloc(m->size1);
      wflag = 1;      
    } else if (CLASS_OF(argv2[0]) == cgsl_eigen_nonsymm_workspace) {
      v = gsl_vector_complex_alloc(m->size1);
      vflag = 1;
      Data_Get_Struct(argv2[0], gsl_eigen_nonsymm_workspace, w);
    } else {
      rb_raise(rb_eArgError, "Wrong argument type.\n");
    }
    break;
  case 2:
    CHECK_VECTOR_COMPLEX(argv2[0]);
    if (CLASS_OF(argv2[1]) != cgsl_eigen_nonsymm_workspace) {
      rb_raise(rb_eArgError, "argv[1] must be a GSL::Eigen::Nonsymm::Workspace.\n");
    }
    Data_Get_Struct(argv2[0], gsl_vector_complex, v);
    Data_Get_Struct(argv2[1], gsl_eigen_nonsymm_workspace, w);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-2).\n", argc);
  }
//  mtmp = make_matrix_clone(m);
  gsl_eigen_nonsymm(m, v, w);
//  gsl_matrix_free(mtmp);
  if (wflag == 1) gsl_eigen_nonsymm_free(w);
  if (vflag == 1)
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, v);
  else
    return argv2[0];
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_eigen_nonsymm_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na;
  VALUE nary;
  gsl_matrix *A = NULL;
  gsl_eigen_nonsymm_workspace *w = NULL;
  gsl_vector_complex_view vv;
  int shape[1];
  int flagw = 0;
  switch (argc) {
  case 2:
    if (!NA_IsNArray(argv[0])) 
      rb_raise(rb_eTypeError, "wrong argument type %s (NArray expected)", 
	       rb_class2name(CLASS_OF(argv[0])));
    GetNArray(argv[0], na);
    if (na->rank < 2) rb_raise(rb_eRuntimeError, "rank >= 2 required");
    if (na->shape[0] != na->shape[1])
      rb_raise(rb_eRuntimeError, "square matrix required");
    A = gsl_matrix_alloc(na->shape[1], na->shape[0]);
    memcpy(A->data, (double*) na->ptr, sizeof(double)*A->size1*A->size2);
    if (CLASS_OF(argv[1]) != cgsl_eigen_nonsymm_workspace)
      rb_raise(rb_eTypeError, 
	       "argv[1]:  wrong argument type %s (Eigen::Symm::Workspace expected", 
	       rb_class2name(CLASS_OF(argv[1])));
    Data_Get_Struct(argv[1], gsl_eigen_nonsymm_workspace, w);
    flagw = 0;
    break;
  case 1:
    if (!NA_IsNArray(argv[0])) 
      rb_raise(rb_eTypeError, "wrong argument type %s (NArray expected)", 
	       rb_class2name(CLASS_OF(argv[0])));
    GetNArray(argv[0], na);
    if (na->rank < 2) rb_raise(rb_eRuntimeError, "rank >= 2 required");
    if (na->shape[0] != na->shape[1])
      rb_raise(rb_eRuntimeError, "square matrix required");
    A = gsl_matrix_alloc(na->shape[1], na->shape[0]);
    memcpy(A->data, (double*) na->ptr, sizeof(double)*A->size1*A->size2);
    w = gsl_eigen_nonsymm_alloc(A->size1);
    flagw = 1;
    break;
  default:
    rb_raise(rb_eArgError, "matrix not given");
    break;
  }
  shape[0] = A->size1;
  nary = na_make_object(NA_DCOMPLEX, 1, shape, cNVector);
  vv = gsl_vector_complex_view_array(NA_PTR_TYPE(nary,double*), A->size1);
  gsl_eigen_nonsymm(A, &vv.vector, w);
  /*  gsl_sort_vector(v);*/
  gsl_matrix_free(A);
  if (flagw == 1) gsl_eigen_nonsymm_free(w);
  return nary;
}
#endif

static VALUE rb_gsl_eigen_nonsymm_Z(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m, *Z;
  gsl_vector_complex *v;
  gsl_eigen_nonsymm_workspace *w;
  int vflag = 0, wflag = 0;
  int istart = 0;
  VALUE *argv2, vv, ZZ;

  if (MATRIX_P(obj)) {
    Data_Get_Struct(obj, gsl_matrix, m);
    argv2 = argv;
    istart = 0;
  } else {
    if (argc < 1) rb_raise(rb_eArgError, "Wrong number of arguments.\n");
    Data_Get_Struct(argv[0], gsl_matrix, m);
    istart = 1;
    argv2 = argv + 1;
  }

  switch (argc-istart) {
  case 0:
    v = gsl_vector_complex_alloc(m->size1);
    Z = gsl_matrix_alloc(m->size1, m->size2);
    w = gsl_eigen_nonsymm_alloc(m->size1);
    vflag = 1;
    wflag = 1;
    break;
  case 1:
    if (CLASS_OF(argv2[0]) == cgsl_eigen_nonsymm_workspace) {
      v = gsl_vector_complex_alloc(m->size1);
      Z = gsl_matrix_alloc(m->size1, m->size2);
      vflag = 1;
      Data_Get_Struct(argv2[0], gsl_eigen_nonsymm_workspace, w);
    } else {
      rb_raise(rb_eArgError, "Wrong argument type.\n");
    }
    break;
  case 3:
    CHECK_VECTOR_COMPLEX(argv2[0]);
    CHECK_MATRIX(argv2[1]);
    if (CLASS_OF(argv2[2]) != cgsl_eigen_nonsymm_workspace) {
      rb_raise(rb_eArgError, "argv[1] must be a GSL::Eigen::Nonsymm::Workspace.\n");
    }
    Data_Get_Struct(argv2[0], gsl_vector_complex, v);
    Data_Get_Struct(argv2[1], gsl_matrix, Z);
    Data_Get_Struct(argv2[2], gsl_eigen_nonsymm_workspace, w);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-2).\n", argc);
  }
//  mtmp = make_matrix_clone(m);
  gsl_eigen_nonsymm_Z(m, v, Z, w);
//  gsl_matrix_free(mtmp);

  if (wflag == 1) gsl_eigen_nonsymm_free(w);
  if (vflag == 1) {
    vv = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, v);
    ZZ = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Z);
  } else {
    vv = argv2[0];
    ZZ = argv2[1];
  }
  return rb_ary_new3(2, vv, ZZ);
}

static VALUE rb_gsl_eigen_nonsymmv_alloc(VALUE klass, VALUE nn)
{
  size_t n;
  gsl_eigen_nonsymmv_workspace *w = NULL;
  n = (size_t) FIX2UINT(nn);
  w = gsl_eigen_nonsymmv_alloc(n);
  return Data_Wrap_Struct(cgsl_eigen_nonsymmv_workspace, 0, gsl_eigen_nonsymmv_free, w);
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_eigen_nonsymmv_narray(int argc, VALUE *argv, VALUE obj);
#endif

static VALUE rb_gsl_eigen_nonsymmv(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m;
  gsl_vector_complex *v = NULL;
  gsl_matrix_complex *evec = NULL;
  gsl_eigen_nonsymmv_workspace *w;
  int vflag = 0, wflag = 0;
  int istart = 0;
  VALUE *argv2;

#ifdef HAVE_NARRAY_H
  if (NA_IsNArray(obj)) return rb_gsl_eigen_nonsymmv_narray(argc, argv, obj);
  if (argc >= 1 && NA_IsNArray(argv[0])) 
    return rb_gsl_eigen_nonsymmv_narray(argc, argv, obj);
#endif

  if (MATRIX_P(obj)) {
    Data_Get_Struct(obj, gsl_matrix, m);
    argv2 = argv;
    istart = 0;
  } else {
    if (argc < 1) rb_raise(rb_eArgError, "Wrong number of arguments.\n");
    Data_Get_Struct(argv[0], gsl_matrix, m);
    istart = 1;
    argv2 = argv + 1;
  }

  switch (argc-istart) {
  case 0:
    v = gsl_vector_complex_alloc(m->size1);
    evec = gsl_matrix_complex_alloc(m->size1, m->size2);
    w = gsl_eigen_nonsymmv_alloc(m->size1);
    vflag = 1;
    wflag = 1;
    break;
  case 1:
    if (CLASS_OF(argv2[0]) == cgsl_eigen_nonsymmv_workspace) {
      v = gsl_vector_complex_alloc(m->size1);
      evec = gsl_matrix_complex_alloc(m->size1, m->size2);
      vflag = 1;
      Data_Get_Struct(argv2[0], gsl_eigen_nonsymmv_workspace, w);
    } else {
      rb_raise(rb_eArgError, "Wrong argument type.\n");
    }
    break;
  case 2:
    CHECK_VECTOR_COMPLEX(argv2[0]);
    CHECK_MATRIX_COMPLEX(argv2[1]);
    w = gsl_eigen_nonsymmv_alloc(m->size1);
    wflag = 1;
    break;
  case 3:
    CHECK_VECTOR_COMPLEX(argv2[0]);
    CHECK_MATRIX_COMPLEX(argv2[1]);
    if (CLASS_OF(argv2[2]) != cgsl_eigen_nonsymmv_workspace) {
      rb_raise(rb_eArgError, "argv[1] must be a GSL::Eigen::Nonsymm::Workspace.\n");
    }
    Data_Get_Struct(argv2[0], gsl_vector_complex, v);
    Data_Get_Struct(argv2[1], gsl_matrix_complex, evec);
    Data_Get_Struct(argv2[2], gsl_eigen_nonsymmv_workspace, w);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-3).\n", argc);
  }
//  mtmp = make_matrix_clone(m);
  gsl_eigen_nonsymmv(m, v, evec, w);
//  gsl_matrix_free(mtmp);

  if (wflag == 1) gsl_eigen_nonsymmv_free(w);
  if (vflag == 1) {
    return rb_ary_new3(2, 
		       Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, v),
		       Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, evec));
  }  else {
    return rb_ary_new3(2, argv2[0], argv2[1]);
  }
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_eigen_nonsymmv_narray(int argc, VALUE *argv, VALUE obj)
{
  struct NARRAY *na;
  VALUE nary, nvec;
  gsl_matrix *A = NULL;
  gsl_eigen_nonsymmv_workspace *w = NULL;
  gsl_vector_complex_view vv;
  gsl_matrix_complex_view mm;
  int shape[1], shape2[2];
  int flagw = 0;
  switch (argc) {
  case 2:
    if (!NA_IsNArray(argv[0])) 
      rb_raise(rb_eTypeError, "wrong argument type %s (NArray expected)", 
	       rb_class2name(CLASS_OF(argv[0])));
    GetNArray(argv[0], na);
    if (na->rank < 2) rb_raise(rb_eRuntimeError, "rank >= 2 required");
    if (na->shape[0] != na->shape[1])
      rb_raise(rb_eRuntimeError, "square matrix required");
    A = gsl_matrix_alloc(na->shape[1], na->shape[0]);
    memcpy(A->data, (double*) na->ptr, sizeof(double)*A->size1*A->size2);
    if (CLASS_OF(argv[1]) != cgsl_eigen_nonsymmv_workspace)
      rb_raise(rb_eTypeError, 
	       "argv[1]:  wrong argument type %s (Eigen::Symm::Workspace expected", 
	       rb_class2name(CLASS_OF(argv[1])));
    Data_Get_Struct(argv[1], gsl_eigen_nonsymmv_workspace, w);
    flagw = 0;
    break;
  case 1:
    if (!NA_IsNArray(argv[0])) 
      rb_raise(rb_eTypeError, "wrong argument type %s (NArray expected)", 
	       rb_class2name(CLASS_OF(argv[0])));
    GetNArray(argv[0], na);
    if (na->rank < 2) rb_raise(rb_eRuntimeError, "rank >= 2 required");
    if (na->shape[0] != na->shape[1])
      rb_raise(rb_eRuntimeError, "square matrix required");
    A = gsl_matrix_alloc(na->shape[1], na->shape[0]);
    memcpy(A->data, (double*) na->ptr, sizeof(double)*A->size1*A->size2);
    w = gsl_eigen_nonsymmv_alloc(A->size1);
    flagw = 1;
    break;
  default:
    rb_raise(rb_eArgError, "matrix not given");
    break;
  }
  shape[0] = A->size1; shape2[0] = A->size1; shape2[1] = A->size2;
  nary = na_make_object(NA_DCOMPLEX, 1, shape, cNVector);
  vv = gsl_vector_complex_view_array(NA_PTR_TYPE(nary,double*), A->size1);
  nvec = na_make_object(NA_DCOMPLEX, 2, shape2, CLASS_OF(argv[0]));
  mm = gsl_matrix_complex_view_array(NA_PTR_TYPE(nvec,double*), A->size1, A->size2);
  gsl_eigen_nonsymmv(A, &vv.vector, &mm.matrix, w);
  /*  gsl_sort_vector(v);*/
  gsl_matrix_free(A);
  if (flagw == 1) gsl_eigen_nonsymmv_free(w);
  return rb_ary_new3(2, nary, nvec);
}
#endif

static VALUE rb_gsl_eigen_nonsymmv_Z(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m, *Z = NULL;
  gsl_vector_complex *v = NULL;
  gsl_matrix_complex *evec = NULL;
  gsl_eigen_nonsymmv_workspace *w;
  int vflag = 0, wflag = 0;
  int istart = 0;
  VALUE *argv2;

  if (MATRIX_P(obj)) {
    Data_Get_Struct(obj, gsl_matrix, m);
    argv2 = argv;
    istart = 0;
  } else {
    if (argc < 1) rb_raise(rb_eArgError, "Wrong number of arguments.\n");
    Data_Get_Struct(argv[0], gsl_matrix, m);
    istart = 1;
    argv2 = argv + 1;
  }

  switch (argc-istart) {
  case 0:
    v = gsl_vector_complex_alloc(m->size1);
    evec = gsl_matrix_complex_alloc(m->size1, m->size2);
    Z = gsl_matrix_alloc(m->size1, m->size2);
    w = gsl_eigen_nonsymmv_alloc(m->size1);
    vflag = 1;
    wflag = 1;
    break;
  case 1:
    if (CLASS_OF(argv2[0]) == cgsl_eigen_nonsymm_workspace) {
      v = gsl_vector_complex_alloc(m->size1);
      evec = gsl_matrix_complex_alloc(m->size1, m->size2);
      vflag = 1;
      Data_Get_Struct(argv2[0], gsl_eigen_nonsymmv_workspace, w);
    } else {
      rb_raise(rb_eArgError, "Wrong argument type.\n");
    }
    break;
  case 3:
    CHECK_VECTOR_COMPLEX(argv2[0]);
    CHECK_MATRIX_COMPLEX(argv2[1]);
    CHECK_MATRIX(argv2[2]);
    w = gsl_eigen_nonsymmv_alloc(m->size1);
    wflag = 1;
    break;
  case 4:
    CHECK_VECTOR_COMPLEX(argv2[0]);
    CHECK_MATRIX_COMPLEX(argv2[1]);
    CHECK_MATRIX(argv2[2]);
    if (CLASS_OF(argv2[3]) != cgsl_eigen_nonsymm_workspace) {
      rb_raise(rb_eArgError, "argv[1] must be a GSL::Eigen::Nonsymm::Workspace.\n");
    }
    Data_Get_Struct(argv2[0], gsl_vector_complex, v);
    Data_Get_Struct(argv2[1], gsl_matrix_complex, evec);
    Data_Get_Struct(argv2[1], gsl_matrix, Z);
    Data_Get_Struct(argv2[3], gsl_eigen_nonsymmv_workspace, w);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-3).\n", argc);
  }
//  mtmp = make_matrix_clone(m);
  gsl_eigen_nonsymmv_Z(m, v, evec, Z, w);
//  gsl_matrix_free(mtmp);

  if (wflag == 1) gsl_eigen_nonsymmv_free(w);
  if (vflag == 1) {
    return rb_ary_new3(3, 
		       Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, v),
		       Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, evec),
Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Z));
  }  else {
    return rb_ary_new3(2, argv2[0], argv2[1], argv2[2]);
  }
}

static VALUE rb_gsl_eigen_complex_sort2(int argc, VALUE *argv, VALUE obj,
  int (*sortfunc)(gsl_vector_complex*, gsl_matrix_complex*, gsl_eigen_sort_t))
{
  gsl_vector_complex *v = NULL;
  gsl_matrix_complex *m = NULL;
  gsl_eigen_sort_t type = GSL_EIGEN_SORT_ABS_DESC;

  switch (argc) {
  case 3:
    CHECK_FIXNUM(argv[2]);
    type = FIX2INT(argv[2]);
    /* no break, do next */
  case 2:
    if (argv[0] == Qnil) {
      v = NULL;
    } else {
      CHECK_VECTOR_COMPLEX(argv[0]);
      Data_Get_Struct(argv[0], gsl_vector_complex, v);      
    }
    if (argv[1] == Qnil) {
      m = NULL;
    } else {
      CHECK_MATRIX_COMPLEX(argv[1]);
      Data_Get_Struct(argv[1], gsl_matrix_complex, m);  
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
  }
  return INT2FIX((*sortfunc)(v, m, type));
}


static VALUE rb_gsl_eigen_nonsymmv_sort(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_eigen_complex_sort2(argc, argv, obj, gsl_eigen_nonsymmv_sort);  
}

#endif

#ifdef GSL_1_10_LATER

static VALUE rb_gsl_eigen_gensymm_alloc(VALUE klass, VALUE nn)
{
  gsl_eigen_gensymm_workspace *w = NULL;
  CHECK_FIXNUM(nn);
  w = gsl_eigen_gensymm_alloc(FIX2INT(nn));
  return Data_Wrap_Struct(cgensymm, 0, gsl_eigen_gensymm_free, w);
}
static VALUE rb_gsl_eigen_gensymmv_alloc(VALUE klass, VALUE nn)
{
  gsl_eigen_gensymmv_workspace *w = NULL;
  CHECK_FIXNUM(nn);
  w = gsl_eigen_gensymmv_alloc(FIX2INT(nn));
  return Data_Wrap_Struct(cgensymmv, 0, gsl_eigen_gensymmv_free, w);
}
static VALUE rb_gsl_eigen_genherm_alloc(VALUE klass, VALUE nn)
{
  gsl_eigen_genherm_workspace *w = NULL;
  CHECK_FIXNUM(nn);
  w = gsl_eigen_genherm_alloc(FIX2INT(nn));
  return Data_Wrap_Struct(cgenherm, 0, gsl_eigen_genherm_free, w);
}
static VALUE rb_gsl_eigen_genhermv_alloc(VALUE klass, VALUE nn)
{
  gsl_eigen_genhermv_workspace *w = NULL;
  CHECK_FIXNUM(nn);
  w = gsl_eigen_genhermv_alloc(FIX2INT(nn));
  return Data_Wrap_Struct(cgenhermv, 0, gsl_eigen_genhermv_free, w);
}

static int check_argv_gensymm(int argc, VALUE *argv, VALUE obj, gsl_matrix **A, gsl_matrix **B,
	gsl_vector **eval, gsl_eigen_gensymm_workspace **w)
{
	int argc2 = argc;
	int flag = 0;
	if (CLASS_OF(obj) == cgensymm) {
		Data_Get_Struct(obj, gsl_eigen_gensymm_workspace, *w);
	} else {
		if (rb_obj_is_kind_of(argv[argc-1], cgensymm)) {
			Data_Get_Struct(argv[argc-1], gsl_eigen_gensymm_workspace, *w);
			argc2 = argc-1;
		} else {
			/* workspace is not given */
		}
	}
	switch (argc2) {
	case 2:
		CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix, *A);
		Data_Get_Struct(argv[1], gsl_matrix, *B);
		break;
	case 3:
		if (rb_obj_is_kind_of(argv[2], cgensymm)) {
			Data_Get_Struct(argv[2], gsl_eigen_gensymm_workspace, *w);
		} else {
			CHECK_VECTOR(argv[2]);
			Data_Get_Struct(argv[2], gsl_vector, *eval);
		}
		CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix, *A);
		Data_Get_Struct(argv[1], gsl_matrix, *B);
		break;
	default:
		rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
	}
	if (*eval == NULL) {
		*eval = gsl_vector_alloc((*A)->size1);
		flag += 1;
	}
	if (*w == NULL) {
		*w = gsl_eigen_gensymm_alloc((*A)->size1);
		flag += 2;
	}
	return flag;
}

static int check_argv_gensymmv(int argc, VALUE *argv, VALUE obj, gsl_matrix **A, gsl_matrix **B,
	gsl_vector **eval, gsl_matrix **evec, gsl_eigen_gensymmv_workspace **w)
{
	int argc2 = argc;
	int flag = 0;
	if (CLASS_OF(obj) == cgensymmv) {
		Data_Get_Struct(obj, gsl_eigen_gensymmv_workspace, *w);
	} else {
		if (rb_obj_is_kind_of(argv[argc-1], cgensymmv)) {
			Data_Get_Struct(argv[argc-1], gsl_eigen_gensymmv_workspace, *w);
			argc2 = argc-1;
		} else {

			/* workspace is not given */
		}
	}

	switch (argc2) {
	case 2:
		CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix, *A);
		Data_Get_Struct(argv[1], gsl_matrix, *B);
		break;
	case 3:
		if (rb_obj_is_kind_of(argv[2], cgensymmv)) {
			Data_Get_Struct(argv[2], gsl_eigen_gensymmv_workspace, *w);
		} else {
      rb_raise(rb_eTypeError, "Wrong argument type %s (GSL::Eigen::Gensymmv::Workspace expected)",
        rb_class2name(CLASS_OF(argv[2])));
		}
		CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix, *A);
		Data_Get_Struct(argv[1], gsl_matrix, *B);

		break;
	case 5:
		if (rb_obj_is_kind_of(argv[4], cgensymmv)) {
			Data_Get_Struct(argv[4], gsl_eigen_gensymmv_workspace, *w);
		} else {
      rb_raise(rb_eTypeError, "Wrong argument type %s (GSL::Eigen::Gensymmv::Workspace expected)",
        rb_class2name(CLASS_OF(argv[4])));      
    }
		CHECK_VECTOR(argv[2]);
		Data_Get_Struct(argv[2], gsl_vector, *eval);
		CHECK_MATRIX(argv[3]);
		Data_Get_Struct(argv[3], gsl_matrix, *evec);			
		
		CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix, *A);
		Data_Get_Struct(argv[1], gsl_matrix, *B);
		break;		
	default:
		rb_raise(rb_eArgError, "wrong number of arguments (%d for 2, 3 or 5)", argc);
  }
	if (*eval == NULL && *evec == NULL) {
		*eval = gsl_vector_alloc((*A)->size1);
		*evec = gsl_matrix_alloc((*A)->size1, (*A)->size2);		
		flag += 1;
	}
	if (*w == NULL) {
		*w = gsl_eigen_gensymmv_alloc((*A)->size1);
		flag += 2;
	}
	return flag;
}

static int check_argv_genherm(int argc, VALUE *argv, VALUE obj, gsl_matrix_complex **A, gsl_matrix_complex **B,
	gsl_vector **eval, gsl_eigen_genherm_workspace **w)
{
	int argc2 = argc;
	int flag = 0;
	if (CLASS_OF(obj) == cgenherm) {
		Data_Get_Struct(obj, gsl_eigen_genherm_workspace, *w);
	} else {
		if (rb_obj_is_kind_of(argv[argc-1], cgenherm)) {
			Data_Get_Struct(argv[argc-1], gsl_eigen_genherm_workspace, *w);
			argc2 = argc-1;
		} else {
			/* workspace is not given */
		}
	}
	switch (argc2) {
	case 2:
		CHECK_MATRIX_COMPLEX(argv[0]); CHECK_MATRIX_COMPLEX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix_complex, *A);
		Data_Get_Struct(argv[1], gsl_matrix_complex, *B);
		break;
	case 3:
		if (rb_obj_is_kind_of(argv[2], cgenherm)) {
			Data_Get_Struct(argv[2], gsl_eigen_genherm_workspace, *w);
		} else {
			CHECK_VECTOR(argv[2]);
			Data_Get_Struct(argv[2], gsl_vector, *eval);
		}
		CHECK_MATRIX_COMPLEX(argv[0]); CHECK_MATRIX_COMPLEX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix_complex, *A);
		Data_Get_Struct(argv[1], gsl_matrix_complex, *B);
		break;
	default:
		rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
	}
	if (*eval == NULL) {
		*eval = gsl_vector_alloc((*A)->size1);
		flag += 1;
	}
	if (*w == NULL) {
		*w = gsl_eigen_genherm_alloc((*A)->size1);
		flag += 2;
	}
	return flag;
}
static int check_argv_genhermv(int argc, VALUE *argv, VALUE obj, gsl_matrix_complex **A, gsl_matrix_complex **B,
	gsl_vector **eval, gsl_matrix_complex **evec, gsl_eigen_genhermv_workspace **w)
{
	int argc2 = argc;
	int flag = 0;
	if (CLASS_OF(obj) == cgenhermv) {
		Data_Get_Struct(obj, gsl_eigen_genhermv_workspace, *w);
	} else {
		if (rb_obj_is_kind_of(argv[argc-1], cgenhermv)) {
			Data_Get_Struct(argv[argc-1], gsl_eigen_genhermv_workspace, *w);
			argc2 = argc-1;
		} else {

			/* workspace is not given */
		}
	}

	switch (argc2) {
	case 2:
		CHECK_MATRIX_COMPLEX(argv[0]); CHECK_MATRIX_COMPLEX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix_complex, *A);
		Data_Get_Struct(argv[1], gsl_matrix_complex, *B);
		break;
	case 3:
		if (rb_obj_is_kind_of(argv[2], cgenhermv)) {
			Data_Get_Struct(argv[2], gsl_eigen_genhermv_workspace, *w);
		} else {
      rb_raise(rb_eTypeError, "Wrong argument type %s (GSL::Eigen::Genhermv::Workspace expected)",
        rb_class2name(CLASS_OF(argv[2])));
		}
		CHECK_MATRIX_COMPLEX(argv[0]); CHECK_MATRIX_COMPLEX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix_complex, *A);
		Data_Get_Struct(argv[1], gsl_matrix_complex, *B);

		break;
	case 5:
		if (rb_obj_is_kind_of(argv[4], cgenhermv)) {
			Data_Get_Struct(argv[4], gsl_eigen_genhermv_workspace, *w);
		} else {
      rb_raise(rb_eTypeError, "Wrong argument type %s (GSL::Eigen::Genhermv::Workspace expected)",
        rb_class2name(CLASS_OF(argv[4])));      
    }
		CHECK_VECTOR(argv[2]);
		Data_Get_Struct(argv[2], gsl_vector, *eval);
		CHECK_MATRIX_COMPLEX(argv[3]);
		Data_Get_Struct(argv[3], gsl_matrix_complex, *evec);			
		
		CHECK_MATRIX_COMPLEX(argv[0]); CHECK_MATRIX_COMPLEX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix_complex, *A);
		Data_Get_Struct(argv[1], gsl_matrix_complex, *B);
		break;		
	default:
		rb_raise(rb_eArgError, "wrong number of arguments (%d for 2, 3 or 5)", argc);
  }
	if (*eval == NULL && *evec == NULL) {
		*eval = gsl_vector_alloc((*A)->size1);
		*evec = gsl_matrix_complex_alloc((*A)->size1, (*A)->size2);		
		flag += 1;
	}
	if (*w == NULL) {
		*w = gsl_eigen_genhermv_alloc((*A)->size1);
		flag += 2;
	}
	return flag;
}
static VALUE rb_gsl_eigen_gensymm(int argc, VALUE *argv, VALUE obj)
{
	gsl_matrix *A = NULL, *B = NULL;
	gsl_matrix *Atmp = NULL;	
	gsl_vector *eval = NULL;
	gsl_eigen_gensymm_workspace *w = NULL;
	int flag;
	VALUE veval = NULL;
	flag = check_argv_gensymm(argc, argv, obj, &A, &B, &eval, &w);
  Atmp = make_matrix_clone(A);
//  Btmp = make_matrix_clone(B);  
	gsl_eigen_gensymm(Atmp, B, eval, w);
  gsl_matrix_free(Atmp);
//  gsl_matrix_free(Btmp);  
	switch (flag) {
	case 0:
		veval = argv[2];
		break;
	case 1:
		veval = Data_Wrap_Struct(cgsl_eigen_values, 0, gsl_vector_free, eval);
		break;
	case 2:
		veval = argv[2];
		gsl_eigen_gensymm_free(w);
		break;
	case 3:
		veval = Data_Wrap_Struct(cgsl_eigen_values, 0, gsl_vector_free, eval);
		gsl_eigen_gensymm_free(w);
		break;
	}

	return veval;
}

static VALUE rb_gsl_eigen_gensymmv(int argc, VALUE *argv, VALUE obj)
{
	gsl_matrix *A = NULL, *B = NULL;
	gsl_matrix *Atmp = NULL;
	gsl_vector *eval = NULL;
	gsl_matrix *evec = NULL;
	gsl_eigen_gensymmv_workspace *w = NULL;
	int flag;
	VALUE veval = NULL, vevec = NULL;
	flag = check_argv_gensymmv(argc, argv, obj, &A, &B, &eval, &evec, &w);
  Atmp = make_matrix_clone(A);
//  Btmp = make_matrix_clone(B);  	
	gsl_eigen_gensymmv(Atmp, B, eval, evec, w);
  gsl_matrix_free(Atmp);
//  gsl_matrix_free(Btmp);  
    
	switch (flag) {
	case 0:
		veval = argv[2];
		vevec = argv[3];		
		break;
	case 1:
		veval = Data_Wrap_Struct(cgsl_eigen_values, 0, gsl_vector_free, eval);
		vevec = Data_Wrap_Struct(cgsl_eigen_vectors, 0, gsl_matrix_free, evec);		
		break;
	case 2:
		veval = argv[2];
		vevec = argv[3];
		gsl_eigen_gensymmv_free(w);
	  break;
	case 3:
		veval = Data_Wrap_Struct(cgsl_eigen_values, 0, gsl_vector_free, eval);
		vevec = Data_Wrap_Struct(cgsl_eigen_vectors, 0, gsl_matrix_free, evec);	
		gsl_eigen_gensymmv_free(w);
		break;
	}

	return rb_ary_new3(2, veval, vevec);
}

static VALUE rb_gsl_eigen_genherm(int argc, VALUE *argv, VALUE obj)
{
	gsl_matrix_complex *A = NULL, *B = NULL;
	gsl_matrix_complex *Atmp = NULL, *Btmp = NULL;	
	gsl_vector *eval = NULL;
	gsl_eigen_genherm_workspace *w = NULL;
	int flag;
	VALUE veval = NULL;
	flag = check_argv_genherm(argc, argv, obj, &A, &B, &eval, &w);
  Atmp = make_matrix_complex_clone(A);
  Btmp = make_matrix_complex_clone(B);  
	gsl_eigen_genherm(Atmp, Btmp, eval, w);
  gsl_matrix_complex_free(Atmp);
  gsl_matrix_complex_free(Btmp);  
	switch (flag) {
	case 0:
		veval = argv[2];
		break;
	case 1:
		veval = Data_Wrap_Struct(cgsl_eigen_values, 0, gsl_vector_free, eval);
		break;
	case 2:
		veval = argv[2];
		gsl_eigen_genherm_free(w);
		break;
	case 3:
		veval = Data_Wrap_Struct(cgsl_eigen_values, 0, gsl_vector_free, eval);
		gsl_eigen_genherm_free(w);
		break;
	}

	return veval;
}

static VALUE rb_gsl_eigen_genhermv(int argc, VALUE *argv, VALUE obj)
{
	gsl_matrix_complex *A = NULL, *B = NULL;
	gsl_matrix_complex *Atmp = NULL, *Btmp = NULL;		
	gsl_vector *eval = NULL;
	gsl_matrix_complex *evec = NULL;
	gsl_eigen_genhermv_workspace *w = NULL;
	int flag;
	VALUE veval = NULL, vevec = NULL;
	flag = check_argv_genhermv(argc, argv, obj, &A, &B, &eval, &evec, &w);
  Atmp = make_matrix_complex_clone(A);
  Btmp = make_matrix_complex_clone(B);  	
	gsl_eigen_genhermv(Atmp, Btmp, eval, evec, w);
  gsl_matrix_complex_free(Atmp);
  gsl_matrix_complex_free(Btmp);  
    
	switch (flag) {
	case 0:
		veval = argv[2];
		vevec = argv[3];		
		break;
	case 1:
		veval = Data_Wrap_Struct(cgsl_eigen_values, 0, gsl_vector_free, eval);
		vevec = Data_Wrap_Struct(cgsl_eigen_herm_vectors, 0, gsl_matrix_complex_free, evec);		
		break;
	case 2:
		veval = argv[2];
		vevec = argv[3];
		gsl_eigen_genhermv_free(w);
	  break;
	case 3:
		veval = Data_Wrap_Struct(cgsl_eigen_values, 0, gsl_vector_free, eval);
		vevec = Data_Wrap_Struct(cgsl_eigen_herm_vectors, 0, gsl_matrix_complex_free, evec);	
		gsl_eigen_genhermv_free(w);
		break;
	}

	return rb_ary_new3(2, veval, vevec);
}

static VALUE rb_gsl_eigen_gensymmv_sort(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_eigen_real_sort(argc, argv, obj, gsl_eigen_gensymmv_sort);
}
static VALUE rb_gsl_eigen_genhermv_sort(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_eigen_complex_sort(argc, argv, obj, gsl_eigen_genhermv_sort);  
}

static VALUE rb_gsl_eigen_gen_alloc(VALUE klass, VALUE n)
{
  gsl_eigen_gen_workspace *w;
  w = gsl_eigen_gen_alloc(FIX2INT(n));
  return Data_Wrap_Struct(cgenw, 0, gsl_eigen_gen_free, w);
}

static VALUE rb_gsl_eigen_genv_alloc(VALUE klass, VALUE n)
{
  gsl_eigen_genv_workspace *w;
  w = gsl_eigen_genv_alloc(FIX2INT(n));
  return Data_Wrap_Struct(cgenvw, 0, gsl_eigen_genv_free, w);
}

static VALUE rb_gsl_eigen_gen_params(int argc, VALUE *argv, VALUE obj)
{
  gsl_eigen_gen_workspace *w = NULL;
  int istart = 0;
  if (CLASS_OF(obj) == cgenw) {
    Data_Get_Struct(obj, gsl_eigen_gen_workspace, w);
    istart = 0;
  } else {
    if (argc != 4) rb_raise(rb_eArgError, "too few arguments (%d for 3)\n", argc);
    if (CLASS_OF(argv[3]) != cgenw) 
      rb_raise(rb_eTypeError, "Wrong argument type %s (GSL::Eigen::Gen::Workspace expected)",
        rb_class2name(CLASS_OF(argv[3])));
    Data_Get_Struct(argv[3], gsl_eigen_gen_workspace, w);
    istart = 1;
  }
  switch (argc - istart) {
  case 3:
    gsl_eigen_gen_params(FIX2INT(argv[0]), FIX2INT(argv[1]), FIX2INT(argv[2]), w);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments.\n");
  }
  return Qtrue;
}

static int check_argv_gen(int argc, VALUE *argv, VALUE obj, gsl_matrix **A, gsl_matrix **B,
	gsl_vector_complex **alpha, gsl_vector **beta, gsl_eigen_gen_workspace **w)
{
	int argc2 = argc;
	int flag = 0;
	if (CLASS_OF(obj) == cgenw) {
		Data_Get_Struct(obj, gsl_eigen_gen_workspace, *w);
	} else {
		if (rb_obj_is_kind_of(argv[argc-1], cgenw)) {
			Data_Get_Struct(argv[argc-1], gsl_eigen_gen_workspace, *w);
			argc2 = argc-1;
		} else {

			/* workspace is not given */
		}
	}

	switch (argc2) {
	case 2:
		CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix, *A);
		Data_Get_Struct(argv[1], gsl_matrix, *B);
		break;
	case 3:
		if (rb_obj_is_kind_of(argv[2], cgenw)) {
			Data_Get_Struct(argv[2], gsl_eigen_gen_workspace, *w);
		} else {
      rb_raise(rb_eTypeError, "Wrong argument type %s (GSL::Eigen::Gen::Workspace expected)",
        rb_class2name(CLASS_OF(argv[2])));
		}
		CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix, *A);
		Data_Get_Struct(argv[1], gsl_matrix, *B);
		break;
	case 5:
		if (rb_obj_is_kind_of(argv[4], cgenw)) {
			Data_Get_Struct(argv[4], gsl_eigen_gen_workspace, *w);
		} else {
      rb_raise(rb_eTypeError, "Wrong argument type %s (GSL::Eigen::Gen::Workspace expected)",
        rb_class2name(CLASS_OF(argv[4])));      
    }
		CHECK_VECTOR_COMPLEX(argv[2]);
		Data_Get_Struct(argv[2], gsl_vector_complex, *alpha);
		CHECK_VECTOR(argv[3]);
		Data_Get_Struct(argv[3], gsl_vector, *beta);			
		
		CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix, *A);
		Data_Get_Struct(argv[1], gsl_matrix, *B);
		break;		
	default:
		rb_raise(rb_eArgError, "wrong number of arguments (%d for 2, 3 or 5)", argc);
  }
	if (*alpha == NULL && *beta == NULL) {
		*alpha = gsl_vector_complex_alloc((*A)->size1);
		*beta = gsl_vector_alloc((*A)->size1);
		flag += 1;
	}
	if (*w == NULL) {
		*w = gsl_eigen_gen_alloc((*A)->size1);
		flag += 2;
	}
	return flag;
}

static VALUE rb_gsl_eigen_gen(int argc, VALUE *argv, VALUE obj)
{
	gsl_matrix *A = NULL, *B = NULL;
//	gsl_matrix *Atmp = NULL, *Btmp = NULL;		
	gsl_vector_complex *alpha = NULL;
	gsl_vector *beta = NULL;
	gsl_eigen_gen_workspace *w = NULL;
	int flag;
	VALUE valpha = NULL, vbeta = NULL;
	flag = check_argv_gen(argc, argv, obj, &A, &B, &alpha, &beta, &w);
//  Atmp = make_matrix_clone(A);
//  Btmp = make_matrix_clone(B);  	
	gsl_eigen_gen(A, B, alpha, beta, w);
//  gsl_matrix_free(Atmp);
//  gsl_matrix_free(Btmp);  
    
	switch (flag) {
	case 0:
		valpha = argv[2];
		vbeta = argv[3];		
		break;
	case 1:
		valpha = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, alpha);
		vbeta = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, beta);		
		break;
	case 2:
		valpha = argv[2];
		vbeta = argv[3];
		gsl_eigen_gen_free(w);
	  break;
	case 3:
		valpha = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, alpha);
		vbeta = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, beta);	
		gsl_eigen_gen_free(w);
		break;
	}

	return rb_ary_new3(2, valpha, vbeta);
}

static VALUE rb_gsl_eigen_gen_QZ(int argc, VALUE *argv, VALUE obj)
{
	gsl_matrix *A = NULL, *B = NULL;
//	gsl_matrix *Atmp = NULL, *Btmp = NULL;		
	gsl_vector_complex *alpha = NULL;
	gsl_vector *beta = NULL;
	gsl_matrix *Q, *Z;
	gsl_eigen_gen_workspace *w = NULL;
	int flag;
	VALUE valpha = NULL, vbeta = NULL, vQ, vZ;
	flag = check_argv_gen(argc, argv, obj, &A, &B, &alpha, &beta, &w);
/*  Atmp = make_matrix_clone(A);
  Btmp = make_matrix_clone(B);  	*/
  Q = gsl_matrix_alloc(A->size1, A->size2);
  Z = gsl_matrix_alloc(A->size1, A->size2);  
//	gsl_eigen_gen_QZ(Atmp, Btmp, alpha, beta, Q, Z, w);
	gsl_eigen_gen_QZ(A, B, alpha, beta, Q, Z, w);	
/*  gsl_matrix_free(Atmp);
  gsl_matrix_free(Btmp);  */
    
	switch (flag) {
	case 0:
		valpha = argv[2];
		vbeta = argv[3];		
		break;
	case 1:
		valpha = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, alpha);
		vbeta = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, beta);		
		break;
	case 2:
		valpha = argv[2];
		vbeta = argv[3];
		gsl_eigen_gen_free(w);
	  break;
	case 3:
		valpha = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, alpha);
		vbeta = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, beta);	
		gsl_eigen_gen_free(w);
		break;
	}
  vQ = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Q);
  vZ = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Z);  
	return rb_ary_new3(4, valpha, vbeta, vQ, vZ);
}

static int check_argv_genv(int argc, VALUE *argv, VALUE obj, gsl_matrix **A, gsl_matrix **B,
	gsl_vector_complex **alpha, gsl_vector **beta, gsl_matrix_complex **evec, gsl_eigen_genv_workspace **w)
{
	int argc2 = argc;
	int flag = 0;
	if (CLASS_OF(obj) == cgenvw) {
		Data_Get_Struct(obj, gsl_eigen_genv_workspace, *w);
	} else {
		if (rb_obj_is_kind_of(argv[argc-1], cgenvw)) {
			Data_Get_Struct(argv[argc-1], gsl_eigen_genv_workspace, *w);
			argc2 = argc-1;
		} else {

			/* workspace is not given */
		}
	}

	switch (argc2) {
	case 2:
		CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix, *A);
		Data_Get_Struct(argv[1], gsl_matrix, *B);
		break;
	case 3:
		if (rb_obj_is_kind_of(argv[2], cgenvw)) {
			Data_Get_Struct(argv[2], gsl_eigen_genv_workspace, *w);
		} else {
      rb_raise(rb_eTypeError, "Wrong argument type %s (GSL::Eigenv::Gen::Workspace expected)",
        rb_class2name(CLASS_OF(argv[2])));
		}
		CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix, *A);
		Data_Get_Struct(argv[1], gsl_matrix, *B);
		break;
	case 6:
		if (rb_obj_is_kind_of(argv[4], cgenvw)) {
			Data_Get_Struct(argv[4], gsl_eigen_genv_workspace, *w);
		} else {
      rb_raise(rb_eTypeError, "Wrong argument type %s (GSL::Eigenv::Gen::Workspace expected)",
        rb_class2name(CLASS_OF(argv[4])));      
    }
		CHECK_VECTOR_COMPLEX(argv[2]);
		Data_Get_Struct(argv[2], gsl_vector_complex, *alpha);
		CHECK_VECTOR(argv[3]);
		Data_Get_Struct(argv[3], gsl_vector, *beta);			
		CHECK_MATRIX_COMPLEX(argv[3]);
		Data_Get_Struct(argv[4], gsl_matrix_complex, *evec);					
		
		CHECK_MATRIX(argv[0]); CHECK_MATRIX(argv[1]);
		Data_Get_Struct(argv[0], gsl_matrix, *A);
		Data_Get_Struct(argv[1], gsl_matrix, *B);
		break;		
	default:
		rb_raise(rb_eArgError, "wrong number of arguments (%d for 2, 3 or 6)", argc);
  }
	if (*alpha == NULL && *beta == NULL) {
		*alpha = gsl_vector_complex_alloc((*A)->size1);
		*beta = gsl_vector_alloc((*A)->size1);
		*evec = gsl_matrix_complex_alloc((*A)->size1, (*A)->size2);		
		flag += 1;
	}
	if (*w == NULL) {
		*w = gsl_eigen_genv_alloc((*A)->size1);
		flag += 2;
	}
	return flag;
}

static VALUE rb_gsl_eigen_genv(int argc, VALUE *argv, VALUE obj)
{
	gsl_matrix *A = NULL, *B = NULL;
//	gsl_matrix *Atmp = NULL, *Btmp = NULL;		
	gsl_vector_complex *alpha = NULL;
	gsl_vector *beta = NULL;
	gsl_matrix_complex *evec = NULL;
	gsl_eigen_genv_workspace *w = NULL;
	int flag;
	VALUE valpha = NULL, vbeta = NULL, vevec = NULL;
	flag = check_argv_genv(argc, argv, obj, &A, &B, &alpha, &beta, &evec, &w);
//  Atmp = make_matrix_clone(A);
//  Btmp = make_matrix_clone(B);  	
//	gsl_eigen_genv(Atmp, Btmp, alpha, beta, evec, w);
	gsl_eigen_genv(A, B, alpha, beta, evec, w);	
//  gsl_matrix_free(Atmp);
//  gsl_matrix_free(Btmp);  
    
	switch (flag) {
	case 0:
		valpha = argv[2];
		vbeta = argv[3];
		vevec = argv[4];
		break;
	case 1:
		valpha = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, alpha);
		vbeta = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, beta);
		vevec = Data_Wrap_Struct(cgsl_eigen_herm_vectors, 0, gsl_matrix_complex_free, evec);		
		break;
	case 2:
		valpha = argv[2];
		vbeta = argv[3];
		vevec = argv[4];		
		gsl_eigen_genv_free(w);
	  break;
	case 3:
		valpha = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, alpha);
		vbeta = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, beta);
		vevec = Data_Wrap_Struct(cgsl_eigen_herm_vectors, 0, gsl_matrix_complex_free, evec);
		gsl_eigen_genv_free(w);
		break;
	}

	return rb_ary_new3(3, valpha, vbeta, vevec);
}


static VALUE rb_gsl_eigen_genv_QZ(int argc, VALUE *argv, VALUE obj)
{
	gsl_matrix *A = NULL, *B = NULL;
//	gsl_matrix *Atmp = NULL, *Btmp = NULL;		
	gsl_vector_complex *alpha = NULL;
	gsl_vector *beta = NULL;
	gsl_matrix_complex *evec = NULL;
	gsl_matrix *Q, *Z;
	gsl_eigen_genv_workspace *w = NULL;
	int flag;
	VALUE valpha = NULL, vbeta = NULL, vevec = NULL, vQ, vZ;
	flag = check_argv_genv(argc, argv, obj, &A, &B, &alpha, &beta, &evec, &w);
/*  Atmp = make_matrix_clone(A);
  Btmp = make_matrix_clone(B);  	*/
  Q = gsl_matrix_alloc(A->size1, A->size2);
  Z = gsl_matrix_alloc(A->size1, A->size2);  
//	gsl_eigen_genv_QZ(Atmp, Btmp, alpha, beta, evec, Q, Z, w);
	gsl_eigen_genv_QZ(A, B, alpha, beta, evec, Q, Z, w);	
/*  gsl_matrix_free(Atmp);
  gsl_matrix_free(Btmp);  */
    
	switch (flag) {
	case 0:
		valpha = argv[2];
		vbeta = argv[3];
		vevec = argv[4];
		break;
	case 1:
		valpha = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, alpha);
		vbeta = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, beta);
		vevec = Data_Wrap_Struct(cgsl_eigen_herm_vectors, 0, gsl_matrix_complex_free, evec);		
		break;
	case 2:
		valpha = argv[2];
		vbeta = argv[3];
		vevec = argv[4];		
		gsl_eigen_genv_free(w);
	  break;
	case 3:
		valpha = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, alpha);
		vbeta = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, beta);
		vevec = Data_Wrap_Struct(cgsl_eigen_herm_vectors, 0, gsl_matrix_complex_free, evec);
		gsl_eigen_genv_free(w);
		break;
	}
  vQ = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Q);
  vZ = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, Z);  
	return rb_ary_new3(5, valpha, vbeta, vevec, vQ, vZ);
}

static VALUE rb_gsl_eigen_genv_sort(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector_complex *alpha = NULL;
  gsl_vector *beta = NULL;
  gsl_matrix_complex *evec = NULL;
  gsl_eigen_sort_t type = GSL_EIGEN_SORT_VAL_DESC;

  switch (argc) {
  case 4:
    CHECK_FIXNUM(argv[3]);
    type = FIX2INT(argv[3]);
    /* no break, do next */
  case 3:
    if (argv[0] == Qnil) {
      alpha = NULL;
    } else {
      CHECK_VECTOR_COMPLEX(argv[0]);
      Data_Get_Struct(argv[0], gsl_vector_complex, alpha);
    }
    if (argv[1] == Qnil) {
      beta = NULL;
    } else {
      CHECK_VECTOR(argv[1]);
      Data_Get_Struct(argv[1], gsl_vector, beta);        
    }
   if (argv[2] == Qnil) {
      evec = NULL;
    } else {
      CHECK_MATRIX_COMPLEX(argv[2]);
      Data_Get_Struct(argv[2], gsl_matrix_complex, evec);        
    }    
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 3 or 4)", argc);
  }
  return INT2FIX(gsl_eigen_genv_sort(alpha, beta, evec, type));
}
#endif

void Init_gsl_eigen(VALUE module)
{
  VALUE mgsl_eigen;
  VALUE mgsl_eigen_symm;
  VALUE mgsl_eigen_symmv;
  VALUE mgsl_eigen_herm;
  VALUE mgsl_eigen_hermv;
#ifdef HAVE_EIGEN_FRANCIS
  VALUE mgsl_eigen_francis;
#endif
#ifdef GSL_1_9_LATER
  VALUE mgsl_eigen_nonsymmv;
  VALUE mgsl_eigen_nonsymm;
#endif

  mgsl_eigen = rb_define_module_under(module, "Eigen");
  mgsl_eigen_symm = rb_define_module_under(mgsl_eigen, "Symm");
  mgsl_eigen_symmv = rb_define_module_under(mgsl_eigen, "Symmv");

  mgsl_eigen_herm = rb_define_module_under(mgsl_eigen, "Herm");
  mgsl_eigen_hermv = rb_define_module_under(mgsl_eigen, "Hermv");

  cgsl_eigen_values = rb_define_class_under(mgsl_eigen, "EigenValues",
					    cgsl_vector);
  cgsl_eigen_vectors = rb_define_class_under(mgsl_eigen, "EigenVectors",
					    cgsl_matrix);
  cgsl_eigen_vector = rb_define_class_under(mgsl_eigen, "EigenVector",
					    cgsl_vector);
  cgsl_eigen_herm_vectors = rb_define_class_under(mgsl_eigen, "ComplexEigenVectors",
						  cgsl_matrix_complex);
  cgsl_eigen_vector_complex = rb_define_class_under(mgsl_eigen, "ComplexEigenVector",
						  cgsl_vector_complex);
  cgsl_eigen_symm_workspace = rb_define_class_under(mgsl_eigen_symm, 
						     "Workspace", cGSL_Object);
  cgsl_eigen_symmv_workspace = rb_define_class_under(mgsl_eigen_symmv, 
						     "Workspace", cGSL_Object);
  cgsl_eigen_herm_workspace = rb_define_class_under(mgsl_eigen_herm, 
						     "Workspace", cGSL_Object);

  cgsl_eigen_hermv_workspace = rb_define_class_under(mgsl_eigen_hermv, 
						     "Workspace", cGSL_Object);

  rb_define_singleton_method(cgsl_eigen_symm_workspace, "alloc", 
			     rb_gsl_eigen_symm_alloc, 1);
  rb_define_singleton_method(cgsl_eigen_symmv_workspace, "alloc", 
			     rb_gsl_eigen_symmv_alloc, 1);
  rb_define_singleton_method(cgsl_eigen_herm_workspace, "alloc", 
			     rb_gsl_eigen_herm_alloc, 1);
  rb_define_singleton_method(cgsl_eigen_hermv_workspace, "alloc", 
			     rb_gsl_eigen_hermv_alloc, 1);

  rb_define_singleton_method(mgsl_eigen_symm, "alloc", 
			     rb_gsl_eigen_symm_alloc, 1);
  rb_define_singleton_method(mgsl_eigen_symmv, "alloc", 
			     rb_gsl_eigen_symmv_alloc, 1);
  rb_define_singleton_method(mgsl_eigen_herm, "alloc", 
			     rb_gsl_eigen_herm_alloc, 1);
  rb_define_singleton_method(mgsl_eigen_hermv, "alloc", 
			     rb_gsl_eigen_hermv_alloc, 1);
			     
  rb_define_module_function(mgsl_eigen, "symm",
			     rb_gsl_eigen_symm, -1);
  rb_define_module_function(mgsl_eigen, "symmv",
			     rb_gsl_eigen_symmv, -1);
  rb_define_module_function(mgsl_eigen, "herm",
			     rb_gsl_eigen_herm, -1);
  rb_define_module_function(mgsl_eigen, "hermv",
			     rb_gsl_eigen_hermv, -1);

  rb_define_module_function(module, "eigen_symm",
			     rb_gsl_eigen_symm, -1);
  rb_define_module_function(module, "eigen_symmv",
			     rb_gsl_eigen_symmv, -1);
  rb_define_module_function(module, "eigen_herm",
			     rb_gsl_eigen_herm, -1);
  rb_define_module_function(module, "eigen_hermv",
			     rb_gsl_eigen_hermv, -1);

  rb_define_method(cgsl_matrix, "eigen_symm", rb_gsl_eigen_symm, -1);
  rb_define_method(cgsl_matrix, "eigen_symmv", rb_gsl_eigen_symmv, -1);
  rb_define_method(cgsl_matrix_complex, "eigen_herm", rb_gsl_eigen_herm, -1);
  rb_define_method(cgsl_matrix_complex, "eigen_hermv", rb_gsl_eigen_hermv, -1);

  rb_define_method(cgsl_eigen_vectors, "unpack", rb_gsl_eigen_vectors_unpack, 0);
  rb_define_method(cgsl_eigen_herm_vectors, "unpack", rb_gsl_eigen_vectors_complex_unpack, 0);

  rb_gsl_eigen_define_const(module, mgsl_eigen);

  rb_define_module_function(mgsl_eigen, "symmv_sort", 
			     rb_gsl_eigen_symmv_sort, -1);
  rb_define_module_function(mgsl_eigen, "hermv_sort", 
			     rb_gsl_eigen_hermv_sort, -1);
  rb_define_module_function(mgsl_eigen_symmv, "sort", 
			     rb_gsl_eigen_symmv_sort, -1);
  rb_define_module_function(mgsl_eigen_hermv, "sort", 
			     rb_gsl_eigen_hermv_sort, -1);

#ifdef HAVE_EIGEN_FRANCIS
  mgsl_eigen_francis = rb_define_module_under(mgsl_eigen, "francis");
  cgsl_eigen_francis_workspace = rb_define_class_under(mgsl_eigen_francis, 
						    "Workspace", cGSL_Object);
  rb_define_singleton_method(cgsl_eigen_francis_workspace, "alloc",
			     rb_gsl_eigen_francis_alloc, 0);

  rb_define_method(cgsl_matrix, "eigen_francis", rb_gsl_eigen_francis, -1);
  rb_define_module_function(mgsl_eigen, "francis", rb_gsl_eigen_francis, -1);
  rb_define_module_function(module, "eigen_francis", rb_gsl_eigen_francis, -1);
  rb_define_method(cgsl_matrix, "eigen_francis_Z", rb_gsl_eigen_francis_Z, -1);
  rb_define_module_function(mgsl_eigen, "francis_Z", rb_gsl_eigen_francis_Z, -1);
  rb_define_module_function(module, "eigen_francis_Z", rb_gsl_eigen_francis_Z, -1);

  rb_define_method(cgsl_eigen_francis_workspace, "T", rb_gsl_eigen_francis_T, -1);
  rb_define_module_function(mgsl_eigen_francis, "T", rb_gsl_eigen_francis_T, -1);

#endif

#ifdef GSL_1_9_LATER
  mgsl_eigen_nonsymm = rb_define_module_under(mgsl_eigen, "Nonsymm");
  mgsl_eigen_nonsymmv = rb_define_module_under(mgsl_eigen, "Nonsymmv");
  cgsl_eigen_nonsymm_workspace = rb_define_class_under(mgsl_eigen_nonsymm, 
						    "Workspace", cGSL_Object);
  rb_define_singleton_method(cgsl_eigen_nonsymm_workspace, "alloc",
			     rb_gsl_eigen_nonsymm_alloc, 1);
  rb_define_singleton_method(mgsl_eigen_nonsymm, "alloc",
			     rb_gsl_eigen_nonsymm_alloc, 1);
			     
  rb_define_method(cgsl_matrix, "eigen_nonsymm", rb_gsl_eigen_nonsymm, -1);
  rb_define_module_function(mgsl_eigen, "nonsymm", rb_gsl_eigen_nonsymm, -1);
  rb_define_module_function(module, "eigen_nonsymm", rb_gsl_eigen_nonsymm, -1);
  rb_define_method(cgsl_matrix, "eigen_nonsymm_Z", rb_gsl_eigen_nonsymm_Z, -1);
  rb_define_module_function(mgsl_eigen, "nonsymm_Z", rb_gsl_eigen_nonsymm_Z, -1);
  rb_define_module_function(module, "eigen_nonsymm_Z", rb_gsl_eigen_nonsymm_Z, -1);

  rb_define_method(cgsl_eigen_nonsymm_workspace, "params", rb_gsl_eigen_nonsymm_params, -1);
  rb_define_module_function(mgsl_eigen_nonsymm, "params", rb_gsl_eigen_nonsymm_params, -1);

  cgsl_eigen_nonsymmv_workspace = rb_define_class_under(mgsl_eigen_nonsymmv, 
						    "Workspace", cGSL_Object);
  rb_define_singleton_method(cgsl_eigen_nonsymmv_workspace, "alloc",
			     rb_gsl_eigen_nonsymmv_alloc, 1);
  rb_define_singleton_method(mgsl_eigen_nonsymmv, "alloc",
			     rb_gsl_eigen_nonsymmv_alloc, 1);
  rb_define_method(cgsl_matrix, "eigen_nonsymmv", rb_gsl_eigen_nonsymmv, -1);
  rb_define_module_function(mgsl_eigen, "nonsymmv", rb_gsl_eigen_nonsymmv, -1);
  rb_define_module_function(module, "eigen_nonsymmv", rb_gsl_eigen_nonsymmv, -1);
  rb_define_method(cgsl_matrix, "eigen", rb_gsl_eigen_nonsymmv, -1);
  rb_define_alias(cgsl_matrix, "eig", "eigen");
  rb_define_method(cgsl_matrix, "eigen_nonsymmv_Z", rb_gsl_eigen_nonsymmv_Z, -1);
  rb_define_module_function(mgsl_eigen, "nonsymmv_Z", rb_gsl_eigen_nonsymmv_Z, -1);
  rb_define_module_function(module, "eigen_nonsymmv_Z", rb_gsl_eigen_nonsymmv_Z, -1);

  rb_define_module_function(mgsl_eigen, "nonsymmv_sort", 
			     rb_gsl_eigen_nonsymmv_sort, -1);
  rb_define_module_function(mgsl_eigen_nonsymmv, "sort", 
			     rb_gsl_eigen_nonsymmv_sort, -1);
  rb_define_module_function(module, "eigen_nonsymmv_sort", 
			     rb_gsl_eigen_nonsymmv_sort, -1);
#endif

#ifdef GSL_1_10_LATER
  /** gensymm, gensymmv **/
  mgensymm = rb_define_module_under(mgsl_eigen, "Gensymm");
  cgensymm = rb_define_class_under(mgensymm, "Workspace", cGSL_Object);
  mgensymmv = rb_define_module_under(mgsl_eigen, "Gensymmv");
  cgensymmv = rb_define_class_under(mgensymmv, "Workspace", cGSL_Object);

  rb_define_singleton_method(cgensymm, "alloc", rb_gsl_eigen_gensymm_alloc, 1);
  rb_define_singleton_method(cgensymmv, "alloc", rb_gsl_eigen_gensymmv_alloc, 1);
  rb_define_singleton_method(mgensymm, "alloc", rb_gsl_eigen_gensymm_alloc, 1);
  rb_define_singleton_method(mgensymmv, "alloc", rb_gsl_eigen_gensymmv_alloc, 1);

  rb_define_method(cgensymm, "gensymm", rb_gsl_eigen_gensymm, -1);
  rb_define_module_function(module, "eigen_gensymm", rb_gsl_eigen_gensymm, -1);
  rb_define_module_function(mgsl_eigen, "gensymm", rb_gsl_eigen_gensymm, -1);	
  rb_define_module_function(mgensymm, "gensymm", rb_gsl_eigen_gensymm, -1);	

  rb_define_method(cgensymmv, "gensymmv", rb_gsl_eigen_gensymmv, -1);
  rb_define_module_function(module, "eigen_gensymmv", rb_gsl_eigen_gensymmv, -1);
  rb_define_module_function(mgsl_eigen, "gensymmv", rb_gsl_eigen_gensymmv, -1);	
  rb_define_module_function(mgensymmv, "gensymmv", rb_gsl_eigen_gensymmv, -1);

  rb_define_module_function(mgsl_eigen, "gensymmv_sort", 
         rb_gsl_eigen_gensymmv_sort, -1);
  rb_define_module_function(mgensymmv, "sort", 
         rb_gsl_eigen_gensymmv_sort, -1);
  rb_define_module_function(module, "eigen_gensymmv_sort", 
         rb_gsl_eigen_gensymmv_sort, -1);

  /** genherm, genhermv **/
  mgenherm = rb_define_module_under(mgsl_eigen, "Genherm");
  cgenherm = rb_define_class_under(mgenherm, "Workspace", cGSL_Object);
  mgenhermv = rb_define_module_under(mgsl_eigen, "Genhermv");
  cgenhermv = rb_define_class_under(mgenhermv, "Workspace", cGSL_Object);
  
  rb_define_singleton_method(cgenherm, "alloc", rb_gsl_eigen_genherm_alloc, 1);
  rb_define_singleton_method(cgenhermv, "alloc", rb_gsl_eigen_genhermv_alloc, 1);
  rb_define_singleton_method(mgenherm, "alloc", rb_gsl_eigen_genherm_alloc, 1);
  rb_define_singleton_method(mgenhermv, "alloc", rb_gsl_eigen_genhermv_alloc, 1);

  rb_define_method(cgenherm, "genherm", rb_gsl_eigen_genherm, -1);
  rb_define_module_function(module, "eigen_genherm", rb_gsl_eigen_genherm, -1);
  rb_define_module_function(mgsl_eigen, "genherm", rb_gsl_eigen_genherm, -1);	
  rb_define_module_function(mgenherm, "genherm", rb_gsl_eigen_genherm, -1);	

  rb_define_method(cgenhermv, "genhermv", rb_gsl_eigen_genhermv, -1);
  rb_define_module_function(module, "eigen_genhermv", rb_gsl_eigen_genhermv, -1);
  rb_define_module_function(mgsl_eigen, "genhermv", rb_gsl_eigen_genhermv, -1);	
  rb_define_module_function(mgenhermv, "genhermv", rb_gsl_eigen_genhermv, -1);
  
  rb_define_module_function(mgsl_eigen, "genhermv_sort", 
         rb_gsl_eigen_genhermv_sort, -1);
  rb_define_module_function(mgenhermv, "sort", 
         rb_gsl_eigen_genhermv_sort, -1);
  rb_define_module_function(module, "eigen_genhermv_sort", 
         rb_gsl_eigen_genhermv_sort, -1); 
         
  /* gen */
  mgen = rb_define_module_under(mgsl_eigen, "Gen");
  mgenv = rb_define_module_under(mgsl_eigen, "Genv");  
  cgenw = rb_define_class_under(mgen, "Workspace", cGSL_Object);
  cgenvw = rb_define_class_under(mgenv, "Workspace", cGSL_Object);
  rb_define_singleton_method(mgen, "alloc", rb_gsl_eigen_gen_alloc, 1);
  rb_define_singleton_method(cgenw, "alloc", rb_gsl_eigen_gen_alloc, 1);  
  rb_define_singleton_method(mgenv, "alloc", rb_gsl_eigen_genv_alloc, 1);
  rb_define_singleton_method(cgenvw, "alloc", rb_gsl_eigen_genv_alloc, 1);
  
  rb_define_module_function(mgen, "params", rb_gsl_eigen_gen_params, -1);
  rb_define_method(cgenw, "params", rb_gsl_eigen_gen_params, -1);
  rb_define_module_function(mgsl_eigen, "gen_params", rb_gsl_eigen_gen_params, -1);
  
    
  rb_define_module_function(mgen, "gen", rb_gsl_eigen_gen, -1);
  rb_define_module_function(mgsl_eigen, "gen", rb_gsl_eigen_gen, -1);  
  rb_define_method(cgenw, "gen", rb_gsl_eigen_gen, -1);
  
  rb_define_module_function(mgenv, "genv", rb_gsl_eigen_genv, -1);
  rb_define_module_function(mgsl_eigen, "genv", rb_gsl_eigen_genv, -1);  
  rb_define_method(cgenvw, "genv", rb_gsl_eigen_genv, -1);
  
  rb_define_module_function(mgen, "gen_QZ", rb_gsl_eigen_gen_QZ, -1);
  rb_define_module_function(mgsl_eigen, "gen_QZ", rb_gsl_eigen_gen_QZ, -1);  
  rb_define_method(cgenw, "gen_QZ", rb_gsl_eigen_gen_QZ, -1);
  
   rb_define_module_function(mgenv, "genv_QZ", rb_gsl_eigen_genv_QZ, -1);
  rb_define_module_function(mgsl_eigen, "genv_QZ", rb_gsl_eigen_genv_QZ, -1);  
  rb_define_method(cgenvw, "genv_QZ", rb_gsl_eigen_genv_QZ, -1);  
  
  rb_define_module_function(mgsl_eigen, "genv_sort", 
         rb_gsl_eigen_genv_sort, -1);
  rb_define_module_function(mgenv, "sort", 
         rb_gsl_eigen_genv_sort, -1);
  rb_define_module_function(module, "eigen_genv_sort", 
         rb_gsl_eigen_genv_sort, -1);   
#endif

}

