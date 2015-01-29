/*
  wavelet.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_array.h"
#include "include/rb_gsl_common.h"

#ifdef GSL_1_6_LATER
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_wavelet2d.h>
#endif

#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif

#ifndef CHECK_WAVELET
#define CHECK_WAVELET(x) if(!rb_obj_is_kind_of(x,cgsl_wavelet))\
    rb_raise(rb_eTypeError, "wrong argument type (Wavelet expected)");
#endif

#ifndef CHECK_WORKSPACE
#define CHECK_WORKSPACE(x) if(!rb_obj_is_kind_of(x,cgsl_wavelet_workspace))\
    rb_raise(rb_eTypeError, "wrong argument type (Wavelet::Workspace expected)");
#endif

enum RB_GSL_DWT {
  RB_GSL_DWT_COPY,
  RB_GSL_DWT_INPLACE,
};

static VALUE cgsl_wavelet;

enum {
  GSL_WAVELET_DAUBECHIES,
  GSL_WAVELET_DAUBECHIES_CENTERED,
  GSL_WAVELET_HAAR,
  GSL_WAVELET_HAAR_CENTERED,
  GSL_WAVELET_BSPLINE,
  GSL_WAVELET_BSPLINE_CENTERED,
};

#ifdef GSL_1_6_LATER
static const gsl_wavelet_type* rb_gsl_wavelet_get_type(VALUE t);
static VALUE cgsl_wavelet_workspace;
#endif

static VALUE rb_gsl_wavelet_new(VALUE klass, VALUE t, VALUE m)
{
#ifdef GSL_1_6_LATER
  const gsl_wavelet_type *T;
  size_t member;
  gsl_wavelet *w = NULL;
  CHECK_FIXNUM(m);
  T = rb_gsl_wavelet_get_type(t);
  member = FIX2INT(m);
  w = gsl_wavelet_alloc(T, member);
  if (w == NULL) rb_raise(rb_eNoMemError, "gsl_wavelet_alloc failed");
  return Data_Wrap_Struct(klass, 0, gsl_wavelet_free, w);
#else
  rb_raise(rb_eNotImpError, "Wavelet transforms not supported in GSL-%s, use GSL-1.6 or later", GSL_VERSION);
  return Qnil;
#endif
}

#ifdef GSL_1_6_LATER
static const gsl_wavelet_type* rb_gsl_wavelet_get_type_str(char *name);
static const gsl_wavelet_type* rb_gsl_wavelet_get_type_int(int t);
static const gsl_wavelet_type* rb_gsl_wavelet_get_type(VALUE t)
{
  const gsl_wavelet_type *T;
  switch (TYPE(t)) {
  case T_STRING:
    T = rb_gsl_wavelet_get_type_str(STR2CSTR(t));
    break;
  case T_FIXNUM:
    T = rb_gsl_wavelet_get_type_int(FIX2INT(t));
    break;
  default:
    rb_raise(rb_eTypeError, 
       "wrong type of argument %s (String or Fixnum expected)",
       rb_class2name(CLASS_OF(t)));
    break;
  }
  return T;
}

static const gsl_wavelet_type* rb_gsl_wavelet_get_type_str(char *name)
{
  const gsl_wavelet_type *T;
  if (str_tail_grep(name, "daubechies") == 0)
    T = gsl_wavelet_daubechies;
  else if (str_tail_grep(name, "daubechies_centered") == 0)
    T = gsl_wavelet_daubechies_centered;
  else if (str_tail_grep(name, "haar") == 0)
    T = gsl_wavelet_haar;
  else if (str_tail_grep(name, "haar_centered") == 0)
    T = gsl_wavelet_haar_centered;
  else if (str_tail_grep(name, "bspline") == 0)
    T = gsl_wavelet_bspline;
  else if (str_tail_grep(name, "bspline_centered") == 0)
    T = gsl_wavelet_bspline_centered;
  else
    rb_raise(rb_eArgError, "unknown type %s", name);
  return T;
}

static const gsl_wavelet_type* rb_gsl_wavelet_get_type_int(int t)
{
  const gsl_wavelet_type *T;
  switch (t) {
  case GSL_WAVELET_DAUBECHIES:
    T = gsl_wavelet_daubechies;
    break;
  case GSL_WAVELET_DAUBECHIES_CENTERED:
    T = gsl_wavelet_daubechies_centered;
    break;
  case GSL_WAVELET_HAAR:
    T = gsl_wavelet_haar;
    break;
  case GSL_WAVELET_HAAR_CENTERED:
    T = gsl_wavelet_haar_centered;
    break;
  case GSL_WAVELET_BSPLINE:
    T = gsl_wavelet_bspline;
    break;
  case GSL_WAVELET_BSPLINE_CENTERED:
    T = gsl_wavelet_bspline_centered;
    break;
  default:
    rb_raise(rb_eArgError, "unknown type %d", t);
    break;
  }
  return T;
}

static void rb_gsl_wavelet_define_const(VALUE klass);
static void rb_gsl_wavelet_define_const(VALUE klass)
{
  rb_define_const(klass, "DAUBECHIES", INT2FIX(GSL_WAVELET_DAUBECHIES));
  rb_define_const(klass, "DAUBECHIES_CENTERED", INT2FIX(GSL_WAVELET_DAUBECHIES_CENTERED));
  rb_define_const(klass, "HAAR", INT2FIX(GSL_WAVELET_HAAR));
  rb_define_const(klass, "HAAR_CENTERED", INT2FIX(GSL_WAVELET_HAAR_CENTERED));
  rb_define_const(klass, "BSPLINE", INT2FIX(GSL_WAVELET_BSPLINE));
  rb_define_const(klass, "BSPLINE_CENTERED", INT2FIX(GSL_WAVELET_BSPLINE_CENTERED));
  /*****/
  rb_define_const(klass, "FORWARD", INT2FIX(gsl_wavelet_forward));
  rb_define_const(klass, "Forward", INT2FIX(gsl_wavelet_forward));
  rb_define_const(klass, "BACKWARD", INT2FIX(gsl_wavelet_backward));
  rb_define_const(klass, "Backward", INT2FIX(gsl_wavelet_backward));
}

static VALUE rb_gsl_wavelet_name(VALUE ww)
{
  gsl_wavelet *w = NULL;
  Data_Get_Struct(ww, gsl_wavelet, w);
  return rb_str_new2(gsl_wavelet_name(w));
}

static VALUE rb_gsl_wavelet_workspace_new(VALUE klass, VALUE nn)
{
  gsl_wavelet_workspace *wspace = NULL;
  CHECK_FIXNUM(nn);
  wspace = gsl_wavelet_workspace_alloc(FIX2INT(nn));
  if (wspace == NULL) rb_raise(rb_eNoMemError, "gsl_wavelet_workspace_alloc failed");
  return Data_Wrap_Struct(klass, 0, gsl_wavelet_workspace_free, wspace);
}

static VALUE rb_gsl_wavelet2d_trans(int argc, VALUE *argv, VALUE obj,
            int (*trans)(const gsl_wavelet *, 
             gsl_matrix *,
             gsl_wavelet_workspace *),
            int sss);
static VALUE rb_gsl_wavelet2d(int argc, VALUE *argv, VALUE obj,
            int (*trans)(const gsl_wavelet *, 
             gsl_matrix *, 
             gsl_wavelet_direction, 
             gsl_wavelet_workspace *),
            int sss);

static VALUE rb_gsl_wavelet_transform0(int argc, VALUE *argv, VALUE obj,
               int sss)
{
  gsl_wavelet *w = NULL;
  gsl_vector *v = NULL, *vnew;
  gsl_wavelet_direction dir = gsl_wavelet_forward;
  gsl_wavelet_workspace *work = NULL;
  int itmp, flag = 0;
  // local variable "status" declared and set, but never used
  //int status;
  double *ptr1, *ptr2;
  size_t n, stride;
  int naflag = 0;
  VALUE ary, ret;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na1 = NULL;
#endif

  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 2) rb_raise(rb_eArgError, "too few arguments");
    CHECK_WAVELET(argv[0]);

    if (MATRIX_P(argv[1])) {
      return rb_gsl_wavelet2d(argc, argv, obj,
            gsl_wavelet2d_transform_matrix, sss);
    }
    if (VECTOR_P(argv[1])) {
      Data_Get_Struct(argv[0], gsl_wavelet, w);
      Data_Get_Struct(argv[1], gsl_vector, v);
      ret = argv[1];
      ptr1 = v->data;
      n = v->size;
      stride = v->stride;
#ifdef HAVE_NARRAY_H
    } else if (NA_IsNArray(argv[1])) {
      GetNArray(argv[1], na1);
      ret = argv[1];
      ptr1 = (double*) na1->ptr;
      n = na1->total;
      naflag = 1;
      stride = 1;
#endif
    } else {
      rb_raise(rb_eTypeError, "wrong argument type (Vector expected)");
    }
    itmp = 2;
    break;
  default:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments");

    if (MATRIX_P(argv[0])) {
      return rb_gsl_wavelet2d(argc, argv, obj,
            gsl_wavelet2d_transform_matrix, sss);
    }
    if (VECTOR_P(obj)) {
      CHECK_WAVELET(argv[0]);
      Data_Get_Struct(argv[0], gsl_wavelet, w);
      Data_Get_Struct(obj, gsl_vector, v);
      ret = obj;
      ptr1 = v->data;
      n = v->size;
      stride = v->stride;
    } else if (VECTOR_P(argv[0])) {

      CHECK_WAVELET(obj);
      Data_Get_Struct(obj, gsl_wavelet, w);
      Data_Get_Struct(argv[0], gsl_vector, v);
      ret = argv[0];
      ptr1 = v->data;
      n = v->size;
      stride = v->stride;
#ifdef HAVE_NARRAY_H
    } else if (NA_IsNArray(obj)) {
      CHECK_WAVELET(argv[0]);
      Data_Get_Struct(argv[0], gsl_wavelet, w);
      GetNArray(obj, na1);
      ret = obj;
      ptr1 = (double*) na1->ptr;
      n = na1->total;
      naflag = 1;
      stride = 1;
    } else if (NA_IsNArray(argv[0])) {
      CHECK_WAVELET(obj);
      Data_Get_Struct(obj, gsl_wavelet, w);
      GetNArray(argv[0], na1);
      ret = argv[0];
      ptr1 = (double*) na1->ptr;
      n = na1->total;
      naflag = 1;
      stride = 1;
#endif
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    itmp = 1;
    break;
  }
  switch (argc - itmp) {
  case 2:
    CHECK_FIXNUM(argv[itmp]);
    CHECK_WORKSPACE(argv[itmp+1]);
    dir = FIX2INT(argv[itmp]);
    Data_Get_Struct(argv[itmp+1], gsl_wavelet_workspace, work);
    break;
  case 1:
    if (TYPE(argv[itmp]) == T_FIXNUM) {
      dir = FIX2INT(argv[itmp]);
      work = gsl_wavelet_workspace_alloc(v->size);
      flag = 1;
    } else if (rb_obj_is_kind_of(argv[itmp], cgsl_wavelet_workspace)) {
      Data_Get_Struct(argv[itmp], gsl_wavelet_workspace, work);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  case 0:
    work = gsl_wavelet_workspace_alloc(v->size);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments");
    break;
  }
  if (naflag == 0) {
    if (sss == RB_GSL_DWT_COPY) { 
      vnew = gsl_vector_alloc(v->size);
      gsl_vector_memcpy(vnew, v);
      ary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
      ptr2 = vnew->data;
    } else {
      ary = ret;
      ptr2 = ptr1;
    }
  } else {
#ifdef HAVE_NARRAY_H
    if (sss == RB_GSL_DWT_COPY) {
      ary = na_make_object(NA_DFLOAT, na1->rank, na1->shape, cNArray);
      ptr2 = NA_PTR_TYPE(ary, double*);
      memcpy(ptr2, ptr1, sizeof(double)*n);
    } else {
      ary = ret;
      ptr2 = ptr1;
    }
#endif
  }
  /*status =*/ gsl_wavelet_transform(w, ptr2, stride, n, dir, work);
  if (flag) gsl_wavelet_workspace_free(work);
  return ary;
}

static VALUE rb_gsl_wavelet_transform(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet_transform0(argc, argv, obj, RB_GSL_DWT_COPY);
}

static VALUE rb_gsl_wavelet_transform2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet_transform0(argc, argv, obj, RB_GSL_DWT_INPLACE);
}

static VALUE rb_gsl_wavelet_trans(int argc, VALUE *argv, VALUE obj,
          int (*trans)(const gsl_wavelet *, 
                 double *, size_t, size_t, 
                 gsl_wavelet_workspace *),
          int sss)
{
  gsl_wavelet *w = NULL;
  gsl_vector *v = NULL, *vnew;
  gsl_wavelet_workspace *work = NULL;
  int itmp, flag = 0, naflag = 0;
  // local variable "status" declared and set, but never used
  //int status;
  double *ptr1 = NULL, *ptr2 = NULL;
  size_t n, stride;
  VALUE ary = Qnil, ret;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
#endif
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 2) rb_raise(rb_eArgError, "too few arguments");
    CHECK_WAVELET(argv[0]);

    if (MATRIX_P(argv[1])) {
      if (trans == gsl_wavelet_transform_forward) {
  return rb_gsl_wavelet2d_trans(argc, argv, obj, 
              gsl_wavelet2d_transform_matrix_forward, sss);
      } else {
  return rb_gsl_wavelet2d_trans(argc, argv, obj, 
              gsl_wavelet2d_transform_matrix_inverse, sss);
      }
    }
    if (VECTOR_P(argv[1])) {
      Data_Get_Struct(argv[0], gsl_wavelet, w);
      Data_Get_Struct(argv[1], gsl_vector, v);
      ret = argv[1];
      ptr1 = v->data;
      n = v->size;
      stride = v->stride;
#ifdef HAVE_NARRAY_H
    } else if (NA_IsNArray(argv[1])) {
      GetNArray(argv[1], na);
      ret = argv[1];
      ptr1 = (double*) na->ptr;
      n = na->total;
      naflag = 1;
      stride = 1;
#endif
    } else {
      rb_raise(rb_eTypeError, "wrong argument type (Vector expected)");
    }
    itmp = 2;
    break;
  default:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments");

    if (MATRIX_P(argv[0])) {
      if (trans == gsl_wavelet_transform_forward) {
  return rb_gsl_wavelet2d_trans(argc, argv, obj, 
              gsl_wavelet2d_transform_matrix_forward, sss);
      } else {
  return rb_gsl_wavelet2d_trans(argc, argv, obj, 
              gsl_wavelet2d_transform_matrix_inverse, sss);
      }
    }
    if (VECTOR_P(obj)) {
      CHECK_WAVELET(argv[0]);
      Data_Get_Struct(argv[0], gsl_wavelet, w);
      Data_Get_Struct(obj, gsl_vector, v);
      ret = obj;
      ptr1 = v->data;
      n = v->size;
      stride = v->stride;
    } else if (VECTOR_P(argv[0])) {
      CHECK_WAVELET(obj);
      Data_Get_Struct(obj, gsl_wavelet, w);
      Data_Get_Struct(argv[0], gsl_vector, v);
      ret = argv[0];
      ptr1 = v->data;
      n = v->size;
      stride = v->stride;
#ifdef HAVE_NARRAY_H
    } else if (NA_IsNArray(obj)) {
      CHECK_WAVELET(argv[0]);
      Data_Get_Struct(argv[0], gsl_wavelet, w);
      GetNArray(obj, na);
      ret = obj;
      ptr1 = (double*) na->ptr;
      n = na->total;
      naflag = 1;
      stride = 1;
    } else if (NA_IsNArray(argv[0])) {
      CHECK_WAVELET(obj);
      Data_Get_Struct(obj, gsl_wavelet, w);
      GetNArray(argv[0], na);
      ret = argv[0];
      ptr1 = (double*) na->ptr;
      n = na->total;
      naflag = 1;
      stride = 1;
#endif
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    itmp = 1;
    break;
  }
  switch (argc - itmp) {
  case 1:
    CHECK_WORKSPACE(argv[itmp]);
    Data_Get_Struct(argv[itmp], gsl_wavelet_workspace, work);
    break;
  case 0:
    work = gsl_wavelet_workspace_alloc(v->size);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments");
    break;
  }
  if (naflag == 0) {
    if (sss == RB_GSL_DWT_COPY) {
      vnew = gsl_vector_alloc(v->size);
      gsl_vector_memcpy(vnew, v);
      ary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
      ptr2 = vnew->data;
    } else {
      ptr2 = ptr1;
      ary = ret;
    }
  } else {
#ifdef HAVA_NARRAY_H
    if (sss == RB_GSL_DWT_COPY) {
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, cNArray);
      ptr2 = NA_PTR_TYPE(ary, double*);
      memcpy(ptr2, ptr1, sizeof(double)*n);
    } else {
      ptr2 = ptr1;
      ary = ret;
    }
#endif
  }
  /*status =*/ (*trans)(w, ptr2, stride, n, work);
  if (flag) gsl_wavelet_workspace_free(work);
  return ary;
}

static VALUE rb_gsl_wavelet_transform_forward(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet_trans(argc, argv, obj, gsl_wavelet_transform_forward,
            RB_GSL_DWT_COPY);
}

static VALUE rb_gsl_wavelet_transform_inverse(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet_trans(argc, argv, obj, gsl_wavelet_transform_inverse,
            RB_GSL_DWT_COPY);
}

static VALUE rb_gsl_wavelet_transform_forward2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet_trans(argc, argv, obj, gsl_wavelet_transform_forward,
            RB_GSL_DWT_INPLACE);
}

static VALUE rb_gsl_wavelet_transform_inverse2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet_trans(argc, argv, obj, gsl_wavelet_transform_inverse,
            RB_GSL_DWT_INPLACE);
}

static VALUE rb_gsl_wavelet2d(int argc, VALUE *argv, VALUE obj,
            int (*trans)(const gsl_wavelet *, 
             gsl_matrix *, 
             gsl_wavelet_direction, 
             gsl_wavelet_workspace *),
            int sss)
{
  gsl_wavelet *w = NULL;
  gsl_matrix *m = NULL, *mnew;
  gsl_wavelet_direction dir = gsl_wavelet_forward;
  gsl_wavelet_workspace *work = NULL;
  VALUE ary, ret;
  int itmp, flag = 0;
  // local variable "status" declared and set, but never used
  //int status;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 2) rb_raise(rb_eArgError, "too few arguments");
    CHECK_WAVELET(argv[0]);
    CHECK_MATRIX(argv[1]);
    ret = argv[1];
    Data_Get_Struct(argv[0], gsl_wavelet, w);
    Data_Get_Struct(argv[1], gsl_matrix, m);
    itmp = 2;
    break;
  default:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments");
    if (MATRIX_P(obj)) {
      CHECK_WAVELET(argv[0]);
      Data_Get_Struct(argv[0], gsl_wavelet, w);
      Data_Get_Struct(obj, gsl_matrix, m);
      ret = obj;
    } else {
      CHECK_MATRIX(argv[0]);
      Data_Get_Struct(obj, gsl_wavelet, w);
      Data_Get_Struct(argv[0], gsl_matrix, m);
      ret = argv[0];
    }
    itmp = 1;
    break;
  }
  switch (argc - itmp) {
  case 2:
    CHECK_FIXNUM(argv[itmp]);
    CHECK_WORKSPACE(argv[itmp+1]);
    dir = FIX2INT(argv[itmp]);
    Data_Get_Struct(argv[itmp+1], gsl_wavelet_workspace, work);
    break;
  case 1:
    if (TYPE(argv[itmp]) == T_FIXNUM) {
      dir = FIX2INT(argv[itmp]);
      work = gsl_wavelet_workspace_alloc(m->size1);
      flag = 1;
    } else if (rb_obj_is_kind_of(argv[itmp], cgsl_wavelet_workspace)) {
      Data_Get_Struct(argv[itmp], gsl_wavelet_workspace, work);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  case 0:
    work = gsl_wavelet_workspace_alloc(m->size1);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments");
    break;
  }
  if (sss == RB_GSL_DWT_COPY) {
    mnew = make_matrix_clone(m);
    ary = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
  } else {
    mnew = m;
    ary = ret;
  }
  /*status =*/ (*trans)(w, mnew, dir, work);
  if (flag) gsl_wavelet_workspace_free(work);
  return ary;
}

static VALUE rb_gsl_wavelet2d_transform_matrix(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet2d(argc, argv, obj, gsl_wavelet2d_transform_matrix,
        RB_GSL_DWT_COPY);
}

static VALUE rb_gsl_wavelet2d_transform_matrix2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet2d(argc, argv, obj, gsl_wavelet2d_transform_matrix,
        RB_GSL_DWT_INPLACE);
}

static VALUE rb_gsl_wavelet2d_trans(int argc, VALUE *argv, VALUE obj,
            int (*trans)(const gsl_wavelet *, 
             gsl_matrix *,
             gsl_wavelet_workspace *),
            int sss)
{
  gsl_wavelet *w = NULL;
  gsl_matrix *m = NULL, *mnew;
  gsl_wavelet_workspace *work = NULL;
  VALUE ary, ret;
  int itmp, flag = 0;
  // local variable "status" declared and set, but never used
  //int status;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 2) rb_raise(rb_eArgError, "too few arguments");
    CHECK_WAVELET(argv[0]);
    CHECK_MATRIX(argv[1]);
    Data_Get_Struct(argv[0], gsl_wavelet, w);
    Data_Get_Struct(argv[1], gsl_matrix, m);
    ret = argv[1];
    itmp = 2;
    break;
  default:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments");
    if (MATRIX_P(obj)) {
      CHECK_WAVELET(argv[0]);
      Data_Get_Struct(argv[0], gsl_wavelet, w);
      Data_Get_Struct(obj, gsl_matrix, m);
      ret = obj;
    } else {
      CHECK_MATRIX(argv[0]);
      Data_Get_Struct(obj, gsl_wavelet, w);
      Data_Get_Struct(argv[0], gsl_matrix, m);
      ret = argv[0];
    }
    itmp = 1;
    break;
  }
  switch (argc - itmp) {
  case 1:
    CHECK_WORKSPACE(argv[itmp]);
    Data_Get_Struct(argv[itmp], gsl_wavelet_workspace, work);
    break;
  case 0:
    work = gsl_wavelet_workspace_alloc(m->size1);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments");
    break;
  }
  if (sss == RB_GSL_DWT_COPY) {
    mnew = make_matrix_clone(m);
    ary = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
  } else {
    mnew = m;
    ary = ret;
  }
  /*status =*/ (*trans)(w, mnew, work);
  if (flag) gsl_wavelet_workspace_free(work);
  return ary;
}

static VALUE rb_gsl_wavelet2d_transform_matrix_forward(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet2d_trans(argc, argv, obj, 
        gsl_wavelet2d_transform_matrix_forward,
        RB_GSL_DWT_COPY);
}

static VALUE rb_gsl_wavelet2d_transform_matrix_forward2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet2d_trans(argc, argv, obj, 
        gsl_wavelet2d_transform_matrix_forward,
        RB_GSL_DWT_INPLACE);
}

static VALUE rb_gsl_wavelet2d_transform_matrix_inverse(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet2d_trans(argc, argv, obj, 
        gsl_wavelet2d_transform_matrix_inverse,
        RB_GSL_DWT_COPY);
}

static VALUE rb_gsl_wavelet2d_transform_matrix_inverse2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet2d_trans(argc, argv, obj, 
        gsl_wavelet2d_transform_matrix_inverse,
        RB_GSL_DWT_INPLACE);
}

/** nstransform **/
static VALUE rb_gsl_wavelet2d_nstransform_matrix(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet2d(argc, argv, obj, gsl_wavelet2d_nstransform_matrix,
        RB_GSL_DWT_COPY);
}

static VALUE rb_gsl_wavelet2d_nstransform_matrix2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet2d(argc, argv, obj, gsl_wavelet2d_nstransform_matrix,
        RB_GSL_DWT_INPLACE);
}

static VALUE rb_gsl_wavelet2d_nstransform_matrix_forward(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet2d_trans(argc, argv, obj, 
        gsl_wavelet2d_nstransform_matrix_forward,
        RB_GSL_DWT_COPY);
}

static VALUE rb_gsl_wavelet2d_nstransform_matrix_forward2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet2d_trans(argc, argv, obj, 
        gsl_wavelet2d_nstransform_matrix_forward,
        RB_GSL_DWT_INPLACE);
}

static VALUE rb_gsl_wavelet2d_nstransform_matrix_inverse(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet2d_trans(argc, argv, obj, 
        gsl_wavelet2d_nstransform_matrix_inverse,
        RB_GSL_DWT_COPY);
}

static VALUE rb_gsl_wavelet2d_nstransform_matrix_inverse2(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_wavelet2d_trans(argc, argv, obj, 
        gsl_wavelet2d_nstransform_matrix_inverse,
        RB_GSL_DWT_INPLACE);
}

#endif

void Init_wavelet(VALUE module)
{
  VALUE cgsl_wavelet2d;

  cgsl_wavelet = rb_define_class_under(module, "Wavelet", cGSL_Object);
  cgsl_wavelet2d = rb_define_class_under(module, "Wavelet2d", cgsl_wavelet);

  rb_define_singleton_method(cgsl_wavelet, "alloc", rb_gsl_wavelet_new, 2);

#ifdef GSL_1_6_LATER
  rb_gsl_wavelet_define_const(cgsl_wavelet);
  rb_define_method(cgsl_wavelet, "name", rb_gsl_wavelet_name, 0);

  cgsl_wavelet_workspace = rb_define_class_under(cgsl_wavelet, "Workspace", cGSL_Object);
  rb_define_singleton_method(cgsl_wavelet_workspace, "alloc", 
           rb_gsl_wavelet_workspace_new, 1);

  /*****/

  rb_define_singleton_method(cgsl_wavelet, "transform", 
           rb_gsl_wavelet_transform, -1);
  rb_define_method(cgsl_wavelet, "transform", rb_gsl_wavelet_transform, -1);
  rb_define_method(cgsl_vector, "wavelet_transform", rb_gsl_wavelet_transform, -1);
  rb_define_singleton_method(cgsl_wavelet, "transform!", 
           rb_gsl_wavelet_transform2, -1);
  rb_define_method(cgsl_wavelet, "transform!", rb_gsl_wavelet_transform2, -1);
  rb_define_method(cgsl_vector, "wavelet_transform!", rb_gsl_wavelet_transform2, -1);

  /**/

  rb_define_singleton_method(cgsl_wavelet, "transform_forward", 
           rb_gsl_wavelet_transform_forward, -1);
  rb_define_method(cgsl_wavelet, "transform_forward", 
       rb_gsl_wavelet_transform_forward, -1);
  rb_define_alias(cgsl_wavelet, "forward", "transform_forward");
  rb_define_method(cgsl_vector, "wavelet_transform_forward", 
       rb_gsl_wavelet_transform_forward, -1);

  rb_define_singleton_method(cgsl_wavelet, "transform_inverse", 
           rb_gsl_wavelet_transform_inverse, -1);
  rb_define_method(cgsl_wavelet, "transform_inverse", 
       rb_gsl_wavelet_transform_inverse, -1);
  rb_define_alias(cgsl_wavelet, "inverse", "transform_inverse");
  rb_define_method(cgsl_vector, "wavelet_transform_inverse", 
       rb_gsl_wavelet_transform_inverse, -1);

  rb_define_singleton_method(cgsl_wavelet, "transform_forward!", 
           rb_gsl_wavelet_transform_forward2, -1);
  rb_define_method(cgsl_wavelet, "transform_forward!", 
       rb_gsl_wavelet_transform_forward2, -1);
  rb_define_alias(cgsl_wavelet, "forward!", "transform_forward!");
  rb_define_method(cgsl_vector, "wavelet_transform_forward!", 
       rb_gsl_wavelet_transform_forward2, -1);

  rb_define_singleton_method(cgsl_wavelet, "transform_inverse!", 
           rb_gsl_wavelet_transform_inverse2, -1);
  rb_define_method(cgsl_wavelet, "transform_inverse!", 
       rb_gsl_wavelet_transform_inverse2, -1);
  rb_define_alias(cgsl_wavelet, "inverse!", "transform_inverse!");
  rb_define_method(cgsl_vector, "wavelet_transform_inverse!", 
       rb_gsl_wavelet_transform_inverse2, -1);
  /***** 2d *****/
  rb_define_singleton_method(cgsl_wavelet, "transform_matrix", 
           rb_gsl_wavelet2d_transform_matrix, -1);
  rb_define_singleton_method(cgsl_wavelet2d, "transform", 
           rb_gsl_wavelet2d_transform_matrix, -1);
  rb_define_method(cgsl_wavelet, "transform_matrix",
       rb_gsl_wavelet2d_transform_matrix, -1);
  rb_define_method(cgsl_wavelet2d, "transform", 
       rb_gsl_wavelet2d_transform_matrix, -1);
  rb_define_method(cgsl_matrix, "wavelet_transform", 
       rb_gsl_wavelet2d_transform_matrix, -1);

  rb_define_singleton_method(cgsl_wavelet, "transform_matrix!", 
           rb_gsl_wavelet2d_transform_matrix2, -1);
  rb_define_singleton_method(cgsl_wavelet2d, "transform!", 
           rb_gsl_wavelet2d_transform_matrix2, -1);
  rb_define_method(cgsl_wavelet, "transform_matrix!",
       rb_gsl_wavelet2d_transform_matrix2, -1);
  rb_define_method(cgsl_wavelet2d, "transform!", 
       rb_gsl_wavelet2d_transform_matrix2, -1);
  rb_define_method(cgsl_matrix, "wavelet_transform!", 
       rb_gsl_wavelet2d_transform_matrix2, -1);
  /**/

  rb_define_singleton_method(cgsl_wavelet, "transform_matrix_forward", 
           rb_gsl_wavelet2d_transform_matrix_forward, -1);
  rb_define_singleton_method(cgsl_wavelet2d, "transform_forward", 
           rb_gsl_wavelet2d_transform_matrix_forward, -1);
  rb_define_method(cgsl_wavelet, "transform_matrix_forward", 
       rb_gsl_wavelet2d_transform_matrix_forward, -1);
  rb_define_method(cgsl_wavelet2d, "transform_forward", 
       rb_gsl_wavelet2d_transform_matrix_forward, -1);
  rb_define_method(cgsl_matrix, "wavelet_transform_forward", 
       rb_gsl_wavelet2d_transform_matrix_forward, -1);

  rb_define_singleton_method(cgsl_wavelet, "transform_matrix_forward!", 
           rb_gsl_wavelet2d_transform_matrix_forward2, -1);
  rb_define_singleton_method(cgsl_wavelet2d, "transform_forward!", 
           rb_gsl_wavelet2d_transform_matrix_forward2, -1);
  rb_define_method(cgsl_wavelet, "transform_matrix_forward!", 
       rb_gsl_wavelet2d_transform_matrix_forward2, -1);
  rb_define_method(cgsl_wavelet2d, "transform_forward!", 
       rb_gsl_wavelet2d_transform_matrix_forward2, -1);
  rb_define_method(cgsl_matrix, "wavelet_transform_forward!", 
       rb_gsl_wavelet2d_transform_matrix_forward2, -1);

  rb_define_singleton_method(cgsl_wavelet, "transform_matrix_inverse", 
           rb_gsl_wavelet2d_transform_matrix_inverse, -1);
  rb_define_singleton_method(cgsl_wavelet2d, "transform_inverse", 
           rb_gsl_wavelet2d_transform_matrix_inverse, -1);
  rb_define_method(cgsl_wavelet, "transform_matrix_inverse", 
       rb_gsl_wavelet2d_transform_matrix_inverse, -1);
  rb_define_method(cgsl_wavelet2d, "transform_inverse", 
       rb_gsl_wavelet2d_transform_matrix_inverse, -1);
  rb_define_method(cgsl_matrix, "wavelet_transform_inverse", 
       rb_gsl_wavelet2d_transform_matrix_inverse, -1);

  rb_define_singleton_method(cgsl_wavelet, "transform_matrix_inverse!", 
           rb_gsl_wavelet2d_transform_matrix_inverse2, -1);
  rb_define_singleton_method(cgsl_wavelet2d, "transform_inverse!", 
           rb_gsl_wavelet2d_transform_matrix_inverse2, -1);
  rb_define_method(cgsl_wavelet, "transform_matrix_inverse!", 
       rb_gsl_wavelet2d_transform_matrix_inverse2, -1);
  rb_define_method(cgsl_wavelet2d, "transform_inverse!", 
       rb_gsl_wavelet2d_transform_matrix_inverse2, -1);
  rb_define_method(cgsl_matrix, "wavelet_transform_inverse!", 
       rb_gsl_wavelet2d_transform_matrix_inverse2, -1);

  /** nstransform **/
  rb_define_singleton_method(cgsl_wavelet, "nstransform_matrix", 
           rb_gsl_wavelet2d_nstransform_matrix, -1);
  rb_define_singleton_method(cgsl_wavelet2d, "nstransform", 
           rb_gsl_wavelet2d_nstransform_matrix, -1);
  rb_define_method(cgsl_wavelet, "nstransform_matrix", 
       rb_gsl_wavelet2d_nstransform_matrix, -1);
  rb_define_method(cgsl_wavelet2d, "nstransform", 
       rb_gsl_wavelet2d_nstransform_matrix, -1);

  rb_define_singleton_method(cgsl_wavelet, "nstransform_matrix!", 
           rb_gsl_wavelet2d_nstransform_matrix2, -1);
  rb_define_singleton_method(cgsl_wavelet2d, "nstransform!", 
           rb_gsl_wavelet2d_nstransform_matrix2, -1);
  rb_define_method(cgsl_wavelet, "nstransform_matrix!", 
       rb_gsl_wavelet2d_nstransform_matrix2, -1);
  rb_define_method(cgsl_wavelet2d, "nstransform!", 
       rb_gsl_wavelet2d_nstransform_matrix2, -1);
  /**/

  rb_define_singleton_method(cgsl_wavelet, "nstransform_matrix_forward", 
           rb_gsl_wavelet2d_nstransform_matrix_forward, -1);
  rb_define_singleton_method(cgsl_wavelet2d, "nstransform_forward", 
           rb_gsl_wavelet2d_nstransform_matrix_forward, -1);
  rb_define_method(cgsl_wavelet, "nstransform_matrix_forward", 
       rb_gsl_wavelet2d_nstransform_matrix_forward, -1);
  rb_define_method(cgsl_wavelet2d, "nstransform_forward", 
       rb_gsl_wavelet2d_nstransform_matrix_forward, -1);

  rb_define_singleton_method(cgsl_wavelet, "nstransform_matrix_forward!", 
           rb_gsl_wavelet2d_nstransform_matrix_forward2, -1);
  rb_define_singleton_method(cgsl_wavelet2d, "nstransform_forward!", 
           rb_gsl_wavelet2d_nstransform_matrix_forward2, -1);
  rb_define_method(cgsl_wavelet, "nstransform_matrix_forward!", 
       rb_gsl_wavelet2d_nstransform_matrix_forward2, -1);
  rb_define_method(cgsl_wavelet2d, "nstransform_forward!", 
       rb_gsl_wavelet2d_nstransform_matrix_forward2, -1);

  rb_define_singleton_method(cgsl_wavelet, "nstransform_matrix_inverse", 
           rb_gsl_wavelet2d_nstransform_matrix_inverse, -1);
  rb_define_singleton_method(cgsl_wavelet2d, "nstransform_inverse", 
           rb_gsl_wavelet2d_nstransform_matrix_inverse, -1);
  rb_define_method(cgsl_wavelet, "nstransform_matrix_inverse", 
       rb_gsl_wavelet2d_nstransform_matrix_inverse, -1);
  rb_define_method(cgsl_wavelet2d, "nstransform_inverse", 
       rb_gsl_wavelet2d_nstransform_matrix_inverse, -1);

  rb_define_singleton_method(cgsl_wavelet, "nstransform_matrix_inverse!", 
           rb_gsl_wavelet2d_nstransform_matrix_inverse2, -1);
  rb_define_singleton_method(cgsl_wavelet2d, "nstransform_inverse!", 
           rb_gsl_wavelet2d_nstransform_matrix_inverse2, -1);
  rb_define_method(cgsl_wavelet, "nstransform_matrix_inverse!", 
       rb_gsl_wavelet2d_nstransform_matrix_inverse2, -1);
  rb_define_method(cgsl_wavelet2d, "nstransform_inverse!", 
       rb_gsl_wavelet2d_nstransform_matrix_inverse2, -1);

#endif
}

#undef CHECK_WAVELET
#undef CHECK_WORKSPACE

