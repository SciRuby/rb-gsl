/*
  fft.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_fft.h"

VALUE mgsl_fft;
VALUE cgsl_fft_wavetable;
VALUE cgsl_fft_complex_wavetable, cgsl_fft_complex_workspace;
VALUE cgsl_fft_real_wavetable, cgsl_fft_halfcomplex_wavetable;
VALUE cgsl_fft_real_workspace;

extern VALUE cgsl_vector_int;
extern VALUE cgsl_vector_complex;

static void GSL_FFT_Workspace_free(GSL_FFT_Workspace *space);

static VALUE rb_gsl_fft_complex_wavetable_new(VALUE klass, VALUE n)
{
  CHECK_FIXNUM(n);
  return Data_Wrap_Struct(cgsl_fft_complex_wavetable, 0,
        gsl_fft_complex_wavetable_free,
        gsl_fft_complex_wavetable_alloc(FIX2INT(n)));
}

static VALUE rb_gsl_fft_real_wavetable_new(VALUE klass, VALUE n)
{
  CHECK_FIXNUM(n);
  return Data_Wrap_Struct(klass, 0, gsl_fft_real_wavetable_free,
        gsl_fft_real_wavetable_alloc(FIX2INT(n)));
}

static VALUE rb_gsl_fft_halfcomplex_wavetable_new(VALUE klass, VALUE n)
{
  CHECK_FIXNUM(n);
  return Data_Wrap_Struct(klass, 0, gsl_fft_halfcomplex_wavetable_free,
        gsl_fft_halfcomplex_wavetable_alloc(FIX2INT(n)));

}

static void GSL_FFT_Wavetable_free(GSL_FFT_Wavetable *table)
{
  gsl_fft_complex_wavetable_free((gsl_fft_complex_wavetable *) table);
}

static VALUE rb_gsl_fft_complex_workspace_new(VALUE klass, VALUE n)
{
  CHECK_FIXNUM(n);
  return Data_Wrap_Struct(klass, 0, gsl_fft_complex_workspace_free,
        gsl_fft_complex_workspace_alloc(FIX2INT(n)));
}

static VALUE rb_gsl_fft_real_workspace_new(VALUE klass, VALUE n)
{
  CHECK_FIXNUM(n);
  return Data_Wrap_Struct(klass, 0, gsl_fft_real_workspace_free,
        gsl_fft_real_workspace_alloc(FIX2INT(n)));
}

static void GSL_FFT_Workspace_free(GSL_FFT_Workspace *space)
{
  gsl_fft_complex_workspace_free((gsl_fft_complex_workspace *) space);
}

// The FFT methods used to allow passing stride and n values as optional
// parameters to control which elements get transformed.  This created problems
// for Views which can have their own stride, so support for stride and n
// parameters to the transform methods is being dropped.  This method used to
// be called to determine the stride and n values to use based on the
// parameters and/or the vector itself (depending on how many parameters were
// passed). Now this function is somewhat unneceesary, but to simplify the code
// refactoring, it has been left in place for the time being.  Eventually it
// can be refactored away completely.
static VALUE get_complex_stride_n(VALUE obj,
             gsl_vector_complex **vin,
             gsl_complex_packed_array *data, size_t *stride, size_t *n)
{
  gsl_vector_complex *v = NULL;

  // obj must be a GSL::Vector::Complex
  CHECK_VECTOR_COMPLEX(obj);
  Data_Get_Struct(obj, gsl_vector_complex, v);

  if(vin) *vin = v;
  *data = (gsl_complex_packed_array) v->data;
  *stride = v->stride;
  *n = v->size;
  return obj;
}

static VALUE rb_fft_complex_radix2(VALUE obj,
           int (*trans)(gsl_complex_packed_array,
            size_t, size_t), int flag)
{
  size_t stride, n;
  gsl_complex_packed_array data;
  gsl_vector_complex *vin, *vout;
  VALUE ary;
  ary = get_complex_stride_n(obj, &vin, &data, &stride, &n);
  if (flag == RB_GSL_FFT_COPY) {
    vout = gsl_vector_complex_alloc(n);
    gsl_vector_complex_memcpy(vout, vin);
    (*trans)(vout->data, vout->stride /*1*/, vout->size /*n*/);
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vout);
  } else { /* in-place */
    (*trans)(data, stride, n);
    return ary;
  }
}

static VALUE rb_gsl_fft_complex_radix2_forward(VALUE obj)
{
  return rb_fft_complex_radix2(obj, gsl_fft_complex_radix2_forward,
             RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_radix2_forward2(VALUE obj)
{
  return rb_fft_complex_radix2(obj, gsl_fft_complex_radix2_forward,
             RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_complex_radix2_transform(VALUE obj, VALUE val_sign)
{
  size_t stride, n;
  gsl_complex_packed_array data;
  gsl_fft_direction sign;
  gsl_vector_complex *vin, *vout;
  sign = NUM2INT(val_sign);
  get_complex_stride_n(obj, &vin, &data, &stride, &n);
  vout = gsl_vector_complex_alloc(n);
  gsl_vector_complex_memcpy(vout, vin);
  gsl_fft_complex_radix2_transform(vout->data, vout->stride /*1*/, vout->size /*n*/, sign);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vout);
}

static VALUE rb_gsl_fft_complex_radix2_transform2(VALUE obj, VALUE val_sign)
{
  size_t stride, n;
  gsl_complex_packed_array data;
  gsl_fft_direction sign;
  VALUE ary;
  sign = NUM2INT(val_sign);
  ary = get_complex_stride_n(obj, NULL, &data, &stride, &n);
  gsl_fft_complex_radix2_transform(data, stride, n, sign);
  return ary;
}

static VALUE rb_gsl_fft_complex_radix2_backward(VALUE obj)
{
  return rb_fft_complex_radix2(obj, gsl_fft_complex_radix2_backward,
             RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_radix2_inverse(VALUE obj)
{
  return rb_fft_complex_radix2(obj, gsl_fft_complex_radix2_inverse,
             RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_radix2_dif_forward(VALUE obj)
{
  return rb_fft_complex_radix2(obj,
             gsl_fft_complex_radix2_dif_forward,
             RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_radix2_backward2(VALUE obj)
{
  return rb_fft_complex_radix2(obj, gsl_fft_complex_radix2_backward,
             RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_complex_radix2_inverse2(VALUE obj)
{
  return rb_fft_complex_radix2(obj, gsl_fft_complex_radix2_inverse,
             RB_GSL_FFT_INPLACE);
}


static VALUE rb_gsl_fft_complex_radix2_dif_forward2(VALUE obj)
{
  return rb_fft_complex_radix2(obj,
             gsl_fft_complex_radix2_dif_forward,
             RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_complex_radix2_dif_transform(VALUE obj, VALUE val_sign)
{
  size_t stride, n;
  gsl_complex_packed_array data;
  gsl_vector_complex *vin, *vout;
  gsl_fft_direction sign;
  sign = NUM2INT(val_sign);
  get_complex_stride_n(obj, &vin, &data, &stride, &n);
  vout = gsl_vector_complex_alloc(n);
  gsl_vector_complex_memcpy(vout, vin);
  gsl_fft_complex_radix2_dif_transform(vout->data, vout->stride /*1*/, vout->size /*n*/, sign);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vout);
}

/* in-place */
static VALUE rb_gsl_fft_complex_radix2_dif_transform2(VALUE obj, VALUE val_sign)
{
  size_t stride, n;
  gsl_complex_packed_array data;
  gsl_fft_direction sign;
  VALUE ary;
  sign = NUM2INT(val_sign);
  ary = get_complex_stride_n(obj, NULL, &data, &stride, &n);
  gsl_fft_complex_radix2_dif_transform(data, stride, n, sign);
  return ary;
}

static VALUE rb_gsl_fft_complex_radix2_dif_backward(VALUE obj)
{
  return rb_fft_complex_radix2(obj,
             gsl_fft_complex_radix2_dif_backward,
             RB_GSL_FFT_COPY);
}
static VALUE rb_gsl_fft_complex_radix2_dif_inverse(VALUE obj)
{
  return rb_fft_complex_radix2(obj,
             gsl_fft_complex_radix2_dif_inverse,
             RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_radix2_dif_backward2(VALUE obj)
{
  return rb_fft_complex_radix2(obj,
             gsl_fft_complex_radix2_dif_backward,
             RB_GSL_FFT_INPLACE);
}
static VALUE rb_gsl_fft_complex_radix2_dif_inverse2(VALUE obj)
{
  return rb_fft_complex_radix2(obj,
             gsl_fft_complex_radix2_dif_inverse,
             RB_GSL_FFT_INPLACE);
}

static VALUE rb_GSL_FFT_Wavetable_n(VALUE obj)
{
  GSL_FFT_Wavetable *table;
  Data_Get_Struct(obj, GSL_FFT_Wavetable, table);
  return INT2FIX(table->n);
}

static VALUE rb_GSL_FFT_Wavetable_nf(VALUE obj)
{
  GSL_FFT_Wavetable *table;
  Data_Get_Struct(obj, GSL_FFT_Wavetable, table);
  return INT2FIX(table->nf);
}

static VALUE rb_GSL_FFT_Wavetable_factor(VALUE obj)
{
  GSL_FFT_Wavetable *table;
  gsl_vector_int *v;
  size_t i;
  Data_Get_Struct(obj, GSL_FFT_Wavetable, table);
  v = gsl_vector_int_alloc(table->nf);
  for (i = 0; i < table->nf; i++) gsl_vector_int_set(v, i, table->factor[i]);
  return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, v);
}

enum {
  NONE_OF_TWO = 0,
  ALLOC_SPACE = 1,
  ALLOC_TABLE = 2,
  BOTH_OF_TWO = 3,
} FFTComplexStructAllocFlag;

static void gsl_fft_free(int flag, GSL_FFT_Wavetable *table,
       GSL_FFT_Workspace *space);

// Parse argc, argv.  obj must be GSL::Vector::Complex.
// This can be simplified at some point.
// See comments preceding get_complex_stride_n()
static int gsl_fft_get_argv_complex(int argc, VALUE *argv, VALUE obj,
          gsl_vector_complex ** vin,
          gsl_complex_packed_array *data, size_t *stride,
          size_t *n, gsl_fft_complex_wavetable **table,
          gsl_fft_complex_workspace **space)
{
  int flag = NONE_OF_TWO, flagtmp, i, itmp = argc, itmp2 = 0, ccc;
  int flagw = 0;

  CHECK_VECTOR_COMPLEX(obj);

  ccc = argc;
  flagtmp = 0;
  flagw = 0;
  for (i = argc-1; i >= itmp2; i--) {
    if (rb_obj_is_kind_of(argv[i], cgsl_fft_complex_workspace)) {
      Data_Get_Struct(argv[i], gsl_fft_complex_workspace, *space);
      flagtmp = 1;
      flagw = 1;
      itmp = i;
      ccc--;
      break;
    }
  }
  flagtmp = 0;
  for (i = itmp-1; i >= itmp2; i--) {
    if (rb_obj_is_kind_of(argv[i], cgsl_fft_complex_wavetable)) {
      Data_Get_Struct(argv[i], gsl_fft_complex_wavetable, *table);
      flagtmp = 1;
      ccc--;
      break;
    }
  }
  get_complex_stride_n(obj, vin, data, stride, n);
  if (flagw == 0) {
    *space = gsl_fft_complex_workspace_alloc(*n);
    flag += ALLOC_SPACE;
  }
  if (flagtmp == 0) {
    *table = gsl_fft_complex_wavetable_alloc(*n);
    flag += ALLOC_TABLE;
  }
  if (*table == NULL) {
    rb_raise(rb_eRuntimeError, "something wrong with wavetable");
  }
  if (*space == NULL) {
    rb_raise(rb_eRuntimeError, "something wrong with workspace");
  }
  return flag;
}

// Parse argc, argv.  obj must be GSL::Vector of real data
static int gsl_fft_get_argv_real(int argc, VALUE *argv, VALUE obj,
           double **ptr, size_t *stride,
           size_t *n, gsl_fft_real_wavetable **table,
           gsl_fft_real_workspace **space, int *naflag)
{
  int flag = NONE_OF_TWO, flagtmp, i, itmp = argc, itmp2 = 0, ccc;
  int flagw = 0;
  *naflag = 0;

  *ptr = get_ptr_double3(obj, n, stride, naflag);

  ccc = argc;
  flagtmp = 0;
  flagw = 0;
  for (i = argc-1; i >= itmp2; i--) {
    if (rb_obj_is_kind_of(argv[i], cgsl_fft_real_workspace)) {
      Data_Get_Struct(argv[i], gsl_fft_real_workspace, *space);
      flagtmp = 1;
      flagw = 1;
      itmp = i;
      ccc--;
      break;
    }
  }
  flagtmp = 0;
  for (i = itmp-1; i >= itmp2; i--) {
    if (rb_obj_is_kind_of(argv[i], cgsl_fft_real_wavetable)) {
      Data_Get_Struct(argv[i], gsl_fft_real_wavetable, *table);
      flagtmp = 1;
      ccc--;
      break;
    }
  }
  if (flagw == 0) {
    *space = gsl_fft_real_workspace_alloc(*n);
    flag += ALLOC_SPACE;
  }
  if (flagtmp == 0) {
    *table = gsl_fft_real_wavetable_alloc(*n);
    flag += ALLOC_TABLE;
  }
  if (*table == NULL) {
    rb_raise(rb_eRuntimeError, "something wrong with wavetable");
  }
  if (*space == NULL) {
    rb_raise(rb_eRuntimeError, "something wrong with workspace");
  }
  return flag;
}

// Parse argc, argv.  obj must be GSL::Vector of halfcomplex data
static int gsl_fft_get_argv_halfcomplex(int argc, VALUE *argv, VALUE obj,
         double **ptr, size_t *stride,
          size_t *n, gsl_fft_halfcomplex_wavetable **table,
          gsl_fft_real_workspace **space, int *naflag)
{
  int flag = NONE_OF_TWO, flagtmp, i, itmp = argc, itmp2 = 0, ccc;
  int flagw = 0;

  *ptr = get_ptr_double3(obj, n, stride, naflag);

  ccc = argc;
  flagtmp = 0;
  flagw = 0;
  for (i = argc-1; i >= itmp2; i--) {
    if (rb_obj_is_kind_of(argv[i], cgsl_fft_real_workspace)) {
      Data_Get_Struct(argv[i], gsl_fft_real_workspace, *space);
      flagtmp = 1;
      flagw = 1;
      itmp = i;
      ccc--;
      break;
    }
  }
  flagtmp = 0;
  for (i = itmp-1; i >= itmp2; i--) {
    if (rb_obj_is_kind_of(argv[i], cgsl_fft_halfcomplex_wavetable)) {
      Data_Get_Struct(argv[i], gsl_fft_halfcomplex_wavetable, *table);
      flagtmp = 1;
      ccc--;
      break;
    }
  }
  if (flagw == 0) {
    *space = gsl_fft_real_workspace_alloc(*n);
    flag += ALLOC_SPACE;
  }
  if (flagtmp == 0) {
    *table = gsl_fft_halfcomplex_wavetable_alloc(*n);
    flag += ALLOC_TABLE;
  }
  if (*table == NULL) {
    rb_raise(rb_eRuntimeError, "something wrong with wavetable");
  }
  if (*space == NULL) {
    rb_raise(rb_eRuntimeError, "something wrong with workspace");
  }
  return flag;
}

static void gsl_fft_free(int flag, GSL_FFT_Wavetable *table,
       GSL_FFT_Workspace *space)
{
  switch (flag) {
  case ALLOC_TABLE:
    GSL_FFT_Wavetable_free(table);
    break;
  case ALLOC_SPACE:
    GSL_FFT_Workspace_free(space);
    break;
  case BOTH_OF_TWO:
    GSL_FFT_Wavetable_free(table);
    GSL_FFT_Workspace_free(space);
    break;
  default:
    /* never happens */
    break;
  }
}

static VALUE rb_fft_complex_trans(int argc, VALUE *argv, VALUE obj,
          int (*transform)(gsl_complex_packed_array,
               size_t, size_t,
               const gsl_fft_complex_wavetable *,
               gsl_fft_complex_workspace *),
          int sss)
{
  int flag = 0;
  // local variable "status" was defined and set, but never used
  //int status;
  size_t stride, n;
  gsl_complex_packed_array data;
  gsl_vector_complex *vin, *vout;
  gsl_fft_complex_wavetable *table = NULL;
  gsl_fft_complex_workspace *space = NULL;
  flag = gsl_fft_get_argv_complex(argc, argv, obj, &vin, &data, &stride, &n, &table, &space);
  if (sss == RB_GSL_FFT_COPY) {
    vout = gsl_vector_complex_alloc(n);
    gsl_vector_complex_memcpy(vout, vin);
    /*status =*/ (*transform)(vout->data, vout->stride /*1*/, vout->size /*n*/, table, space);
    gsl_fft_free(flag, (GSL_FFT_Wavetable *) table, (GSL_FFT_Workspace *) space);
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vout);
  } else {    /* in-place */
    /*status =*/ (*transform)(data, stride, n, table, space);
    gsl_fft_free(flag, (GSL_FFT_Wavetable *) table, (GSL_FFT_Workspace *) space);
    return obj;
  }
}

static VALUE rb_gsl_fft_complex_forward(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_trans(argc, argv, obj, gsl_fft_complex_forward,
            RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_forward2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_trans(argc, argv, obj, gsl_fft_complex_forward,
            RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_complex_transform(int argc, VALUE *argv, VALUE obj)
{
  int flag = 0;
  // local variable "status" was defined and set, but never used
  //int status;
  size_t stride, n;
  gsl_vector_complex *vin, *vout;
  gsl_fft_direction sign;
  gsl_complex_packed_array data;
  gsl_fft_complex_wavetable *table = NULL;
  gsl_fft_complex_workspace *space = NULL;
  CHECK_FIXNUM(argv[argc-1]);
  sign = FIX2INT(argv[argc-1]);
  flag = gsl_fft_get_argv_complex(argc-1, argv, obj, &vin, &data, &stride, &n, &table, &space);
  vout = gsl_vector_complex_alloc(n);
  gsl_vector_complex_memcpy(vout, vin);
  /*status =*/ gsl_fft_complex_transform(vout->data, stride, n, table, space, sign);
  gsl_fft_free(flag, (GSL_FFT_Wavetable *) table, (GSL_FFT_Workspace *) space);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vout);
}

/* in-place */
static VALUE rb_gsl_fft_complex_transform2(int argc, VALUE *argv, VALUE obj)
{
  int flag = 0;
  // local variable "status" was defined and set, but never used
  //int status;
  size_t stride, n;
  gsl_fft_direction sign;
  gsl_complex_packed_array data;
  gsl_fft_complex_wavetable *table = NULL;
  gsl_fft_complex_workspace *space = NULL;
  CHECK_FIXNUM(argv[argc-1]);
  sign = FIX2INT(argv[argc-1]);
  flag = gsl_fft_get_argv_complex(argc-1, argv, obj, NULL, &data, &stride, &n, &table, &space);
  /*status =*/ gsl_fft_complex_transform(data, stride, n, table, space, sign);
  gsl_fft_free(flag, (GSL_FFT_Wavetable *) table, (GSL_FFT_Workspace *) space);
  return obj;
}

static VALUE rb_gsl_fft_complex_backward(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_trans(argc, argv, obj, gsl_fft_complex_backward,
            RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_backward2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_trans(argc, argv, obj, gsl_fft_complex_backward,
            RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_complex_inverse(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_trans(argc, argv, obj, gsl_fft_complex_inverse,
            RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_inverse2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_trans(argc, argv, obj, gsl_fft_complex_inverse,
            RB_GSL_FFT_INPLACE);
}

// The FFT methods used to allow passing stride and n values as optional
// parameters to control which elements get transformed.  This created problems
// for Views which can have their own stride, so support for stride and n
// parameters to the transform methods is being dropped.  This method used to
// be called to determine the stride and n values to use based on the
// parameters and/or the vector itself (depending on how many parameters were
// passed). Now this function is somewhat unneceesary, but to simplify the code
// refactoring, it has been left in place for the time being.  Eventually it
// can be refactored away completely.
//
// obj must be GSL::Vector of real or halfcomplex data
static VALUE get_ptr_stride_n(VALUE obj,
           double **ptr, size_t *stride, size_t *n, int *flag)
{
  *flag = 0;
  *ptr = get_ptr_double3(obj, n, stride, flag);
  return obj;
}

static VALUE rb_fft_radix2(VALUE obj,
         int (*trans)(double [], size_t, size_t),
         int sss)
{
  size_t stride, n;
  gsl_vector *vnew;
  gsl_vector_view vv;
  double *ptr1, *ptr2;
  int flag;
#ifdef HAVE_NARRAY_H
  int shape[1];
#endif
  VALUE ary;
  get_ptr_stride_n(obj, &ptr1, &stride, &n, &flag);
  if (flag == 0) {
    if (sss == RB_GSL_FFT_COPY) {
      vnew = gsl_vector_alloc(n);
      vv.vector.data = ptr1;
      vv.vector.stride = stride;
      vv.vector.size = n;
      gsl_vector_memcpy(vnew, &vv.vector);
      ptr2 = vnew->data;
      stride = 1;
      ary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else {
      ary = obj;
      ptr2 = ptr1;
    }
#ifdef HAVE_NARRAY_H
  } else if (flag == 1) {
    if (sss == RB_GSL_FFT_COPY) {
      shape[0] = n;
      ary = na_make_object(NA_DFLOAT, 1, shape, cNArray);
      ptr2 = NA_PTR_TYPE(ary, double*);
      memcpy(ptr2, ptr1, sizeof(double)*n);
      stride = 1;
    } else {
      ary = obj;
      ptr2 = NA_PTR_TYPE(ary, double*);
    }
#endif
  } else {
    rb_raise(rb_eRuntimeError, "something wrong");
  }
  (*trans)(ptr2, stride, n);
  return ary;
}

static VALUE rb_gsl_fft_real_radix2_transform(VALUE obj)
{
  return rb_fft_radix2(obj, gsl_fft_real_radix2_transform,
           RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_real_radix2_transform2(VALUE obj)
{
  return rb_fft_radix2(obj, gsl_fft_real_radix2_transform,
           RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_halfcomplex_radix2_inverse(VALUE obj)
{
  return rb_fft_radix2(obj, gsl_fft_halfcomplex_radix2_inverse,
           RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_halfcomplex_radix2_inverse2(VALUE obj)
{
  return rb_fft_radix2(obj, gsl_fft_halfcomplex_radix2_inverse,
           RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_halfcomplex_radix2_backward(VALUE obj)
{
  return rb_fft_radix2(obj, gsl_fft_halfcomplex_radix2_backward,
           RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_halfcomplex_radix2_backward2(VALUE obj)
{
  return rb_fft_radix2(obj, gsl_fft_halfcomplex_radix2_backward,
           RB_GSL_FFT_INPLACE);
}

/*****/

static VALUE rb_fft_real_trans(int argc, VALUE *argv, VALUE obj,
             int (*trans)(double [], size_t, size_t,
              const gsl_fft_real_wavetable *,
              gsl_fft_real_workspace *),
             int sss)
{
  int flag = 0, naflag = 0;
  // local variable "status" was defined and set, but never used
  //int status;
  size_t stride, n;
  gsl_vector *vnew;
  gsl_vector_view vv;
  double *ptr1, *ptr2;
#ifdef HAVE_NARRAY_H
  int shape[1];
#endif
  gsl_fft_real_wavetable *table = NULL;
  gsl_fft_real_workspace *space = NULL;
  VALUE ary;
  flag = gsl_fft_get_argv_real(argc, argv, obj, &ptr1, &stride, &n, &table, &space, &naflag);
  if (naflag == 0) {
    if (sss == RB_GSL_FFT_COPY) {
      vnew = gsl_vector_alloc(n);
      vv.vector.data = ptr1;
      vv.vector.stride = stride;
      vv.vector.size = n;
      gsl_vector_memcpy(vnew, &vv.vector);
      ptr2 = vnew->data;
      stride = 1;
      ary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else {
      ptr2 = ptr1;
      ary = obj;
    }
#ifdef HAVE_NARRAY_H
  } else if (naflag == 1) {
    if (sss == RB_GSL_FFT_COPY) {
      shape[0] = n;
      ary = na_make_object(NA_DFLOAT, 1, shape, cNArray);
      ptr2 = NA_PTR_TYPE(ary, double*);
      memcpy(ptr2, ptr1, sizeof(double)*n);
      stride = 1;
    } else {
      ptr2 = ptr1;
      ary = obj;
    }
#endif
  } else {
    rb_raise(rb_eRuntimeError, "something wrong");
  }
  /*status =*/ (*trans)(ptr2, stride, n, table, space);
  gsl_fft_free(flag, (GSL_FFT_Wavetable *) table, (GSL_FFT_Workspace *) space);
  return ary;
}

static VALUE rb_gsl_fft_real_transform(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_real_trans(argc, argv, obj, gsl_fft_real_transform,
         RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_real_transform2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_real_trans(argc, argv, obj, gsl_fft_real_transform,
         RB_GSL_FFT_INPLACE);
}

static VALUE rb_fft_halfcomplex_trans(int argc, VALUE *argv, VALUE obj,
              int (*trans)(double [], size_t, size_t,
               const gsl_fft_halfcomplex_wavetable *, gsl_fft_real_workspace *),
              int sss)
{
  int flag = 0, naflag = 0;
  // local variable "status" was defined and set, but never used
  //int status;
  size_t stride, n;
  gsl_vector *vnew;
  gsl_vector_view vv;
  double *ptr1, *ptr2;
#ifdef HAVE_NARRAY_H
  int shape[1];
#endif
  gsl_fft_halfcomplex_wavetable *table = NULL;
  gsl_fft_real_workspace *space = NULL;
  VALUE ary;
  flag = gsl_fft_get_argv_halfcomplex(argc, argv, obj, &ptr1, &stride, &n,
         &table, &space, &naflag);
  if (naflag == 0) {
    if (sss == RB_GSL_FFT_COPY) {
      vnew = gsl_vector_alloc(n);
      vv.vector.data = ptr1;
      vv.vector.stride = stride;
      vv.vector.size = n;
      gsl_vector_memcpy(vnew, &vv.vector);
      ptr2 = vnew->data;
      stride = 1;
      ary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else {
      ptr2 = ptr1;
      ary = obj;
    }
#ifdef HAVE_NARRAY_H
  } else if (naflag == 1) {
    if (sss == RB_GSL_FFT_COPY) {
      shape[0] = n;
      ary = na_make_object(NA_DFLOAT, 1, shape, cNArray);
      ptr2 = NA_PTR_TYPE(ary, double*);
      memcpy(ptr2, ptr1, sizeof(double)*n);
      stride = 1;
    } else {
      ptr2 = ptr1;
      ary = obj;
    }
#endif
  } else {
    rb_raise(rb_eRuntimeError, "something wrong");
  }
  /*status =*/ (*trans)(ptr2, stride, n, table, space);
  gsl_fft_free(flag, (GSL_FFT_Wavetable *) table, (GSL_FFT_Workspace *) space);
  return ary;
}

static VALUE rb_gsl_fft_halfcomplex_transform(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_halfcomplex_trans(argc, argv, obj,
          gsl_fft_halfcomplex_transform,
          RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_halfcomplex_transform2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_halfcomplex_trans(argc, argv, obj,
          gsl_fft_halfcomplex_transform,
          RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_halfcomplex_backward(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_halfcomplex_trans(argc, argv, obj,
          gsl_fft_halfcomplex_backward,
          RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_halfcomplex_backward2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_halfcomplex_trans(argc, argv, obj,
          gsl_fft_halfcomplex_backward,
          RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_halfcomplex_inverse(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_halfcomplex_trans(argc, argv, obj,
          gsl_fft_halfcomplex_inverse,
          RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_halfcomplex_inverse2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_halfcomplex_trans(argc, argv, obj,
          gsl_fft_halfcomplex_inverse,
          RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_real_unpack(VALUE obj)
{
  gsl_vector *v;
  gsl_vector_complex *vout;

  CHECK_VECTOR(obj);
  Data_Get_Struct(obj, gsl_vector, v);

  vout = gsl_vector_complex_alloc(v->size);
  gsl_fft_real_unpack(v->data, (gsl_complex_packed_array) vout->data, v->stride, v->size);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vout);
}

static VALUE rb_gsl_fft_halfcomplex_unpack(VALUE obj)
{
  gsl_vector *v;
  gsl_vector_complex *vout;

  CHECK_VECTOR(obj);
  Data_Get_Struct(obj, gsl_vector, v);

  vout = gsl_vector_complex_alloc(v->size);
  gsl_fft_halfcomplex_unpack(v->data, (gsl_complex_packed_array) vout->data, v->stride, v->size);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vout);
}

/* Convert a halfcomplex data to Numerical Recipes style */
static VALUE rb_gsl_fft_halfcomplex_to_nrc(VALUE obj)
{
  gsl_vector *v, *vnew;
  size_t i, k;

  CHECK_VECTOR(obj);
  Data_Get_Struct(obj, gsl_vector, v);

  vnew = gsl_vector_alloc(v->size);
  gsl_vector_set(vnew, 0, gsl_vector_get(v, 0));  /* DC */
  gsl_vector_set(vnew, 1, gsl_vector_get(v, v->size/2));  /* Nyquist freq */
  for (i = 2, k = 1; i < vnew->size; i+=2, k++) {
    gsl_vector_set(vnew, i, gsl_vector_get(v, k));
    gsl_vector_set(vnew, i+1, -gsl_vector_get(v, v->size-k));
  }
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
}

static VALUE rb_gsl_fft_halfcomplex_amp_phase(VALUE obj)
{
  gsl_vector *v;
  gsl_vector *amp, *phase;
  double re, im;
  VALUE vamp, vphase;
  size_t i;
  CHECK_VECTOR(obj);
  Data_Get_Struct(obj, gsl_vector, v);
  amp = gsl_vector_alloc(v->size/2);
  phase = gsl_vector_alloc(v->size/2);
  gsl_vector_set(amp, 0, gsl_vector_get(v, 0));
  gsl_vector_set(phase, 0, 0);
  gsl_vector_set(amp, amp->size-1, gsl_vector_get(v, v->size-1));
  gsl_vector_set(phase, phase->size-1, 0);
  for (i = 1; i < v->size-1; i+=2) {
    re = gsl_vector_get(v, i);
    im = gsl_vector_get(v, i+1);
    gsl_vector_set(amp, i/2+1, sqrt(re*re + im*im));
    gsl_vector_set(phase, i/2+1, atan2(im, re));
  }
  vamp = Data_Wrap_Struct(VECTOR_ROW_COL(obj), 0, gsl_vector_free, amp);
  vphase = Data_Wrap_Struct(VECTOR_ROW_COL(obj), 0, gsl_vector_free, phase);
  return rb_ary_new3(2, vamp, vphase);
}

void Init_gsl_fft(VALUE module)
{
  mgsl_fft = rb_define_module_under(module, "FFT");

  /*****/

  rb_define_const(mgsl_fft, "Forward", INT2FIX(gsl_fft_forward));
  rb_define_const(mgsl_fft, "FORWARD", INT2FIX(gsl_fft_forward));
  rb_define_const(mgsl_fft, "Backward", INT2FIX(gsl_fft_backward));
  rb_define_const(mgsl_fft, "BACKWARD", INT2FIX(gsl_fft_backward));

  /* Transforms for complex vectors */
  rb_define_method(cgsl_vector_complex, "radix2_forward",
       rb_gsl_fft_complex_radix2_forward, 0);
  rb_define_method(cgsl_vector_complex, "radix2_transform",
       rb_gsl_fft_complex_radix2_transform, 1);
  rb_define_method(cgsl_vector_complex, "radix2_backward",
       rb_gsl_fft_complex_radix2_backward, 0);
  rb_define_method(cgsl_vector_complex, "radix2_inverse",
       rb_gsl_fft_complex_radix2_inverse, 0);
  rb_define_method(cgsl_vector_complex, "radix2_dif_forward",
       rb_gsl_fft_complex_radix2_dif_forward, 0);
  rb_define_method(cgsl_vector_complex, "radix2_dif_transform",
       rb_gsl_fft_complex_radix2_dif_transform, 1);
  rb_define_method(cgsl_vector_complex, "radix2_dif_backward",
       rb_gsl_fft_complex_radix2_dif_backward, 0);
  rb_define_method(cgsl_vector_complex, "radix2_dif_inverse",
       rb_gsl_fft_complex_radix2_dif_inverse, 0);

  /* In-place radix-2 transforms for complex vectors */
  rb_define_method(cgsl_vector_complex, "radix2_forward!",
       rb_gsl_fft_complex_radix2_forward2, 0);
  rb_define_method(cgsl_vector_complex, "radix2_transform!",
       rb_gsl_fft_complex_radix2_transform2, 1);
  rb_define_method(cgsl_vector_complex, "radix2_backward!",
       rb_gsl_fft_complex_radix2_backward2, 0);
  rb_define_method(cgsl_vector_complex, "radix2_inverse!",
       rb_gsl_fft_complex_radix2_inverse2, 0);
  rb_define_method(cgsl_vector_complex, "radix2_dif_forward!",
       rb_gsl_fft_complex_radix2_dif_forward2, 0);
  rb_define_method(cgsl_vector_complex, "radix2_dif_transform!",
       rb_gsl_fft_complex_radix2_dif_transform2, 1);
  rb_define_method(cgsl_vector_complex, "radix2_dif_backward!",
       rb_gsl_fft_complex_radix2_dif_backward2, 0);
  rb_define_method(cgsl_vector_complex, "radix2_dif_inverse!",
       rb_gsl_fft_complex_radix2_dif_inverse2, 0);

  // class GSL::FFT::Wavetable < GSL::Object
  //
  // Useful since some functionality is shared among
  // GSL::FFT::Complex::Wavetable
  // GSL::FFT::Real::Wavetable
  // GSL::FFT::HalfComplex::Wavetable
  cgsl_fft_wavetable = rb_define_class_under(mgsl_fft, "Wavetable", cGSL_Object);
  // No alloc
  // TODO Make GSL::FFT::Wavetable#initialize private?
  rb_define_method(cgsl_fft_wavetable, "n",
           rb_GSL_FFT_Wavetable_n, 0);
  rb_define_method(cgsl_fft_wavetable, "nf",
           rb_GSL_FFT_Wavetable_nf, 0);
  rb_define_method(cgsl_fft_wavetable, "factor",
       rb_GSL_FFT_Wavetable_factor, 0);

  // class GSL::FFT::ComplexWavetable < GSL::FFT::Wavetable
  cgsl_fft_complex_wavetable = rb_define_class_under(mgsl_fft, "ComplexWavetable",
                 cgsl_fft_wavetable);
  rb_define_singleton_method(cgsl_fft_complex_wavetable, "alloc",
           rb_gsl_fft_complex_wavetable_new, 1);

  // class GSL::FFT::ComplexWorkspace < GSL::Object
  cgsl_fft_complex_workspace = rb_define_class_under(mgsl_fft, "ComplexWorkspace",
                 cGSL_Object);
  rb_define_singleton_method(cgsl_fft_complex_workspace, "alloc",
           rb_gsl_fft_complex_workspace_new, 1);

  rb_define_method(cgsl_vector_complex, "forward", rb_gsl_fft_complex_forward, -1);
  rb_define_method(cgsl_vector_complex, "transform", rb_gsl_fft_complex_transform, -1);
  rb_define_method(cgsl_vector_complex, "backward", rb_gsl_fft_complex_backward, -1);
  rb_define_method(cgsl_vector_complex, "inverse", rb_gsl_fft_complex_inverse, -1);

  rb_define_method(cgsl_vector_complex, "forward!", rb_gsl_fft_complex_forward2, -1);
  rb_define_method(cgsl_vector_complex, "transform!", rb_gsl_fft_complex_transform2, -1);
  rb_define_method(cgsl_vector_complex, "backward!", rb_gsl_fft_complex_backward2, -1);
  rb_define_method(cgsl_vector_complex, "inverse!", rb_gsl_fft_complex_inverse2, -1);

  /*****/

  // TODO Do these method names need the "real_" and "halfcomplex_" prefixes?
  rb_define_method(cgsl_vector, "real_radix2_transform",
       rb_gsl_fft_real_radix2_transform, 0);
  rb_define_alias(cgsl_vector, "radix2_transform", "real_radix2_transform");
  rb_define_alias(cgsl_vector, "radix2_forward", "real_radix2_transform");
  rb_define_method(cgsl_vector, "real_radix2_inverse",
       rb_gsl_fft_halfcomplex_radix2_inverse, 0);
  rb_define_alias(cgsl_vector, "radix2_inverse", "real_radix2_inverse");
  rb_define_alias(cgsl_vector, "halfcomplex_radix2_inverse",
      "real_radix2_inverse");
  rb_define_method(cgsl_vector, "real_radix2_backward",
       rb_gsl_fft_halfcomplex_radix2_backward, 0);
  rb_define_alias(cgsl_vector, "radix2_backward", "real_radix2_backward");
  rb_define_alias(cgsl_vector, "halfcomplex_radix2_backward",
      "real_radix2_backward");

  // TODO Do these method names need the "real_" and "halfcomplex_" prefixes?
  rb_define_method(cgsl_vector, "real_radix2_transform!",
       rb_gsl_fft_real_radix2_transform2, 0);
  rb_define_alias(cgsl_vector, "radix2_transform!", "real_radix2_transform!");
  rb_define_alias(cgsl_vector, "radix2_forward!", "real_radix2_transform!");
  rb_define_method(cgsl_vector, "real_radix2_inverse!",
       rb_gsl_fft_halfcomplex_radix2_inverse2, 0);
  rb_define_alias(cgsl_vector, "radix2_inverse!", "real_radix2_inverse!");
  rb_define_alias(cgsl_vector, "halfcomplex_radix2_inverse!",
      "real_radix2_inverse!");
  rb_define_method(cgsl_vector, "real_radix2_backward!",
       rb_gsl_fft_halfcomplex_radix2_backward2, 0);
  rb_define_alias(cgsl_vector, "radix2_backward!", "real_radix2_backward!");
  rb_define_alias(cgsl_vector, "halfcomplex_radix2_backward!",
      "real_radix2_backward!");

  /*****/

  // class GSL::FFT::RealWavetable < GSL::FFT::Wavetable
  cgsl_fft_real_wavetable = rb_define_class_under(mgsl_fft, "RealWavetable",
              cgsl_fft_wavetable);
  rb_define_singleton_method(cgsl_fft_real_wavetable, "alloc",
           rb_gsl_fft_real_wavetable_new, 1);

  // class GSL::FFT::HalfComplexWavetable < GSL::FFT::Wavetable
  cgsl_fft_halfcomplex_wavetable = rb_define_class_under(mgsl_fft,
               "HalfComplexWavetable", cgsl_fft_wavetable);
  rb_define_singleton_method(cgsl_fft_halfcomplex_wavetable, "alloc",
           rb_gsl_fft_halfcomplex_wavetable_new, 1);

  /*****/

  // class GSL::FFT::RealWorkspace < GSL::Object
  cgsl_fft_real_workspace = rb_define_class_under(mgsl_fft, "RealWorkspace",
              cGSL_Object);
  rb_define_singleton_method(cgsl_fft_real_workspace, "alloc",
           rb_gsl_fft_real_workspace_new, 1);

  /*****/

  // TODO Do these method names need the "real_" and "halfcomplex_" prefixes?
  rb_define_method(cgsl_vector, "real_transform", rb_gsl_fft_real_transform, -1);
  rb_define_alias(cgsl_vector, "transform", "real_transform");
  rb_define_alias(cgsl_vector, "forward", "real_transform");
  rb_define_alias(cgsl_vector, "fft_forward", "real_transform");
  rb_define_alias(cgsl_vector, "fft", "real_transform");
  rb_define_method(cgsl_vector, "halfcomplex_transform",
       rb_gsl_fft_halfcomplex_transform, -1);
  rb_define_method(cgsl_vector, "halfcomplex_backward",
       rb_gsl_fft_halfcomplex_backward, -1);
  rb_define_alias(cgsl_vector, "backward", "halfcomplex_backward");
  rb_define_alias(cgsl_vector, "fft_backward", "halfcomplex_backward");
  rb_define_method(cgsl_vector, "halfcomplex_inverse",
       rb_gsl_fft_halfcomplex_inverse, -1);
  rb_define_alias(cgsl_vector, "fft_inverse", "halfcomplex_inverse");
  rb_define_alias(cgsl_vector, "ifft", "halfcomplex_inverse");
  rb_define_alias(cgsl_vector, "inverse", "halfcomplex_inverse");

  rb_define_method(cgsl_vector, "real_transform!", rb_gsl_fft_real_transform2, -1);
  rb_define_alias(cgsl_vector, "transform!", "real_transform!");
  rb_define_alias(cgsl_vector, "forward!", "real_transform!");
  rb_define_alias(cgsl_vector, "fft_forward!", "real_transform!");
  rb_define_alias(cgsl_vector, "fft!", "real_transform!");
  rb_define_method(cgsl_vector, "halfcomplex_transform!",
       rb_gsl_fft_halfcomplex_transform2, -1);
  rb_define_method(cgsl_vector, "halfcomplex_backward!",
       rb_gsl_fft_halfcomplex_backward2, -1);
  rb_define_alias(cgsl_vector, "backward!", "halfcomplex_backward!");
  rb_define_alias(cgsl_vector, "fft_backward!", "halfcomplex_backward!");
  rb_define_method(cgsl_vector, "halfcomplex_inverse!",
       rb_gsl_fft_halfcomplex_inverse2, -1);
  rb_define_alias(cgsl_vector, "fft_inverse!", "halfcomplex_inverse!");
  rb_define_alias(cgsl_vector, "ifft!", "halfcomplex_inverse!");
  rb_define_alias(cgsl_vector, "inverse!", "halfcomplex_inverse!");

  /***/
  rb_define_method(cgsl_vector, "fft_real_unpack", rb_gsl_fft_real_unpack, 0);
  rb_define_alias(cgsl_vector, "real_unpack", "fft_real_unpack");
  rb_define_alias(cgsl_vector, "real_to_complex", "fft_real_unpack");
  rb_define_alias(cgsl_vector, "r_to_c", "fft_real_unpack");

  rb_define_method(cgsl_vector, "fft_halfcomplex_unpack",
       rb_gsl_fft_halfcomplex_unpack, 0);
  rb_define_alias(cgsl_vector, "halfcomplex_unpack", "fft_halfcomplex_unpack");
  rb_define_alias(cgsl_vector, "halfcomplex_to_complex", "fft_halfcomplex_unpack");
  rb_define_alias(cgsl_vector, "hc_to_c", "fft_halfcomplex_unpack");

  /*****/

  rb_define_method(cgsl_vector, "to_nrc_order",
           rb_gsl_fft_halfcomplex_to_nrc, 0);

  rb_define_method(cgsl_vector, "halfcomplex_amp_phase",
           rb_gsl_fft_halfcomplex_amp_phase, 0);
  rb_define_alias(cgsl_vector, "hc_amp_phase", "halfcomplex_amp_phase");
}
