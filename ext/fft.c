/*
  fft.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_config.h"
#include "rb_gsl_fft.h"

VALUE mgsl_fft_complex;
VALUE mgsl_fft, mgsl_fft_real, mgsl_fft_halfcomplex;
VALUE cgsl_cparray;
VALUE cgsl_fft_wavetable, cgsl_fft_workspace;
VALUE cgsl_fft_wavetable_factor;
VALUE cgsl_fft_complex_wavetable, cgsl_fft_complex_workspace;
VALUE cgsl_fft_real_wavetable, cgsl_fft_halfcomplex_wavetable;
VALUE cgsl_fft_real_workspace;
VALUE cgsl_vector_halfcomplex;

static GSL_FFT_Workspace* GSL_FFT_Workspace_new(size_t n);
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

static GSL_FFT_Workspace* GSL_FFT_Workspace_new(size_t n)
{
  CHECK_FIXNUM(n);
  return (GSL_FFT_Workspace *) gsl_fft_complex_workspace_alloc(n);
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
  
VALUE rb_gsl_vector_new(int argc, VALUE *argv, VALUE klass);
VALUE rb_gsl_permutation_new(VALUE klass, VALUE nn);
static void get_stride_n(int argc, VALUE *argv, int argstart, gsl_vector *v,
			 size_t *stride, size_t *n);

static VALUE cparray_get(VALUE obj, VALUE ii)
{
  gsl_vector *v = NULL;
  gsl_complex *c = NULL;
  size_t i;
  CHECK_FIXNUM(ii);
  CHECK_VECTOR(obj);
  Data_Get_Struct(obj, gsl_vector, v);
  i = FIX2INT(ii);
  c = ALLOC(gsl_complex);
  GSL_SET_REAL(c, gsl_vector_get(v, 2*i));
  GSL_SET_IMAG(c, gsl_vector_get(v, 2*i+1));
  return Data_Wrap_Struct(cgsl_complex, 0, free, c); 
}

static VALUE cparray_set(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v = NULL;
  gsl_complex *c = NULL;
  size_t i;
  if (argc < 2) rb_raise(rb_eArgError, "too few arguments");
  CHECK_VECTOR(obj);
  CHECK_FIXNUM(argv[0]);
  i = FIX2INT(argv[0]);
  Data_Get_Struct(obj, gsl_vector, v);
  if (rb_obj_is_kind_of(argv[1], cgsl_complex)) {
    Data_Get_Struct(argv[1], gsl_complex, c);
    gsl_vector_set(v, 2*i, GSL_REAL(*c));
    gsl_vector_set(v, 2*i+1, GSL_IMAG(*c));
    return obj;
  }
  switch (TYPE(argv[1])) {
  case T_ARRAY:
    gsl_vector_set(v, 2*i, NUM2DBL(rb_ary_entry(argv[1], 0)));
    gsl_vector_set(v, 2*i+1, NUM2DBL(rb_ary_entry(argv[1], 1)));
    break;
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    gsl_vector_set(v, 2*i, NUM2DBL(argv[1]));
    break;
  default:
    rb_raise(rb_eTypeError, "wrong type arguments");
    break;
  }
  return obj;
}

static VALUE cparray_real(VALUE obj, VALUE ii)
{
  gsl_vector *v = NULL;
  CHECK_VECTOR(obj);
  CHECK_FIXNUM(ii);
  Data_Get_Struct(obj, gsl_vector, v);
  return rb_float_new(gsl_vector_get(v, 2*FIX2INT(ii)));
}

static VALUE cparray_re(VALUE obj)
{
  gsl_vector *v = NULL;
  gsl_vector_view *vv = NULL;
  Data_Get_Struct(obj, gsl_vector, v);
  vv = gsl_vector_view_alloc();
  *vv = gsl_vector_subvector_with_stride(v, 0, 2, v->size/2);
  return Data_Wrap_Struct(cgsl_vector_view, 0, gsl_vector_view_free, vv);
}

static VALUE cparray_im(VALUE obj)
{
  gsl_vector *v = NULL;
  gsl_vector_view *vv = NULL;
  Data_Get_Struct(obj, gsl_vector, v);
  vv = gsl_vector_view_alloc();
  *vv = gsl_vector_subvector_with_stride(v, 1, 2, v->size/2);
  return Data_Wrap_Struct(cgsl_vector_view, 0, gsl_vector_view_free, vv);
}

static VALUE cparray_set_real(VALUE obj, VALUE ii, VALUE val)
{
  gsl_vector *v = NULL;
  CHECK_VECTOR(obj);
  CHECK_FIXNUM(ii);
  Need_Float(val);
  Data_Get_Struct(obj, gsl_vector, v);
  gsl_vector_set(v, 2*FIX2INT(ii), NUM2DBL(val));
  return obj;
}

static VALUE cparray_imag(VALUE obj, VALUE ii)
{
  gsl_vector *v = NULL;
  CHECK_VECTOR(obj);
  CHECK_FIXNUM(ii);
  Data_Get_Struct(obj, gsl_vector, v);
  return rb_float_new(gsl_vector_get(v, 2*FIX2INT(ii)+1));
}

static VALUE cparray_set_imag(VALUE obj, VALUE ii, VALUE val)
{
  gsl_vector *v = NULL;
  CHECK_VECTOR(obj);
  CHECK_FIXNUM(ii);
  Need_Float(val);
  Data_Get_Struct(obj, gsl_vector, v);
  gsl_vector_set(v, 2*FIX2INT(ii)+1, NUM2DBL(val));
  return obj;
}

static VALUE get_cpary_stride_n(int argc, VALUE *argv, VALUE obj,
			       gsl_complex_packed_array *data, size_t *stride, size_t *n)
{
  gsl_vector *v = NULL;
  int itmp = 0;
  VALUE ary;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 1)",
			   argc);
    if (obj == mgsl_fft_complex) {
      if (CLASS_OF(argv[0]) != cgsl_cparray)
	rb_raise(rb_eTypeError, "wrong argument type %s (expected PackedArray)",
		 rb_class2name(CLASS_OF(argv[0])));
    }
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, v);
    ary = argv[0];
    itmp = 1;
    break;
  default:
    CHECK_VECTOR(obj);
    Data_Get_Struct(obj, gsl_vector, v);
    ary = obj;
    itmp = 0;
    break;
  }
  *data = (gsl_complex_packed_array) v->data;

  switch (argc-itmp) {
  case 0:
    *stride = v->stride;
    *n = v->size/2;
    break;
  case 1:
    CHECK_FIXNUM(argv[itmp]);
    *n = FIX2INT(argv[itmp]);
    *stride = v->stride;
    break;
  default:
    CHECK_FIXNUM(argv[itmp]);
    CHECK_FIXNUM(argv[itmp+1]);
    *stride = FIX2INT(argv[itmp]);
    *n = FIX2INT(argv[itmp+1]);
    break;
  }
  return ary;
}

static void get_stride_n(int argc, VALUE *argv, int argstart, gsl_vector *v,
			 size_t *stride, size_t *n)
{
  switch (argc-argstart) {
  case 0:
    *stride = v->stride;
    *n = v->size;
    break;
  case 1:
    CHECK_FIXNUM(argv[argstart]);
    *stride = v->stride;
    *n = FIX2INT(argv[argstart]);
    break;
  default:
    CHECK_FIXNUM(argv[argstart]);
    CHECK_FIXNUM(argv[argstart+1]);
    *stride = FIX2INT(argv[argstart]);
    *n = FIX2INT(argv[argstart+1]);
    break;
  }
  return;
}

static VALUE rb_fft_complex_radix2(int argc, VALUE *argv, VALUE obj,
				   int (*trans)(gsl_complex_packed_array, 
						size_t, size_t), int flag)
{
  size_t stride, n;
  gsl_complex_packed_array data;
  gsl_vector *v;
  VALUE ary;
  ary = get_cpary_stride_n(argc, argv, obj, &data, &stride, &n);
  if (flag == RB_GSL_FFT_COPY) {
    v = gsl_vector_alloc(2*n);
    memcpy(v->data, data, sizeof(double)*2*n);
    (*trans)(v->data, stride, n);
    return Data_Wrap_Struct(cgsl_cparray, 0, gsl_vector_free, v);
  } else { /* in-place */
    (*trans)(data, stride, n);
    return ary;
  }
}

static VALUE rb_gsl_fft_complex_radix2_forward(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_radix2(argc, argv, obj, gsl_fft_complex_radix2_forward,
			       RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_radix2_forward2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_radix2(argc, argv, obj, gsl_fft_complex_radix2_forward,
			       RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_complex_radix2_transform(int argc, VALUE *argv, VALUE obj)
{
  size_t stride, n;
  gsl_complex_packed_array data;
  gsl_fft_direction sign;
  gsl_vector *v;
  CHECK_FIXNUM(argv[argc-1]);
  sign = FIX2INT(argv[argc-1]);
  get_cpary_stride_n(argc-1, argv, obj, &data, &stride, &n);
  v = gsl_vector_alloc(2*n);
  memcpy(v->data, data, sizeof(double)*2*n);
  gsl_fft_complex_radix2_transform(v->data, stride, n, sign);
  return Data_Wrap_Struct(cgsl_cparray, 0, gsl_vector_free, v);
}

static VALUE rb_gsl_fft_complex_radix2_transform2(int argc, VALUE *argv, VALUE obj)
{
  size_t stride, n;
  gsl_complex_packed_array data;
  gsl_fft_direction sign;
  VALUE ary;
  CHECK_FIXNUM(argv[argc-1]);
  sign = FIX2INT(argv[argc-1]);
  ary = get_cpary_stride_n(argc-1, argv, obj, &data, &stride, &n);
  gsl_fft_complex_radix2_transform(data, stride, n, sign);
  return ary;
}

static VALUE rb_gsl_fft_complex_radix2_backward(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_radix2(argc, argv, obj, gsl_fft_complex_radix2_backward,
			       RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_radix2_inverse(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_radix2(argc, argv, obj, gsl_fft_complex_radix2_inverse,
			       RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_radix2_dif_forward(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_radix2(argc, argv, obj, 
			       gsl_fft_complex_radix2_dif_forward,
			       RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_radix2_backward2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_radix2(argc, argv, obj, gsl_fft_complex_radix2_backward,
			       RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_complex_radix2_inverse2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_radix2(argc, argv, obj, gsl_fft_complex_radix2_inverse,
			       RB_GSL_FFT_INPLACE);
}


static VALUE rb_gsl_fft_complex_radix2_dif_forward2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_radix2(argc, argv, obj, 
			       gsl_fft_complex_radix2_dif_forward,
			       RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_complex_radix2_dif_transform(int argc, VALUE *argv, VALUE obj)
{
  size_t stride, n;
  gsl_complex_packed_array data;
  gsl_vector *v;
  gsl_fft_direction sign;
  CHECK_FIXNUM(argv[argc-1]);
  sign = FIX2INT(argv[argc-1]);
  get_cpary_stride_n(argc-1, argv, obj, &data, &stride, &n);
  v = gsl_vector_alloc(2*n);
  memcpy(v->data, data, sizeof(double)*2*n);
  gsl_fft_complex_radix2_dif_transform(v->data, stride, n, sign);
  return Data_Wrap_Struct(cgsl_cparray, 0, gsl_vector_free, v);
}

/* in-place */
static VALUE rb_gsl_fft_complex_radix2_dif_transform2(int argc, VALUE *argv, VALUE obj)
{
  size_t stride, n;
  gsl_complex_packed_array data;
  gsl_fft_direction sign;
  VALUE ary;
  CHECK_FIXNUM(argv[argc-1]);
  sign = FIX2INT(argv[argc-1]);
  ary = get_cpary_stride_n(argc-1, argv, obj, &data, &stride, &n);
  gsl_fft_complex_radix2_dif_transform(data, stride, n, sign);
  return ary;
}

static VALUE rb_gsl_fft_complex_radix2_dif_backward(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_radix2(argc, argv, obj, 
			       gsl_fft_complex_radix2_dif_backward,
			       RB_GSL_FFT_COPY);
}
static VALUE rb_gsl_fft_complex_radix2_dif_inverse(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_radix2(argc, argv, obj, 
			       gsl_fft_complex_radix2_dif_inverse,
			       RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_complex_radix2_dif_backward2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_radix2(argc, argv, obj, 
			       gsl_fft_complex_radix2_dif_backward,
			       RB_GSL_FFT_INPLACE);
}
static VALUE rb_gsl_fft_complex_radix2_dif_inverse2(int argc, VALUE *argv, VALUE obj)
{
  return rb_fft_complex_radix2(argc, argv, obj, 
			       gsl_fft_complex_radix2_dif_inverse,
			       RB_GSL_FFT_INPLACE);
}

static VALUE rb_GSL_FFT_Wavetable_new(VALUE klass, VALUE n)
{
  GSL_FFT_Wavetable *table;
  CHECK_FIXNUM(n);
  table = (GSL_FFT_Wavetable *) gsl_fft_complex_wavetable_alloc(FIX2INT(n));
  return Data_Wrap_Struct(klass, 0, GSL_FFT_Wavetable_free, table);
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
  gsl_permutation *p;
  size_t i;
  Data_Get_Struct(obj, GSL_FFT_Wavetable, table);
  p = gsl_permutation_alloc(64);
  for (i = 0; i < table->nf; i++) p->data[i] = table->factor[i];
  return Data_Wrap_Struct(cgsl_fft_wavetable_factor, 0, gsl_permutation_free, p);
}

static VALUE rb_GSL_FFT_Workspace_new(VALUE klass, VALUE n)
{
  GSL_FFT_Workspace *space;
  CHECK_FIXNUM(n);
  space = GSL_FFT_Workspace_new(FIX2INT(n));
  return Data_Wrap_Struct(klass, 0, GSL_FFT_Workspace_free, space);
}

enum {
  NONE_OF_TWO = 0,
  ALLOC_SPACE = 1,
  ALLOC_TABLE = 2,
  BOTH_OF_TWO = 3,
} FFTComplexStructAllocFlag;

static void gsl_fft_free(int flag, GSL_FFT_Wavetable *table,
			 GSL_FFT_Workspace *space);
static int gsl_fft_get_argv(int argc, VALUE *argv, VALUE obj,
			    gsl_complex_packed_array *data, size_t *stride,
			    size_t *n, gsl_fft_complex_wavetable **table,
			    gsl_fft_complex_workspace **space);

static int gsl_fft_get_argv(int argc, VALUE *argv, VALUE obj,
			    gsl_complex_packed_array *data, size_t *stride,
			    size_t *n, gsl_fft_complex_wavetable **table,
			    gsl_fft_complex_workspace **space)
{
  int flag = NONE_OF_TWO, flagtmp, i, itmp = argc, itmp2 = 0, ccc;
  int flagw = 0;
  gsl_vector *v = NULL;

  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 1)",
			   argc);
    if (obj == mgsl_fft_complex) {
      if (CLASS_OF(argv[0]) != cgsl_cparray)
	rb_raise(rb_eTypeError, "wrong argument type %s (expected PackedArray)",
		 rb_class2name(CLASS_OF(argv[0])));
    }
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, v);
    itmp2 = 1;
    break;
  default:
    CHECK_VECTOR(obj);
    Data_Get_Struct(obj, gsl_vector, v);
    break;
  }
  ccc = argc;
  flagtmp = 0;
  flagw = 0;
  for (i = argc-1; i >= itmp2; i--) {
    if (rb_obj_is_kind_of(argv[i], cgsl_fft_workspace)
	|| rb_obj_is_kind_of(argv[i], cgsl_fft_complex_workspace)) {
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
    if (rb_obj_is_kind_of(argv[i], cgsl_fft_wavetable)
	|| rb_obj_is_kind_of(argv[i], cgsl_fft_complex_wavetable)) {
      Data_Get_Struct(argv[i], gsl_fft_complex_wavetable, *table);
      flagtmp = 1;
      ccc--;
      break;
    }
  }
  get_cpary_stride_n(ccc, argv, obj, data, stride, n);
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

static int gsl_fft_get_argv2(int argc, VALUE *argv, VALUE obj,
			     double **ptr, size_t *stride,
			     size_t *n, gsl_fft_real_wavetable **table,
			     gsl_fft_real_workspace **space, int *naflag)
{
  int flag = NONE_OF_TWO, flagtmp, i, itmp = argc, itmp2 = 0, ccc;
  int flagw = 0;
  *naflag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 1)",
			   argc);
    if (obj == mgsl_fft_complex) {
      if (CLASS_OF(argv[0]) != cgsl_cparray)
	rb_raise(rb_eTypeError, "wrong argument type %s (expected PackedArray)",
		 rb_class2name(CLASS_OF(argv[0])));
    }
    *ptr = get_ptr_double3(argv[0], n, stride, naflag);
    itmp2 = 1;
    break;
  default:
    *ptr = get_ptr_double3(obj, n, stride, naflag);
    break;
  }
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

static int gsl_fft_get_argv3(int argc, VALUE *argv, VALUE obj,
			   double **ptr, size_t *stride,
			    size_t *n, gsl_fft_halfcomplex_wavetable **table,
			    gsl_fft_real_workspace **space, int *naflag)
{
  int flag = NONE_OF_TWO, flagtmp, i, itmp = argc, itmp2 = 0, ccc;
  int flagw = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 1)",
			   argc);
    if (obj == mgsl_fft_complex) {
      if (CLASS_OF(argv[0]) != cgsl_cparray)
	rb_raise(rb_eTypeError, "wrong argument type %s (expected PackedArray)",
		 rb_class2name(CLASS_OF(argv[0])));
    }
    *ptr = get_ptr_double3(argv[0], n, stride, naflag);
    itmp2 = 1;
    break;
  default:
    *ptr = get_ptr_double3(obj, n, stride, naflag);
    break;
  }
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

static VALUE rb_gsl_fft_getary(int argc, VALUE *argv, VALUE obj)
{
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    return argv[0];
    break;
  default:
    return obj;
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
  int flag = 0, status;
  size_t stride, n;
  gsl_complex_packed_array data;
  gsl_vector *v;
  gsl_fft_complex_wavetable *table = NULL;
  gsl_fft_complex_workspace *space = NULL;
  flag = gsl_fft_get_argv(argc, argv, obj, &data, &stride, &n, &table, &space);
  if (sss == RB_GSL_FFT_COPY) {
    v = gsl_vector_alloc(2*n);
    memcpy(v->data, data, sizeof(double)*2*n);
    status = (*transform)(v->data, stride, n, table, space);
    gsl_fft_free(flag, (GSL_FFT_Wavetable *) table, (GSL_FFT_Workspace *) space);
    return Data_Wrap_Struct(cgsl_cparray, 0, gsl_vector_free, v);
  } else {    /* in-place */
    status = (*transform)(data, stride, n, table, space);
    gsl_fft_free(flag, (GSL_FFT_Wavetable *) table, (GSL_FFT_Workspace *) space);
    return rb_gsl_fft_getary(argc, argv, obj);
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
  int flag = 0, status;
  size_t stride, n;
  gsl_vector *v;
  gsl_fft_direction sign;
  gsl_complex_packed_array data;
  gsl_fft_complex_wavetable *table = NULL;
  gsl_fft_complex_workspace *space = NULL;
  CHECK_FIXNUM(argv[argc-1]);
  sign = FIX2INT(argv[argc-1]);
  flag = gsl_fft_get_argv(argc-1, argv, obj, &data, &stride, &n, &table, &space);
  v = gsl_vector_alloc(2*n);
  memcpy(v->data, data, sizeof(double)*2*n);
  status = gsl_fft_complex_transform(v->data, stride, n, table, space, sign);
  gsl_fft_free(flag, (GSL_FFT_Wavetable *) table, (GSL_FFT_Workspace *) space);
  return Data_Wrap_Struct(cgsl_cparray, 0, gsl_vector_free, v);
}

/* in-place */
static VALUE rb_gsl_fft_complex_transform2(int argc, VALUE *argv, VALUE obj)
{
  int flag = 0, status;
  size_t stride, n;
  gsl_fft_direction sign;
  gsl_complex_packed_array data;
  gsl_fft_complex_wavetable *table = NULL;
  gsl_fft_complex_workspace *space = NULL;
  CHECK_FIXNUM(argv[argc-1]);
  sign = FIX2INT(argv[argc-1]);
  flag = gsl_fft_get_argv(argc-1, argv, obj, &data, &stride, &n, &table, &space);
  status = gsl_fft_complex_transform(data, stride, n, table, space, sign);
  gsl_fft_free(flag, (GSL_FFT_Wavetable *) table, (GSL_FFT_Workspace *) space);
  return rb_gsl_fft_getary(argc, argv, obj);
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

static VALUE get_ptr_stride_n(int argc, VALUE *argv, VALUE obj,
			     double **ptr, size_t *stride, size_t *n, int *flag)
{
  int itmp = 0;
  VALUE ary;
  *flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 1)",
			   argc);
    if (obj == mgsl_fft_complex) {
      if (CLASS_OF(argv[0]) != cgsl_cparray)
	rb_raise(rb_eTypeError, "wrong argument type %s (expected PackedArray)",
		 rb_class2name(CLASS_OF(argv[0])));
    }
    *ptr = get_ptr_double3(argv[0], n, stride, flag);
    itmp = 1;
    ary = argv[0];
    break;
  default:
    *ptr = get_ptr_double3(obj, n, stride, flag);
    ary = obj;
    itmp = 0;
    break;
  }
  switch (argc-itmp) {
  case 0:
    /* do nothing */
    break;
  case 1:
    CHECK_FIXNUM(argv[itmp]);
    *n = FIX2INT(argv[itmp]);
    break;
  default:
    CHECK_FIXNUM(argv[itmp]);
    CHECK_FIXNUM(argv[itmp+1]);
    *stride = FIX2INT(argv[itmp]);
    *n = FIX2INT(argv[itmp+1]);
    break;
  }
  return ary;
}

static VALUE rb_fft_radix2(int argc, VALUE *argv, VALUE obj,
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
  get_ptr_stride_n(argc, argv, obj, &ptr1, &stride, &n, &flag);
  if (flag == 0) {
    if (sss == RB_GSL_FFT_COPY) {
      vnew = gsl_vector_alloc(n);
      vv.vector.data = ptr1;
      vv.vector.stride = stride;
      vv.vector.size = n;
      gsl_vector_memcpy(vnew, &vv.vector);
      ptr2 = vnew->data;
      ary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else {
      ary = rb_gsl_fft_getary(argc, argv, obj);
      ptr2 = ptr1;
    }
#ifdef HAVE_NARRAY_H
  } else if (flag == 1) {
    if (sss == RB_GSL_FFT_COPY) {
      shape[0] = n;
      ary = na_make_object(NA_DFLOAT, 1, shape, cNArray);
      ptr2 = NA_PTR_TYPE(ary, double*);
      memcpy(ptr2, ptr1, sizeof(double)*n);
    } else {
      ary = rb_gsl_fft_getary(argc, argv, obj);
      ptr2 = NA_PTR_TYPE(ary, double*);
    }
#endif
  } else {
    rb_raise(rb_eRuntimeError, "something wrong");
  }
  (*trans)(ptr2, stride, n);
  return ary;
}

static void rb_fft_radix2_check_arg(int argc, VALUE obj)
{
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 1 && argc != 3) rb_raise(rb_eArgError, 
					 "wrong number of arguments (%d for 1 or 3)", argc);
    break;
  default:
    /*    CHECK_VECTOR(obj);*/
    if (argc != 0 && argc != 2) rb_raise(rb_eArgError, 
					 "wrong number of arguments (%d for 0 or 2)", argc);
    break;
  }
  return;
}

static VALUE rb_gsl_fft_real_radix2_transform(int argc, VALUE *argv, VALUE obj)
{
  rb_fft_radix2_check_arg(argc, obj);
  return rb_fft_radix2(argc, argv, obj, gsl_fft_real_radix2_transform,
		       RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_real_radix2_transform2(int argc, VALUE *argv, VALUE obj)
{
  rb_fft_radix2_check_arg(argc, obj);
  return rb_fft_radix2(argc, argv, obj, gsl_fft_real_radix2_transform,
		       RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_halfcomplex_radix2_inverse(int argc, VALUE *argv, VALUE obj)
{
  rb_fft_radix2_check_arg(argc, obj);
  return rb_fft_radix2(argc, argv, obj, gsl_fft_halfcomplex_radix2_inverse,
		       RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_halfcomplex_radix2_inverse2(int argc, VALUE *argv, VALUE obj)
{
  rb_fft_radix2_check_arg(argc, obj);
  return rb_fft_radix2(argc, argv, obj, gsl_fft_halfcomplex_radix2_inverse,
		       RB_GSL_FFT_INPLACE);
}

static VALUE rb_gsl_fft_halfcomplex_radix2_backward(int argc, VALUE *argv, VALUE obj)
{
  rb_fft_radix2_check_arg(argc, obj);
  return rb_fft_radix2(argc, argv, obj, gsl_fft_halfcomplex_radix2_backward,
		       RB_GSL_FFT_COPY);
}

static VALUE rb_gsl_fft_halfcomplex_radix2_backward2(int argc, VALUE *argv, VALUE obj)
{
  rb_fft_radix2_check_arg(argc, obj);
  return rb_fft_radix2(argc, argv, obj, gsl_fft_halfcomplex_radix2_backward,
		       RB_GSL_FFT_INPLACE);
}

/*****/

static VALUE rb_fft_real_trans(int argc, VALUE *argv, VALUE obj,
			       int (*trans)(double [], size_t, size_t, 
					    const gsl_fft_real_wavetable *, 
					    gsl_fft_real_workspace *),
			       int sss)
{
  int flag = 0, status, naflag = 0;
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
  flag = gsl_fft_get_argv2(argc, argv, obj, &ptr1, &stride, &n, &table, &space, &naflag);
  if (naflag == 0) {
    if (sss == RB_GSL_FFT_COPY) {
      vnew = gsl_vector_alloc(n);
      vv.vector.data = ptr1;
      vv.vector.stride = stride;
      vv.vector.size = n;
      gsl_vector_memcpy(vnew, &vv.vector);
      ptr2 = vnew->data;
      ary = Data_Wrap_Struct(cgsl_vector_halfcomplex, 0, gsl_vector_free, vnew);
    } else {
      ptr2 = ptr1;
      ary = rb_gsl_fft_getary(argc, argv, obj);
    }
#ifdef HAVE_NARRAY_H
  } else if (naflag == 1) {
    if (sss == RB_GSL_FFT_COPY) {
      shape[0] = n;
      ary = na_make_object(NA_DFLOAT, 1, shape, cNArray);
      ptr2 = NA_PTR_TYPE(ary, double*);
      memcpy(ptr2, ptr1, sizeof(double)*n);
    } else {
      ptr2 = ptr1;
      ary = rb_gsl_fft_getary(argc, argv, obj);
    }
#endif
  } else {
    rb_raise(rb_eRuntimeError, "something wrong");
  }
  status = (*trans)(ptr2, stride, n, table, space);
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
  int flag = 0, status, naflag = 0;
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
  flag = gsl_fft_get_argv3(argc, argv, obj, &ptr1, &stride, &n, 
			   &table, &space, &naflag);
  if (naflag == 0) {
    if (sss == RB_GSL_FFT_COPY) {
      vnew = gsl_vector_alloc(n);
      vv.vector.data = ptr1;
      vv.vector.stride = stride;
      vv.vector.size = n;
      gsl_vector_memcpy(vnew, &vv.vector);
      ptr2 = vnew->data;
      ary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else {
      ptr2 = ptr1;
      ary = rb_gsl_fft_getary(argc, argv, obj);
    }
#ifdef HAVE_NARRAY_H
  } else if (naflag == 1) {
    if (sss == RB_GSL_FFT_COPY) {
      shape[0] = n;
      ary = na_make_object(NA_DFLOAT, 1, shape, cNArray);
      ptr2 = NA_PTR_TYPE(ary, double*);
      memcpy(ptr2, ptr1, sizeof(double)*n);
    } else {
      ptr2 = ptr1;
      ary = rb_gsl_fft_getary(argc, argv, obj);
    }
#endif
  } else {
    rb_raise(rb_eRuntimeError, "something wrong");
  }
  status = (*trans)(ptr2, stride, n, table, space);
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

static VALUE rb_gsl_fft_real_unpack(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v, *cpary;
  size_t stride, n;

  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 1)",
			   argc);
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, v);
    get_stride_n(argc-1, argv, 1, v, &stride, &n);
    break;
  default:
    CHECK_VECTOR(obj);
    Data_Get_Struct(obj, gsl_vector, v);
    get_stride_n(argc, argv, 0, v, &stride, &n);
    break;
  }
  cpary = gsl_vector_alloc(2*n);
  gsl_fft_real_unpack(v->data, (gsl_complex_packed_array) cpary->data, stride, n);
  return Data_Wrap_Struct(cgsl_cparray, 0, gsl_vector_free, cpary);
}

static VALUE rb_gsl_fft_halfcomplex_unpack(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v, *cpary;
  size_t stride, n;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 1)",
			   argc);
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, v);
    get_stride_n(argc-1, argv, 1, v, &stride, &n);
    break;
  default:
    Data_Get_Struct(obj, gsl_vector, v);
    get_stride_n(argc, argv, 0, v, &stride, &n);
    break;
  }
  cpary = gsl_vector_alloc(2*n);
  gsl_fft_halfcomplex_unpack(v->data, (gsl_complex_packed_array) cpary->data, stride, n);
  return Data_Wrap_Struct(cgsl_cparray, 0, gsl_vector_free, cpary);
}

/* Convert a halfcomplex data to Numerical Recipes style */
static VALUE rb_gsl_fft_halfcomplex_to_nrc(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v, *vnew;
  size_t i, k;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
			    argc);
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, v);
    break;
  default:
    if (argc != 0) rb_raise(rb_eArgError, "wrong number of arguments (%d for 0)",
			    argc);
    CHECK_VECTOR(obj);
    Data_Get_Struct(obj, gsl_vector, v);
    break;
  }
  vnew = gsl_vector_alloc(v->size);
  gsl_vector_set(vnew, 0, gsl_vector_get(v, 0));  /* DC */
  gsl_vector_set(vnew, 1, gsl_vector_get(v, v->size/2));  /* Nyquist freq */
  for (i = 2, k = 1; i < vnew->size; i+=2, k++) {
    gsl_vector_set(vnew, i, gsl_vector_get(v, k));
    gsl_vector_set(vnew, i+1, -gsl_vector_get(v, v->size-k));
  }
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
}

static VALUE rb_gsl_fft_halfcomplex_amp_phase(VALUE module, VALUE obj)
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
  VALUE mgsl_fft_complex_radix2, mgsl_fft_real_radix2;
  VALUE mgsl_fft_halfcomplex_radix2;

  mgsl_fft = rb_define_module_under(module, "FFT");
  mgsl_fft_complex = rb_define_module_under(mgsl_fft, "Complex");
  mgsl_fft_complex_radix2 = rb_define_module_under(mgsl_fft_complex, "Radix2");
  mgsl_fft_real = rb_define_module_under(mgsl_fft, "Real");
  mgsl_fft_real_radix2 = rb_define_module_under(mgsl_fft_real, "Radix2");
  mgsl_fft_halfcomplex = rb_define_module_under(mgsl_fft, "HalfComplex");
  mgsl_fft_halfcomplex_radix2 = rb_define_module_under(mgsl_fft_halfcomplex, "Radix2");

  cgsl_vector_halfcomplex = rb_define_class_under(cgsl_vector, "HalfComplex",
						  cgsl_vector);

  /*****/

  rb_define_const(mgsl_fft, "Forward", INT2FIX(forward));
  rb_define_const(mgsl_fft, "FORWARD", INT2FIX(forward));
  rb_define_const(mgsl_fft, "Backward", INT2FIX(backward));
  rb_define_const(mgsl_fft, "BACKWARD", INT2FIX(backward));

  cgsl_cparray = rb_define_class_under(mgsl_fft_complex, "PackedArray", cgsl_vector);

  rb_define_singleton_method(cgsl_cparray, "alloc", rb_gsl_vector_new, -1);

  rb_define_method(cgsl_cparray, "get", cparray_get, 1);
  rb_define_alias(cgsl_cparray, "[]", "get");
  rb_define_method(cgsl_cparray, "real", cparray_real, 1);
  rb_define_method(cgsl_cparray, "re", cparray_re, 0);
  rb_define_method(cgsl_cparray, "imag", cparray_imag, 1);
  rb_define_method(cgsl_cparray, "im", cparray_im, 0);
  rb_define_method(cgsl_cparray, "set", cparray_set, -1);
  rb_define_alias(cgsl_cparray, "[]=", "set");
  rb_define_method(cgsl_cparray, "set_real", cparray_set_real, 2);
  rb_define_method(cgsl_cparray, "set_imag", cparray_set_imag, 2);

  rb_define_method(cgsl_cparray, "radix2_forward", 
		   rb_gsl_fft_complex_radix2_forward, -1);
  rb_define_method(cgsl_cparray, "radix2_transform", 
		   rb_gsl_fft_complex_radix2_transform, -1);
  rb_define_method(cgsl_cparray, "radix2_backward", 
		   rb_gsl_fft_complex_radix2_backward, -1);
  rb_define_method(cgsl_cparray, "radix2_inverse", 
		   rb_gsl_fft_complex_radix2_inverse, -1);
  rb_define_method(cgsl_cparray, "radix2_dif_forward", 
		   rb_gsl_fft_complex_radix2_dif_forward, -1);
  rb_define_method(cgsl_cparray, "radix2_dif_transform", 
		   rb_gsl_fft_complex_radix2_dif_transform, -1);
  rb_define_method(cgsl_cparray, "radix2_dif_backward", 
		   rb_gsl_fft_complex_radix2_dif_backward, -1);
  rb_define_method(cgsl_cparray, "radix2_dif_inverse", 
		   rb_gsl_fft_complex_radix2_dif_inverse, -1);

  rb_define_method(cgsl_cparray, "radix2_forward!", 
		   rb_gsl_fft_complex_radix2_forward2, -1);
  rb_define_method(cgsl_cparray, "radix2_transform!",
		   rb_gsl_fft_complex_radix2_transform2, -1);
  rb_define_method(cgsl_cparray, "radix2_backward!", 
		   rb_gsl_fft_complex_radix2_backward2, -1);
  rb_define_method(cgsl_cparray, "radix2_inverse!", 
		   rb_gsl_fft_complex_radix2_inverse2, -1);
  rb_define_method(cgsl_cparray, "radix2_dif_forward!", 
		   rb_gsl_fft_complex_radix2_dif_forward2, -1);
  rb_define_method(cgsl_cparray, "radix2_dif_transform!",
		   rb_gsl_fft_complex_radix2_dif_transform2, -1);
  rb_define_method(cgsl_cparray, "radix2_dif_backward!", 
		   rb_gsl_fft_complex_radix2_dif_backward2, -1);
  rb_define_method(cgsl_cparray, "radix2_dif_inverse!", 
		   rb_gsl_fft_complex_radix2_dif_inverse2, -1);

  cgsl_fft_wavetable = rb_define_class_under(mgsl_fft, "Wavetable", cGSL_Object);
  rb_define_singleton_method(cgsl_fft_wavetable, "alloc",
			     rb_GSL_FFT_Wavetable_new, 1);
  rb_define_method(cgsl_fft_wavetable, "n",
			     rb_GSL_FFT_Wavetable_n, 0);
  rb_define_method(cgsl_fft_wavetable, "nf",
			     rb_GSL_FFT_Wavetable_nf, 0);
  rb_define_method(cgsl_fft_wavetable, "factor",
		   rb_GSL_FFT_Wavetable_factor, 0);
  cgsl_fft_wavetable_factor = rb_define_class_under(mgsl_fft, 
						    "WavetableFactor", 
						    cgsl_permutation);

  cgsl_fft_complex_wavetable = rb_define_class_under(mgsl_fft_complex, "Wavetable",
						     cgsl_fft_wavetable);
  rb_define_singleton_method(cgsl_fft_complex_wavetable, "alloc",
			     rb_gsl_fft_complex_wavetable_new, 1);

  cgsl_fft_workspace = rb_define_class_under(mgsl_fft, "Workspace", cGSL_Object);
  rb_define_singleton_method(cgsl_fft_workspace, "alloc",
			     rb_GSL_FFT_Workspace_new, 1);

  cgsl_fft_complex_workspace = rb_define_class_under(mgsl_fft_complex, "Workspace", 
						     cgsl_fft_workspace);
  rb_define_singleton_method(cgsl_fft_complex_workspace, "alloc",
			     rb_gsl_fft_complex_workspace_new, 1);

  rb_define_method(cgsl_cparray, "forward", rb_gsl_fft_complex_forward, -1);
  rb_define_method(cgsl_cparray, "transform", rb_gsl_fft_complex_transform, -1);
  rb_define_method(cgsl_cparray, "backward", rb_gsl_fft_complex_backward, -1);
  rb_define_method(cgsl_cparray, "inverse", rb_gsl_fft_complex_inverse, -1);

  rb_define_method(cgsl_cparray, "forward!", rb_gsl_fft_complex_forward2, -1);
  rb_define_method(cgsl_cparray, "transform!", rb_gsl_fft_complex_transform2, -1);
  rb_define_method(cgsl_cparray, "backward!", rb_gsl_fft_complex_backward2, -1);
  rb_define_method(cgsl_cparray, "inverse!", rb_gsl_fft_complex_inverse2, -1);

  /*****/

  rb_define_method(cgsl_vector, "real_radix2_transform", 
		   rb_gsl_fft_real_radix2_transform, -1);
  rb_define_alias(cgsl_vector, "radix2_transform", "real_radix2_transform");
  rb_define_alias(cgsl_vector, "radix2_forward", "real_radix2_transform");
  rb_define_method(cgsl_vector, "real_radix2_inverse", 
		   rb_gsl_fft_halfcomplex_radix2_inverse, -1);
  rb_define_alias(cgsl_vector, "radix2_inverse", "real_radix2_inverse");
  rb_define_alias(cgsl_vector, "halfcomplex_radix2_inverse", 
		  "real_radix2_inverse");
  rb_define_method(cgsl_vector, "real_radix2_backward", 
		   rb_gsl_fft_halfcomplex_radix2_backward, -1);
  rb_define_alias(cgsl_vector, "radix2_backward", "real_radix2_backward");
  rb_define_alias(cgsl_vector, "halfcomplex_radix2_backward", 
		  "real_radix2_backward");


  rb_define_method(cgsl_vector, "real_radix2_transform!", 
		   rb_gsl_fft_real_radix2_transform2, -1);
  rb_define_alias(cgsl_vector, "radix2_transform!", "real_radix2_transform!");
  rb_define_alias(cgsl_vector, "radix2_forward!", "real_radix2_transform!");
  rb_define_method(cgsl_vector, "real_radix2_inverse!", 
		   rb_gsl_fft_halfcomplex_radix2_inverse2, -1);
  rb_define_alias(cgsl_vector, "radix2_inverse!", "real_radix2_inverse!");
  rb_define_alias(cgsl_vector, "halfcomplex_radix2_inverse!", 
		  "real_radix2_inverse!");
  rb_define_method(cgsl_vector, "real_radix2_backward!", 
		   rb_gsl_fft_halfcomplex_radix2_backward2, -1);
  rb_define_alias(cgsl_vector, "radix2_backward!", "real_radix2_backward!");
  rb_define_alias(cgsl_vector, "halfcomplex_radix2_backward!", 
		  "real_radix2_backward!");

  /*****/
  cgsl_fft_real_wavetable = rb_define_class_under(mgsl_fft_real, "Wavetable", 
						  cgsl_fft_wavetable);
  rb_define_singleton_method(cgsl_fft_real_wavetable, "alloc",
			     rb_gsl_fft_real_wavetable_new, 1);

  cgsl_fft_halfcomplex_wavetable = rb_define_class_under(mgsl_fft_halfcomplex, 
							 "Wavetable", cgsl_fft_wavetable);
  rb_define_singleton_method(cgsl_fft_halfcomplex_wavetable, "alloc",
			     rb_gsl_fft_halfcomplex_wavetable_new, 1);

  /*****/
  cgsl_fft_real_workspace = rb_define_class_under(mgsl_fft_real, "Workspace", 
						  cgsl_fft_workspace);
  rb_define_singleton_method(cgsl_fft_real_workspace, "alloc",
			     rb_gsl_fft_real_workspace_new, 1);

  /*****/
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
  rb_define_method(cgsl_vector, "fft_real_unpack", rb_gsl_fft_real_unpack, -1);
  rb_define_alias(cgsl_vector, "real_unpack", "fft_real_unpack");
  rb_define_method(cgsl_vector, "fft_halfcomplex_unpack", 
		   rb_gsl_fft_halfcomplex_unpack, -1);
  rb_define_alias(cgsl_vector, "halfcomplex_unpack", "fft_halfcomplex_unpack");

  /*****/

  rb_define_singleton_method(mgsl_fft_complex, "radix2_forward", 
			     rb_gsl_fft_complex_radix2_forward, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_transform", 
			     rb_gsl_fft_complex_radix2_transform, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_backward", 
			     rb_gsl_fft_complex_radix2_backward, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_inverse", 
			     rb_gsl_fft_complex_radix2_inverse, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_dif_forward", 
			     rb_gsl_fft_complex_radix2_dif_forward, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_dif_transform", 
			     rb_gsl_fft_complex_radix2_dif_transform, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_dif_backward", 
			     rb_gsl_fft_complex_radix2_dif_backward, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_dif_inverse", 
			     rb_gsl_fft_complex_radix2_dif_inverse, -1);

  rb_define_singleton_method(mgsl_fft_complex_radix2, "forward", 
			     rb_gsl_fft_complex_radix2_forward, -1);
 rb_define_singleton_method(mgsl_fft_complex_radix2, "transform", 
			     rb_gsl_fft_complex_radix2_transform, -1);
  rb_define_singleton_method(mgsl_fft_complex_radix2, "backward", 
			     rb_gsl_fft_complex_radix2_backward, -1);
  rb_define_singleton_method(mgsl_fft_complex_radix2, "inverse", 
			     rb_gsl_fft_complex_radix2_inverse, -1);
  rb_define_singleton_method(mgsl_fft_complex_radix2, "dif_forward", 
			     rb_gsl_fft_complex_radix2_dif_forward, -1);
  rb_define_singleton_method(mgsl_fft_complex_radix2, "dif_transform", 
			     rb_gsl_fft_complex_radix2_dif_transform, -1);
  rb_define_singleton_method(mgsl_fft_complex_radix2, "dif_backward", 
			     rb_gsl_fft_complex_radix2_dif_backward, -1);
  rb_define_singleton_method(mgsl_fft_complex_radix2, "dif_inverse", 
			     rb_gsl_fft_complex_radix2_dif_inverse, -1);

  rb_define_singleton_method(mgsl_fft_complex, "forward", 
			     rb_gsl_fft_complex_forward, -1);
  rb_define_singleton_method(mgsl_fft_complex, "transform", 
			     rb_gsl_fft_complex_transform, -1);
  rb_define_singleton_method(mgsl_fft_complex, "backward", 
			     rb_gsl_fft_complex_backward, -1);
  rb_define_singleton_method(mgsl_fft_complex, "inverse", 
			     rb_gsl_fft_complex_inverse, -1);

  rb_define_singleton_method(mgsl_fft_real, "radix2_transform",
			     rb_gsl_fft_real_radix2_transform, -1);
  rb_define_singleton_method(mgsl_fft_real, "radix2_forward",
			     rb_gsl_fft_real_radix2_transform, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex, "radix2_inverse",
			     rb_gsl_fft_halfcomplex_radix2_inverse, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex, "radix2_backward",
			     rb_gsl_fft_halfcomplex_radix2_backward, -1);

  rb_define_singleton_method(mgsl_fft_real_radix2, "transform",
			     rb_gsl_fft_real_radix2_transform, -1);
  rb_define_singleton_method(mgsl_fft_real_radix2, "forward",
			     rb_gsl_fft_real_radix2_transform, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex_radix2, "inverse",
			     rb_gsl_fft_halfcomplex_radix2_inverse, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex_radix2, "backward",
			     rb_gsl_fft_halfcomplex_radix2_backward, -1);

  rb_define_singleton_method(mgsl_fft_real, "transform",
			     rb_gsl_fft_real_transform, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex, "transform",
			     rb_gsl_fft_halfcomplex_transform, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex, "inverse",
			     rb_gsl_fft_halfcomplex_inverse, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex, "backward",
			     rb_gsl_fft_halfcomplex_backward, -1);

  /***/

  rb_define_singleton_method(mgsl_fft_complex, "radix2_forward!", 
			     rb_gsl_fft_complex_radix2_forward2, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_transform!", 
			     rb_gsl_fft_complex_radix2_transform2, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_backward!", 
			     rb_gsl_fft_complex_radix2_backward2, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_inverse!", 
			     rb_gsl_fft_complex_radix2_inverse2, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_dif_forward!", 
			     rb_gsl_fft_complex_radix2_dif_forward2, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_dif_transform!", 
			     rb_gsl_fft_complex_radix2_dif_transform2, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_dif_backward!", 
			     rb_gsl_fft_complex_radix2_dif_backward2, -1);
  rb_define_singleton_method(mgsl_fft_complex, "radix2_dif_inverse!", 
			     rb_gsl_fft_complex_radix2_dif_inverse2, -1);

  rb_define_singleton_method(mgsl_fft_complex_radix2, "forward!", 
			     rb_gsl_fft_complex_radix2_forward2, -1);
 rb_define_singleton_method(mgsl_fft_complex_radix2, "transform!", 
			     rb_gsl_fft_complex_radix2_transform2, -1);
  rb_define_singleton_method(mgsl_fft_complex_radix2, "backward!", 
			     rb_gsl_fft_complex_radix2_backward2, -1);
  rb_define_singleton_method(mgsl_fft_complex_radix2, "inverse!", 
			     rb_gsl_fft_complex_radix2_inverse2, -1);
  rb_define_singleton_method(mgsl_fft_complex_radix2, "dif_forward!", 
			     rb_gsl_fft_complex_radix2_dif_forward2, -1);
  rb_define_singleton_method(mgsl_fft_complex_radix2, "dif_transform!", 
			     rb_gsl_fft_complex_radix2_dif_transform2, -1);
  rb_define_singleton_method(mgsl_fft_complex_radix2, "dif_backward!", 
			     rb_gsl_fft_complex_radix2_dif_backward2, -1);
  rb_define_singleton_method(mgsl_fft_complex_radix2, "dif_inverse!", 
			     rb_gsl_fft_complex_radix2_dif_inverse2, -1);

  rb_define_singleton_method(mgsl_fft_complex, "forward!", 
			     rb_gsl_fft_complex_forward2, -1);
  rb_define_singleton_method(mgsl_fft_complex, "transform!", 
			     rb_gsl_fft_complex_transform2, -1);
  rb_define_singleton_method(mgsl_fft_complex, "backward!", 
			     rb_gsl_fft_complex_backward2, -1);
  rb_define_singleton_method(mgsl_fft_complex, "inverse!", 
			     rb_gsl_fft_complex_inverse2, -1);

  rb_define_singleton_method(mgsl_fft_real, "radix2_transform!",
			     rb_gsl_fft_real_radix2_transform2, -1);
  rb_define_singleton_method(mgsl_fft_real, "radix2_forward!",
			     rb_gsl_fft_real_radix2_transform2, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex, "radix2_inverse!",
			     rb_gsl_fft_halfcomplex_radix2_inverse2, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex, "radix2_backward!",
			     rb_gsl_fft_halfcomplex_radix2_backward2, -1);

  rb_define_singleton_method(mgsl_fft_real_radix2, "transform!",
			     rb_gsl_fft_real_radix2_transform2, -1);
  rb_define_singleton_method(mgsl_fft_real_radix2, "forward!",
			     rb_gsl_fft_real_radix2_transform2, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex_radix2, "inverse!",
			     rb_gsl_fft_halfcomplex_radix2_inverse2, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex_radix2, "backward!",
			     rb_gsl_fft_halfcomplex_radix2_backward2, -1);

  rb_define_singleton_method(mgsl_fft_real, "transform!",
			     rb_gsl_fft_real_transform2, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex, "transform!",
			     rb_gsl_fft_halfcomplex_transform2, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex, "inverse!",
			     rb_gsl_fft_halfcomplex_inverse2, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex, "backward!",
			     rb_gsl_fft_halfcomplex_backward2, -1);
  /***/
  rb_define_singleton_method(mgsl_fft_real, "unpack", 
			     rb_gsl_fft_real_unpack, -1);
  rb_define_singleton_method(mgsl_fft_halfcomplex, "unpack", 
			     rb_gsl_fft_halfcomplex_unpack, -1);

  /*****/
  rb_define_singleton_method(mgsl_fft_halfcomplex, "to_nrc_order",
			     rb_gsl_fft_halfcomplex_to_nrc, -1);
  rb_define_method(cgsl_vector, "to_nrc_order",
			     rb_gsl_fft_halfcomplex_to_nrc, -1);
			     
  rb_define_singleton_method(mgsl_fft_halfcomplex, "amp_phase",
			     rb_gsl_fft_halfcomplex_amp_phase, 1);

}
