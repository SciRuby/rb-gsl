/*
  histogram.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_config.h"
#include "rb_gsl_histogram.h"
#include "rb_gsl_array.h"
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

VALUE cgsl_histogram;
VALUE cgsl_histogram_range;
VALUE cgsl_histogram_bin;
static VALUE cgsl_histogram_integ;

static VALUE rb_gsl_histogram_alloc_from_file(VALUE klass, VALUE name);
#ifdef GSL_0_9_4_LATER
static VALUE rb_gsl_histogram_alloc(int argc, VALUE *argv, VALUE klass)
{
  gsl_histogram *h = NULL;
  gsl_vector *v;
  double min, max;
  size_t n, size;
  switch (argc) {
  case 1:
    switch (TYPE(argv[0])) {
    case T_FIXNUM:
      n = FIX2INT(argv[0]);
      h = gsl_histogram_alloc(n);
      break;
    case T_ARRAY:
      v = make_cvector_from_rarray(argv[0]);
      h = gsl_histogram_alloc(v->size-1);
      gsl_histogram_set_ranges(h, v->data, v->size);
      gsl_vector_free(v);
      break;
    case T_STRING:
      return rb_gsl_histogram_alloc_from_file(klass, argv[0]);
      break;
    default:
      CHECK_VECTOR(argv[0]);
      Data_Get_Struct(argv[0], gsl_vector, v);
      h = gsl_histogram_alloc(v->size-1);
      gsl_histogram_set_ranges(h, v->data, v->size);
      break;
    }
    return Data_Wrap_Struct(klass, 0, gsl_histogram_free, h);
    break;
  case 3:
    CHECK_FIXNUM(argv[0]);
    Need_Float(argv[1]); Need_Float(argv[2]);
    n = FIX2INT(argv[0]);
    min = NUM2DBL(argv[1]);
    max = NUM2DBL(argv[2]);
    h = gsl_histogram_calloc(n);
    gsl_histogram_set_ranges_uniform(h, min, max);
    return Data_Wrap_Struct(klass, 0, gsl_histogram_free, h);
    break;
  case 2:
    switch (TYPE(argv[0])) {
    case T_FIXNUM:
      CHECK_FIXNUM(argv[0]);
      if (TYPE(argv[1]) != T_ARRAY) {
	rb_raise(rb_eTypeError, "wrong argument type %s (Array expected)", 
		 rb_class2name(CLASS_OF(argv[1])));
      }
      n = FIX2INT(argv[0]);
      min = NUM2DBL(rb_ary_entry(argv[1], 0));
      max = NUM2DBL(rb_ary_entry(argv[1], 1));
      h = gsl_histogram_calloc(n);
      gsl_histogram_set_ranges_uniform(h, min, max);
      return Data_Wrap_Struct(klass, 0, gsl_histogram_free, h);
      break;
    case T_ARRAY:
      CHECK_FIXNUM(argv[1]);
      v = make_cvector_from_rarray(argv[0]);
      size = FIX2INT(argv[1]);
      h = gsl_histogram_calloc(size-1);
      gsl_histogram_set_ranges(h, v->data, size);
      gsl_vector_free(v);
      return Data_Wrap_Struct(klass, 0, gsl_histogram_free, h);
      break;
    default:
      CHECK_VECTOR(argv[0]);
      CHECK_FIXNUM(argv[1]);
      Data_Get_Struct(argv[0], gsl_vector, v);
      size = FIX2INT(argv[1]);
      h = gsl_histogram_calloc(size-1);
      gsl_histogram_set_ranges(h, v->data, size);
      return Data_Wrap_Struct(klass, 0, gsl_histogram_free, h);
      break;
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1, 2, 3)", argc);
    break;
  }

}

static VALUE rb_gsl_histogram_alloc_from_file(VALUE klass, VALUE name)
{
  char filename[1024], buf[1024];
  gsl_histogram *h;
  int nn;
  size_t n, i;
  FILE *fp = NULL;
  double upper;
  strcpy(filename, STR2CHARPTR(name));
  sprintf(buf, "wc %s", filename);
  fp = popen(buf, "r");
  if (fp == NULL) rb_raise(rb_eIOError, "popen failed.");
  if (fgets(buf, 1024, fp) == NULL)
    rb_sys_fail(0);
  pclose(fp);
  sscanf(buf, "%d", &nn);
  n = (size_t) nn;      /* vector length */
  fp = fopen(filename, "r");
  if (fp == NULL) rb_raise(rb_eIOError, "cannot open file %s.", filename);
  h = gsl_histogram_alloc(n);
  i = 0;
  while (fgets(buf, 1024, fp)) {
    sscanf(buf, "%lg %lg %lg", h->range+i, &upper, h->bin+i);
    i++;
  }
  h->range[n] = upper;
  fclose(fp);
  return Data_Wrap_Struct(klass, 0, gsl_histogram_free, h);
  
}

/* initialization + set uniform ranges (equal spacing from min to max) */
static VALUE rb_gsl_histogram_alloc_uniform(int argc, VALUE *argv, VALUE klass)
{
  gsl_histogram *h = NULL;
  double min, max, tmp;
  size_t n;
  switch (argc) {
  case 3:
    CHECK_FIXNUM(argv[0]);
    Need_Float(argv[1]); Need_Float(argv[2]);
    n = FIX2INT(argv[0]);
    min = NUM2DBL(argv[1]);
    max = NUM2DBL(argv[2]);
    break;
  case 2:
    CHECK_FIXNUM(argv[0]);
    n = FIX2INT(argv[0]);
    Check_Type(argv[1], T_ARRAY);
    min = NUM2DBL(rb_ary_entry(argv[1], 0));
    max = NUM2DBL(rb_ary_entry(argv[1], 1));
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  if (min > max) {
    tmp = min;
    min = max;
    max = tmp;
  }
  h = gsl_histogram_alloc(n);
  gsl_histogram_set_ranges_uniform(h, min, max);
  return Data_Wrap_Struct(klass, 0, gsl_histogram_free, h);
}

/* initialization + set ranges with a given spacing from min to max */
static VALUE rb_gsl_histogram_alloc_with_min_max_step(VALUE klass, VALUE vmin, 
					    VALUE vmax, VALUE ss)
{
  gsl_histogram *h = NULL;
  gsl_vector *v = NULL;
  double min, max, tmp, step;
  size_t i, n;
  Need_Float(vmin); Need_Float(vmax); Need_Float(ss);
  min = NUM2DBL(vmin);
  max = NUM2DBL(vmax);
  step = NUM2DBL(ss);
  if (min > max) {
    tmp = min;
    min = max;
    max = tmp;
  }
  n = (int) ((max - min)/step);
  h = gsl_histogram_alloc(n);
  v = gsl_vector_alloc(n + 1);
  for (i = 0; i < n + 1; i++) gsl_vector_set(v, i, min + step*i);
  gsl_histogram_set_ranges(h, v->data, v->size);
  gsl_vector_free(v);
  return Data_Wrap_Struct(klass, 0, gsl_histogram_free, h);
}
#endif

static VALUE rb_gsl_histogram_calloc(VALUE klass, VALUE nn)
{
  gsl_histogram *h = NULL;
  CHECK_FIXNUM(nn);
  h = gsl_histogram_calloc(FIX2INT(nn));
  return Data_Wrap_Struct(klass, 0, gsl_histogram_free, h);
}

static VALUE rb_gsl_histogram_calloc_range(int argc, VALUE *argv,  VALUE klass)
{
  gsl_histogram *h = NULL;
  gsl_vector *v = NULL;
  size_t n;
  switch (argc) {
  case 1:
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, v);
    n = v->size;
    break;
  case 2:
    CHECK_FIXNUM(argv[0]);
    CHECK_VECTOR(argv[1]);
    n = FIX2INT(argv[0]);
    Data_Get_Struct(argv[1], gsl_vector, v);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  h = gsl_histogram_calloc_range(n, v->data);
  return Data_Wrap_Struct(klass, 0, gsl_histogram_free, h);
}

static VALUE rb_gsl_histogram_bins(VALUE obj)
{
  gsl_histogram *h = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  return INT2FIX(gsl_histogram_bins(h));
}

static VALUE rb_gsl_histogram_set_ranges(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h = NULL;
  gsl_vector *v = NULL;
  size_t size;
  Data_Get_Struct(obj, gsl_histogram, h);
  if (argc != 1 && argc != 2) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  if (TYPE(argv[0]) == T_ARRAY) {
    v = make_cvector_from_rarray(argv[0]);
    if (argc == 1) size = v->size;
    else size = FIX2INT(argv[1]);
    gsl_histogram_set_ranges(h, v->data, size);
    gsl_vector_free(v);
  } else {
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, v);
    if (argc == 1) size = v->size;
    else size = FIX2INT(argv[1]);
    gsl_histogram_set_ranges(h, v->data, size);
  }
  return obj;
}

static VALUE rb_gsl_histogram_range(VALUE obj)
{
  gsl_histogram *h = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  v = gsl_vector_view_alloc();
  v->vector.data = h->range;
  v->vector.size = h->n + 1;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_histogram_range, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_histogram_bin(VALUE obj)
{
  gsl_histogram *h = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  v = gsl_vector_view_alloc();
  v->vector.data = h->bin;
  v->vector.size = h->n;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_histogram_bin, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_histogram_set_ranges_uniform(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h = NULL;
  double xmin, xmax;
  switch (argc) {
  case 1:
    Check_Type(argv[0], T_ARRAY);
    xmin = NUM2DBL(rb_ary_entry(argv[0], 0));
    xmax = NUM2DBL(rb_ary_entry(argv[0], 1));
    break;
  case 2:
    xmin = NUM2DBL(argv[0]);
    xmax = NUM2DBL(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  Data_Get_Struct(obj, gsl_histogram, h);
  gsl_histogram_set_ranges_uniform(h, xmin, xmax);
  return obj;
}

/* singleton */
static VALUE rb_gsl_histogram_memcpy(VALUE obj, VALUE vhdest, VALUE vhsrc)
{
  gsl_histogram *hdest = NULL, *hsrc = NULL;
  CHECK_HISTOGRAM(vhdest);
  CHECK_HISTOGRAM(vhsrc);
  Data_Get_Struct(vhdest, gsl_histogram, hdest);
  Data_Get_Struct(vhsrc, gsl_histogram, hsrc);
  gsl_histogram_memcpy(hdest, hsrc);
  return vhdest;
}

static VALUE rb_gsl_histogram_clone(VALUE obj)
{
  gsl_histogram *hsrc = NULL, *hnew = NULL;
  Data_Get_Struct(obj, gsl_histogram, hsrc);
  hnew = gsl_histogram_clone(hsrc);
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_histogram_free, hnew);
}

static VALUE rb_gsl_histogram_accumulate(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h = NULL;
  gsl_vector *v;
  gsl_vector_int *vi;
  size_t i;
  double weight = 1;
#ifdef HAVE_NARRAY_H
  double *ptr;
  size_t size, stride;
#endif
  switch (argc) {
  case 2:
    Need_Float(argv[1]);
    weight = NUM2DBL(argv[1]);
    break;
  case 1:
    weight = 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  Data_Get_Struct(obj, gsl_histogram, h);
  if (TYPE(argv[0]) == T_ARRAY) {
    //    for (i = 0; i < RARRAY(argv[0])->len; i++)
    for (i = 0; (int) i < RARRAY_LEN(argv[0]); i++)
      gsl_histogram_accumulate(h, NUM2DBL(rb_ary_entry(argv[0], i)), weight);
  } else if (VECTOR_P(argv[0])) {
    Data_Get_Struct(argv[0], gsl_vector, v);
    for (i = 0; i < v->size; i++)
      gsl_histogram_accumulate(h, gsl_vector_get(v, i), weight);
  } else if (VECTOR_INT_P(argv[0])) {
    Data_Get_Struct(argv[0], gsl_vector_int, vi);
    for (i = 0; i < vi->size; i++)
      gsl_histogram_accumulate(h, (double)gsl_vector_int_get(vi, i), weight);
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(argv[0])) {
    ptr = get_vector_ptr(argv[0], &stride, &size);
    for (i = 0; i < size; i++) 
      gsl_histogram_accumulate(h, ptr[i], weight);
#endif
  } else {
    gsl_histogram_accumulate(h, NUM2DBL(argv[0]), weight);
  }
  return argv[0];
}

static VALUE rb_gsl_histogram_accumulate2(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h = NULL;
  double weight = 1;
  double x;
  switch (argc) {
  case 2:
    Need_Float(argv[1]);
    weight = NUM2DBL(argv[1]);
    /* no break; */
  case 1:
    Need_Float(argv[0]);
    x = NUM2DBL(argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  Data_Get_Struct(obj, gsl_histogram, h);
  if (x < h->range[0]) x = h->range[0] + 4*GSL_DBL_EPSILON;
  if (x > h->range[h->n]) x = h->range[h->n] - 4*GSL_DBL_EPSILON;
  gsl_histogram_accumulate(h, x, weight);
  return argv[0];
}

static VALUE rb_gsl_histogram_get(VALUE obj, VALUE i)
{
  gsl_histogram *h = NULL;
  CHECK_FIXNUM(i);
  Data_Get_Struct(obj, gsl_histogram, h);
  return rb_float_new(gsl_histogram_get(h, FIX2INT(i)));
}

static VALUE rb_gsl_histogram_get_range(VALUE obj, VALUE i)
{
  gsl_histogram *h = NULL;
  double lower, upper;
  CHECK_FIXNUM(i);
  Data_Get_Struct(obj, gsl_histogram, h);
  gsl_histogram_get_range(h, FIX2INT(i), &lower, &upper);
  return rb_ary_new3(2, rb_float_new(lower), rb_float_new(upper));
}

static VALUE rb_gsl_histogram_max(VALUE obj)
{
  gsl_histogram *h = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  return rb_float_new(gsl_histogram_max(h));
}

static VALUE rb_gsl_histogram_min(VALUE obj)
{
  gsl_histogram *h = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  return rb_float_new(gsl_histogram_min(h));
}

static VALUE rb_gsl_histogram_reset(VALUE obj)
{
  gsl_histogram *h = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  gsl_histogram_reset(h);
  return obj;
}

static VALUE rb_gsl_histogram_find(VALUE obj, VALUE x)
{
  gsl_histogram *h = NULL;
  size_t i;
  Need_Float(x);
  Data_Get_Struct(obj, gsl_histogram, h);
  gsl_histogram_find(h, NUM2DBL(x), &i);
  return INT2FIX(i);
}

static VALUE rb_gsl_histogram_max_val(VALUE obj)
{
  gsl_histogram *h = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  return rb_float_new(gsl_histogram_max_val(h));
}

static VALUE rb_gsl_histogram_max_bin(VALUE obj)
{
  gsl_histogram *h = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  return INT2FIX(gsl_histogram_max_bin(h));
}

static VALUE rb_gsl_histogram_min_val(VALUE obj)
{
  gsl_histogram *h = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  return rb_float_new(gsl_histogram_min_val(h));
}

static VALUE rb_gsl_histogram_min_bin(VALUE obj)
{
  gsl_histogram *h = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  return INT2FIX(gsl_histogram_min_bin(h));
}

static VALUE rb_gsl_histogram_mean(VALUE obj)
{
  gsl_histogram *h = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  return rb_float_new(gsl_histogram_mean(h));
}

static VALUE rb_gsl_histogram_sigma(VALUE obj)
{
  gsl_histogram *h = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  return rb_float_new(gsl_histogram_sigma(h));
}

#ifndef GSL_1_1_LATER
double gsl_histogram_sum(const gsl_histogram * h)
{  
  double sum=0;  
  size_t i=0, n;  
  n=h->n;  
  while(i < n) sum += h->bin[i++];  
  return sum;
}
#endif

static VALUE rb_gsl_histogram_sum(VALUE obj)
{
  gsl_histogram *h = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  if (CLASS_OF(obj) == cgsl_histogram_integ)
    return rb_float_new(gsl_histogram_get(h, h->n-1));
  else
    return rb_float_new(gsl_histogram_sum(h));
}

static VALUE rb_gsl_histogram_normalize_bang(VALUE obj)
{
  gsl_histogram *h = NULL;
  double scale;
  Data_Get_Struct(obj, gsl_histogram, h);
  if (CLASS_OF(obj) == cgsl_histogram_integ)
    scale = 1.0/gsl_histogram_get(h, h->n-1);
  else
    scale = 1.0/gsl_histogram_sum(h);
  gsl_histogram_scale(h, scale);
  return obj;
}

static VALUE rb_gsl_histogram_normalize(VALUE obj)
{
  gsl_histogram *h = NULL, *hnew = NULL;
  Data_Get_Struct(obj, gsl_histogram, h);
  hnew = gsl_histogram_clone(h);
  return rb_gsl_histogram_normalize_bang(Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_histogram_free, hnew));
}

static VALUE rb_gsl_histogram_integral(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h = NULL;
  size_t istart = 0, iend, i = 0;
  double sum = 0.0;
  Data_Get_Struct(obj, gsl_histogram, h);
  switch (argc) {
  case 0:
    return rb_gsl_histogram_sum(obj);
    break;
  case 1:
    CHECK_FIXNUM(argv[0]);
    istart = 0; iend = FIX2INT(argv[0]);
    break;
  case 2:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    istart = FIX2INT(argv[0]);
    iend = FIX2INT(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0-2)", argc);
    break;
  }
  if (iend >= h->n) iend = h->n - 1;
  i = istart;
  while (i <= iend) sum += h->bin[i++];
  return rb_float_new(sum);
}

static VALUE rb_gsl_histogram_equal_bins_p(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h1 = NULL, *h2 = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
			    argc);
    CHECK_HISTOGRAM(argv[0]);
    CHECK_HISTOGRAM(argv[1]);
    Data_Get_Struct(argv[0], gsl_histogram, h1);
    Data_Get_Struct(argv[1], gsl_histogram, h2);
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
			    argc);
    Data_Get_Struct(obj, gsl_histogram, h1);
    CHECK_HISTOGRAM(argv[0]);
    Data_Get_Struct(argv[0], gsl_histogram, h2);
    break;
  }
  return INT2FIX(gsl_histogram_equal_bins_p(h1, h2));
}

static VALUE rb_gsl_histogram_equal_bins_p2(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h1 = NULL, *h2 = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
			    argc);
    CHECK_HISTOGRAM(argv[0]);
    CHECK_HISTOGRAM(argv[1]);
    Data_Get_Struct(argv[0], gsl_histogram, h1);
    Data_Get_Struct(argv[1], gsl_histogram, h2);
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
			    argc);
    Data_Get_Struct(obj, gsl_histogram, h1);
    CHECK_HISTOGRAM(argv[0]);
    Data_Get_Struct(argv[0], gsl_histogram, h2);
    break;
  }
  if (gsl_histogram_equal_bins_p(h1, h2)) return Qtrue;
  else return Qfalse;
}

static VALUE rb_gsl_histogram_add(VALUE obj, VALUE hh2)
{
  gsl_histogram *h1 = NULL, *h2 = NULL, *hnew = NULL;
  Data_Get_Struct(obj, gsl_histogram, h1);
  hnew = gsl_histogram_clone(h1);
  if (HISTOGRAM_P(hh2)) {
    Data_Get_Struct(hh2, gsl_histogram, h2);
    mygsl_histogram_add(hnew, h2);
  } else {
    Need_Float(hh2);
    gsl_histogram_shift(hnew, NUM2DBL(hh2));
  }
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_histogram_free, hnew);
}

static VALUE rb_gsl_histogram_add2(VALUE obj, VALUE hh2)
{
  gsl_histogram *h1 = NULL, *h2 = NULL;
  Data_Get_Struct(obj, gsl_histogram, h1);
  if (HISTOGRAM_P(hh2)) {
    Data_Get_Struct(hh2, gsl_histogram, h2);
    mygsl_histogram_add(h1, h2);
  } else {
    Need_Float(hh2);
    gsl_histogram_shift(h1, NUM2DBL(hh2));
  }
  return obj;
}

static VALUE rb_gsl_histogram_sub(VALUE obj, VALUE hh2)
{
  gsl_histogram *h1 = NULL, *h2 = NULL, *hnew = NULL;
  Data_Get_Struct(obj, gsl_histogram, h1);
  hnew = gsl_histogram_clone(h1);
  if (HISTOGRAM_P(hh2)) {
    Data_Get_Struct(hh2, gsl_histogram, h2);
    mygsl_histogram_sub(hnew, h2);
  } else {
    Need_Float(hh2);
    gsl_histogram_shift(hnew, -NUM2DBL(hh2));
  }
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_histogram_free, hnew);
}

static VALUE rb_gsl_histogram_sub2(VALUE obj, VALUE hh2)
{
  gsl_histogram *h1 = NULL, *h2 = NULL;
  Data_Get_Struct(obj, gsl_histogram, h1);
  if (HISTOGRAM_P(hh2)) {
    Data_Get_Struct(hh2, gsl_histogram, h2);
    mygsl_histogram_sub(h1, h2);
  } else {
    Need_Float(hh2);
    gsl_histogram_shift(h1, -NUM2DBL(hh2));
  }
  return obj;
}

static VALUE rb_gsl_histogram_mul(VALUE obj, VALUE hh2)
{
  gsl_histogram *h1 = NULL, *h2 = NULL, *hnew = NULL;
  Data_Get_Struct(obj, gsl_histogram, h1);
  hnew = gsl_histogram_clone(h1);
  if (HISTOGRAM_P(hh2)) {
    Data_Get_Struct(hh2, gsl_histogram, h2);
    mygsl_histogram_mul(hnew, h2);
  } else {
    Need_Float(hh2);
    gsl_histogram_scale(hnew, NUM2DBL(hh2));
  }
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_histogram_free, hnew);
}

static VALUE rb_gsl_histogram_mul2(VALUE obj, VALUE hh2)
{
  gsl_histogram *h1 = NULL, *h2 = NULL;
  Data_Get_Struct(obj, gsl_histogram, h1);
  if (HISTOGRAM_P(hh2)) {
    Data_Get_Struct(hh2, gsl_histogram, h2);
    mygsl_histogram_mul(h1, h2);
  } else {
    Need_Float(hh2);
    gsl_histogram_scale(h1, NUM2DBL(hh2));
  }
  return obj;
}

static VALUE rb_gsl_histogram_div(VALUE obj, VALUE hh2)
{
  gsl_histogram *h1 = NULL, *h2 = NULL, *hnew = NULL;
  Data_Get_Struct(obj, gsl_histogram, h1);
  hnew = gsl_histogram_clone(h1);
  if (HISTOGRAM_P(hh2)) {
    Data_Get_Struct(hh2, gsl_histogram, h2);
    mygsl_histogram_div(hnew, h2);
  } else {
    Need_Float(hh2);
    gsl_histogram_scale(hnew, 1.0/NUM2DBL(hh2));
  }
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_histogram_free, hnew);
}

static VALUE rb_gsl_histogram_div2(VALUE obj, VALUE hh2)
{
  gsl_histogram *h1 = NULL, *h2 = NULL;
  Data_Get_Struct(obj, gsl_histogram, h1);
  if (HISTOGRAM_P(hh2)) {
    Data_Get_Struct(hh2, gsl_histogram, h2);
    mygsl_histogram_div(h1, h2);
  } else {
    Need_Float(hh2);
    gsl_histogram_scale(h1, 1.0/NUM2DBL(hh2));
  }
  return obj;
}

static VALUE rb_gsl_histogram_scale_bang(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h = NULL;
  double scale;
  Data_Get_Struct(obj, gsl_histogram, h);
  switch (argc) {
  case 0:
    if (CLASS_OF(obj) == cgsl_histogram_integ)
      scale = 1.0/h->bin[h->n-1];
    else
      scale = 1.0/gsl_histogram_sum(h);
    break;
  case 1:
    scale = NUM2DBL(argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  gsl_histogram_scale(h, scale);
  return obj;
}

static VALUE rb_gsl_histogram_scale(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h = NULL, *hnew = NULL;
  double scale;
  Data_Get_Struct(obj, gsl_histogram, h);
  switch (argc) {
  case 0:
    if (CLASS_OF(obj) == cgsl_histogram_integ)
      scale = 1.0/h->bin[h->n-1];
    else
      scale = 1.0/gsl_histogram_sum(h);
    break;
  case 1:
    scale = NUM2DBL(argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  hnew = gsl_histogram_clone(h);
  gsl_histogram_scale(hnew, scale);
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_histogram_free, hnew);
}

static VALUE rb_gsl_histogram_shift(VALUE obj, VALUE shift)
{
  gsl_histogram *h = NULL;
  Need_Float(shift);
  Data_Get_Struct(obj, gsl_histogram, h);
  gsl_histogram_shift(h, NUM2DBL(shift));
  return obj;
}

static VALUE rb_gsl_histogram_shift2(VALUE obj, VALUE shift)
{
  gsl_histogram *h = NULL, *hnew = NULL;
  Need_Float(shift);
  Data_Get_Struct(obj, gsl_histogram, h);
  hnew = gsl_histogram_clone(h);
  gsl_histogram_shift(hnew, NUM2DBL(shift));
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_histogram_free, hnew);
}

static VALUE rb_gsl_histogram_fwrite(VALUE obj, VALUE io)
{
  gsl_histogram *h = NULL;
  FILE *f;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_histogram, h);
  f = rb_gsl_open_writefile(io, &flag);
  status = gsl_histogram_fwrite(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram_fread(VALUE obj, VALUE io)
{
  gsl_histogram *h = NULL;
  FILE *f;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_histogram, h);
  f = rb_gsl_open_readfile(io, &flag);
  status = gsl_histogram_fread(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram_fprintf(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h = NULL;
  FILE *fp;
  int status, flag = 0;

  if (argc != 1 && argc != 3) {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 3)", argc);
  }
  Data_Get_Struct(obj, gsl_histogram, h);
  fp = rb_gsl_open_writefile(argv[0], &flag);
  if (argc == 3) {
    Check_Type(argv[1], T_STRING);
    Check_Type(argv[2], T_STRING);
    status = gsl_histogram_fprintf(fp, h, STR2CSTR(argv[1]), STR2CSTR(argv[2]));
  } else {
    status = gsl_histogram_fprintf(fp, h, "%g", "%g");
  }
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram_printf(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h = NULL;
  int status;
  Data_Get_Struct(obj, gsl_histogram, h);
  if (argc == 2) {
    Check_Type(argv[0], T_STRING);
    Check_Type(argv[1], T_STRING);
    status = gsl_histogram_fprintf(stdout, h, STR2CSTR(argv[0]), STR2CSTR(argv[1]));
  } else {
    status = gsl_histogram_fprintf(stdout, h, "%g", "%g");
  }
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram_fscanf(VALUE obj, VALUE io)
{
  gsl_histogram *h = NULL;
  FILE *fp;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_histogram, h);
  fp = rb_gsl_open_readfile(io, &flag);
  status = gsl_histogram_fscanf(fp, h);
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram_print(VALUE obj)
{
  gsl_histogram *h = NULL;
  int status;
  Data_Get_Struct(obj, gsl_histogram, h);
  status = gsl_histogram_fprintf(stdout, h, "%g", "%g");
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram_pdf_alloc(VALUE klass, VALUE nn)
{
  gsl_histogram_pdf *h = NULL;
  gsl_histogram *h0 = NULL;
#ifdef GSL_0_9_4_LATER
	if (rb_obj_is_kind_of(nn, cgsl_histogram)) {
		Data_Get_Struct(nn, gsl_histogram, h0);
	  h = gsl_histogram_pdf_alloc(h0->n);
	  gsl_histogram_pdf_init(h, h0);
	} else {
	  CHECK_FIXNUM(nn);
    h = gsl_histogram_pdf_alloc(FIX2INT(nn));
  }
#else
  gsl_histogram *hh = NULL;
  Data_Get_Struct(nn, gsl_histogram, hh);
  h = gsl_histogram_pdf_alloc(hh);
#endif
  return Data_Wrap_Struct(klass, 0, gsl_histogram_pdf_free, h);
}

#ifdef GSL_0_9_4_LATER
static VALUE rb_gsl_histogram_pdf_init(VALUE obj, VALUE hh)
{
  gsl_histogram_pdf *p = NULL;
  gsl_histogram *h = NULL;
  CHECK_HISTOGRAM(hh);
  Data_Get_Struct(obj, gsl_histogram_pdf, p);
  Data_Get_Struct(hh, gsl_histogram, h);
  gsl_histogram_pdf_init(p, h);
  return obj;
}
#endif

static VALUE rb_gsl_histogram_pdf_sample(VALUE obj, VALUE r)
{
  gsl_histogram_pdf *p = NULL;
  Need_Float(r);
  Data_Get_Struct(obj, gsl_histogram_pdf, p);
  return rb_float_new(gsl_histogram_pdf_sample(p, NUM2DBL(r)));
}

static VALUE rb_gsl_histogram_pdf_range(VALUE obj)
{
  gsl_histogram_pdf *h = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_histogram_pdf, h);
  v = gsl_vector_view_alloc(h->n);
  v->vector.data = h->range;
  v->vector.size = h->n + 1;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_histogram_range, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_histogram_pdf_sum(VALUE obj)
{
  gsl_histogram_pdf *h = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_histogram_pdf, h);
  v = gsl_vector_view_alloc(h->n);
  v->vector.data = h->sum;
  v->vector.size = h->n + 1;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_histogram_graph(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_GNU_GRAPH
  gsl_histogram *v = NULL;
  FILE *fp = NULL;
  size_t i;
  char command[1024];
  Data_Get_Struct(obj, gsl_histogram, v);
  switch (argc) {
  case 0:
    strcpy(command, "graph -T X -g 3");
    break;
  case 1:
    make_graphcommand(command, argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    break;
  }
  fp = popen(command, "w");
  if (fp == NULL) rb_raise(rb_eIOError, "GNU graph not found.");
  for (i = 0; i < v->n; i++) {
    fprintf(fp, "%e %e\n%e %e\n", v->range[i], v->bin[i], v->range[i+1], v->bin[i]);
  }
  fflush(fp);
 pclose(fp);
  fp = NULL;
  return Qtrue;
#else
  rb_raise(rb_eNoMethodError, "not implemented");
  return Qfalse;
#endif
}

static VALUE rb_gsl_histogram_plot(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_GNU_GRAPH
  gsl_histogram *v = NULL;
  FILE *fp = NULL;
  size_t i;
  Data_Get_Struct(obj, gsl_histogram, v);
  switch (argc) {
  case 0:
    fp = popen("gnuplot -persist", "w");
    if (fp == NULL) rb_raise(rb_eIOError, "GNU graph not found.");
    fprintf(fp, "plot '-' with fsteps\n");
    break;
  case 1:
    fp = popen("gnuplot -persist", "w");
    if (fp == NULL) rb_raise(rb_eIOError, "GNU graph not found.");
    if (TYPE(argv[0]) == T_STRING) 
      fprintf(fp, "plot '-' %s\n", STR2CSTR(argv[0]));    
    else
      fprintf(fp, "plot '-' with fsteps\n");
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    break;
  }
  for (i = 0; i < v->n; i++) {
    fprintf(fp, "%e %e\n", v->range[i], v->bin[i]);
  }
  fprintf(fp, "e\n");
  fflush(fp);
  pclose(fp);
  fp = NULL;
  return Qtrue;
#else
  rb_raise(rb_eNoMethodError, "not implemented");
  return Qfalse;
#endif
}

struct fit_histogram {
  gsl_histogram *h;
  size_t binstart, binend;
};

static VALUE rb_gsl_histogram_fit_exponential(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h;
  gsl_vector *x, *lny, *w;
  size_t binstart = 0, binend, n, p = 2, dof, i;
  double c0, c1, cov00, cov01, cov11, sumsq, xl, xh;
  Data_Get_Struct(obj, gsl_histogram, h);
  binstart = 0;
  binend = h->n - 1;
  switch (argc) {
  case 2:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    binstart = FIX2INT(argv[0]);
    binend = FIX2INT(argv[1]);
    if (binend >= h->n) binend = h->n - 1;
    break;
  case 0:
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments (%d for 0 or 2)", argc);
    break;
  }
  n = binend - binstart + 1;
  dof = n - p;
  
  x = gsl_vector_alloc(n);
  w = gsl_vector_alloc(n);
  lny = gsl_vector_alloc(n);
  for (i = 0; i < n; i++) {
    if (gsl_histogram_get_range(h, i+binstart, &xl, &xh))
      rb_raise(rb_eIndexError, "wrong index");
    gsl_vector_set(x, i, (xl+xh)/2.0);
    gsl_vector_set(lny, i, log(h->bin[i+binstart]));
    gsl_vector_set(w, i, h->bin[i+binstart]);
  }
  gsl_fit_wlinear(x->data, 1, w->data, 1, lny->data, 1, n, 
		  &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
  gsl_vector_free(lny);
  gsl_vector_free(w);
  gsl_vector_free(x);
  c0 = exp(c0);
  return rb_ary_new3(6, rb_float_new(c0), rb_float_new(c1), 
		     rb_float_new(c0*sqrt(cov00)), rb_float_new(sqrt(cov11)),
		     rb_float_new(sumsq), INT2FIX(dof));
}

static VALUE rb_gsl_histogram_fit_power(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h;
  gsl_vector *lnx, *lny, *w;
  size_t binstart = 0, binend, n, p = 2, dof, i;
  double c0, c1, cov00, cov01, cov11, sumsq, xl, xh;
  Data_Get_Struct(obj, gsl_histogram, h);
  binstart = 0;
  binend = h->n - 1;
  switch (argc) {
  case 2:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    binstart = FIX2INT(argv[0]);
    binend = FIX2INT(argv[1]);
    if (binend >= h->n) binend = h->n - 1;
    break;
  case 0:
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments (%d for 0 or 2)", argc);
    break;
  }
  n = binend - binstart + 1;
  dof = n - p;
  
  lnx = gsl_vector_alloc(n);
  w = gsl_vector_alloc(n);
  lny = gsl_vector_alloc(n);
  for (i = 0; i < n; i++) {
    if (gsl_histogram_get_range(h, i+binstart, &xl, &xh))
      rb_raise(rb_eIndexError, "wrong index");
    gsl_vector_set(lnx, i, (log(xl)+log(xh))/2.0);
    gsl_vector_set(lny, i, log(h->bin[i+binstart]));
    gsl_vector_set(w, i, h->bin[i+binstart]);
  }
  gsl_fit_wlinear(lnx->data, 1, w->data, 1, lny->data, 1, n, 
		  &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
  gsl_vector_free(lny);
  gsl_vector_free(w);
  gsl_vector_free(lnx);
  c0 = exp(c0);
  return rb_ary_new3(6, rb_float_new(c0), rb_float_new(c1), 
		     rb_float_new(c0*sqrt(cov00)), rb_float_new(sqrt(cov11)),
		     rb_float_new(sumsq), INT2FIX(dof));
}

static int Gaussian_f(const gsl_vector *v, void *params, gsl_vector *f);
static int Gaussian_df(const gsl_vector *v, void *params, gsl_matrix * J);
static int Gaussian_f(const gsl_vector *v, void *params, gsl_vector *f)
{
  struct fit_histogram *hh;
  gsl_histogram *h = NULL;
  double amp, mu, var, xl, xh, xi, yi, sqw;
  size_t i, binstart, binend;
  hh = (struct fit_histogram *) params;
  h = hh->h;
  binstart = hh->binstart;
  binend = hh->binend;
  var = gsl_vector_get(v, 0);
  mu = gsl_vector_get(v, 1);
  amp = gsl_vector_get(v, 2);
  for (i = binstart; i <= binend; i++) {
    if (gsl_histogram_get_range(h, i, &xl, &xh))
      rb_raise(rb_eIndexError, "wrong index");
    xi = (xl + xh)/2.0;
    yi = h->bin[i];
    //    sqw = sqrt(yi);
    if (yi >= 1.0) sqw = 1.0/sqrt(yi);
    else sqw = 1.0;
    gsl_vector_set(f, i-binstart, (amp*exp(-(xi - mu)*(xi - mu)/var/2.0) - yi)*sqw);
  }
  return GSL_SUCCESS;
}

static int Gaussian_df(const gsl_vector *v, void *params, gsl_matrix *J)
{
  struct fit_histogram *hh;
  gsl_histogram *h = NULL;
  double amp, mu, var, xl, xh, xi, yi, y, sqw;
  size_t i, binstart, binend;
  hh = (struct fit_histogram *) params;
  h = hh->h;
  binstart = hh->binstart;
  binend = hh->binend;
  var = gsl_vector_get(v, 0);
  mu = gsl_vector_get(v, 1);
  amp = gsl_vector_get(v, 2);
  for (i = binstart; i <= binend; i++) {
    if (gsl_histogram_get_range(h, i, &xl, &xh))
      rb_raise(rb_eIndexError, "wrong index");
    xi = (xl + xh)/2.0;
    yi = h->bin[i];
    //    sqw = sqrt(yi);
    if (yi >= 1.0) sqw = 1.0/sqrt(yi);
    else sqw = 1.0;
    y = exp(-(xi - mu)*(xi - mu)/var/2.0);
    gsl_matrix_set(J, i-binstart, 0, amp*y*(xi - mu)*(xi - mu)/2/var/var*sqw);
    gsl_matrix_set(J, i-binstart, 1, amp*y*(xi - mu)/var*sqw);
    gsl_matrix_set(J, i-binstart, 2, y*sqw);
  }
  return GSL_SUCCESS;
}

static int Gaussian_fdf(const gsl_vector *v, void *params, gsl_vector *f, 
			gsl_matrix *J)
{
  Gaussian_f(v, params, f);
  Gaussian_df(v, params, J);
  return GSL_SUCCESS;
}

static VALUE rb_gsl_histogram_fit_gaussian(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h = NULL;
  struct fit_histogram hh;
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  size_t iter = 0, binstart, binend;
  size_t n, dof;      /* # of data points */
  size_t p = 3;  /* # of fitting parameters */
  gsl_multifit_function_fdf f;
  gsl_matrix *covar = NULL;
  gsl_vector *x = NULL;
  double sigma, mean, height, errs, errm, errh, chi2;
  Data_Get_Struct(obj, gsl_histogram, h);
  binstart = 0;
  binend = h->n - 1;
  switch (argc) {
  case 2:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    binstart = FIX2INT(argv[0]);
    binend = FIX2INT(argv[1]);
    if (binend >= h->n) binend = h->n - 1;
    break;
  case 0:
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments (%d for 0 or 2)", argc);
    break;
  }
  x = gsl_vector_alloc(p);
  gsl_vector_set(x, 0, gsl_pow_2(gsl_histogram_sigma(h)));   /* initial values, var = 1 */
  gsl_vector_set(x, 1, gsl_histogram_mean(h));   /* mu = 0 */
  gsl_vector_set(x, 2, gsl_histogram_max_val(h));   /* amp = 1 */
  hh.h = h;
  hh.binstart = binstart;
  hh.binend = binend;
  n = binend - binstart + 1;

  covar = gsl_matrix_alloc(p, p);

  f.f = Gaussian_f;
  f.df = Gaussian_df;
  f.fdf = Gaussian_fdf;
  f.n = n;
  f.p = p;
  f.params = &hh;

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, n, p);
  gsl_multifit_fdfsolver_set(s, &f, x);

  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(s);
    if (status) break;
    status = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
  } while (status == GSL_CONTINUE);
  sigma = sqrt(gsl_vector_get(s->x, 0));
  mean = gsl_vector_get(s->x, 1);
  height = gsl_vector_get(s->x, 2)*sigma*sqrt(2*M_PI);
  gsl_multifit_covar(s->J, 0.0, covar);
  chi2 = gsl_pow_2(gsl_blas_dnrm2(s->f));   /* not reduced chi-square */
  dof = n - p;
  errs = sqrt(chi2/dof*gsl_matrix_get(covar, 0, 0))/sigma/2;
  errm = sqrt(chi2/dof*gsl_matrix_get(covar, 1, 1));
  errh = sqrt(chi2/dof*gsl_matrix_get(covar, 2, 2));

  gsl_multifit_fdfsolver_free(s);
  gsl_vector_free(x);
  gsl_matrix_free(covar);
  return rb_ary_new3(8, rb_float_new(sigma), rb_float_new(mean),
		     rb_float_new(height), rb_float_new(errs),
		     rb_float_new(errm), rb_float_new(errh),
		     rb_float_new(chi2), INT2FIX(dof));
}

static int Rayleigh_f(const gsl_vector *v, void *params, gsl_vector *f);
static int Rayleigh_df(const gsl_vector *v, void *params, gsl_matrix * J);
static int Rayleigh_f(const gsl_vector *v, void *params, gsl_vector *f)
{
  struct fit_histogram *hh;
  gsl_histogram *h = NULL;
  double amp, var, xl, xh, xi, yi, sqw;
  size_t i, binstart, binend;
  hh = (struct fit_histogram *) params;
  h = hh->h;
  binstart = hh->binstart;
  binend = hh->binend;
  var = gsl_vector_get(v, 0);
  amp = gsl_vector_get(v, 1);
  for (i = binstart; i <= binend; i++) {
    if (gsl_histogram_get_range(h, i, &xl, &xh))
      rb_raise(rb_eIndexError, "wrong index");
    xi = (xl + xh)/2.0;
    yi = h->bin[i];
    sqw = sqrt(yi);
    gsl_vector_set(f, i-binstart, (amp*xi*exp(-xi*xi/var/2.0) - yi)*sqw);
  }
  return GSL_SUCCESS;
}

static int Rayleigh_df(const gsl_vector *v, void *params, gsl_matrix *J)
{
  struct fit_histogram *hh;
  gsl_histogram *h = NULL;
  double amp, var, xl, xh, xi, yi, y, sqw;
  size_t i, binstart, binend;
  hh = (struct fit_histogram *) params;
  h = hh->h;
  binstart = hh->binstart;
  binend = hh->binend;
  var = gsl_vector_get(v, 0);
  amp = gsl_vector_get(v, 1);
  for (i = binstart; i <= binend; i++) {
    if (gsl_histogram_get_range(h, i, &xl, &xh))
      rb_raise(rb_eIndexError, "wrong index");
    xi = (xl + xh)/2.0;
    yi = h->bin[i];
    sqw = sqrt(yi);
    y = xi*exp(-xi*xi/var/2.0);
    gsl_matrix_set(J, i-binstart, 0, amp*y*xi*xi/2/var/var*sqw);
    gsl_matrix_set(J, i-binstart, 1, y*sqw);
  }
  return GSL_SUCCESS;
}

static int Rayleigh_fdf(const gsl_vector *v, void *params, gsl_vector *f, 
			gsl_matrix *J)
{
  Rayleigh_f(v, params, f);
  Rayleigh_df(v, params, J);
  return GSL_SUCCESS;
}


static VALUE rb_gsl_histogram_fit_rayleigh(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h = NULL;
  struct fit_histogram hh;
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  size_t iter = 0, binstart, binend;
  size_t n, dof;      /* # of data points */
  size_t p = 2;  /* # of fitting parameters */
  gsl_multifit_function_fdf f;
  gsl_matrix *covar = NULL;
  gsl_vector *x = NULL;
  double sigma, height, errs, errh, chi2;
  Data_Get_Struct(obj, gsl_histogram, h);
  binstart = 0;
  binend = h->n - 1;
  switch (argc) {
  case 2:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    binstart = FIX2INT(argv[0]);
    binend = FIX2INT(argv[1]);
    if (binend >= h->n) binend = h->n - 1;
    break;
  case 0:
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments (%d for 0 or 2)", argc);
    break;
  }
  x = gsl_vector_alloc(p);
  gsl_vector_set(x, 0, gsl_pow_2(gsl_histogram_sigma(h)));   /* initial values, var = 1 */
  gsl_vector_set(x, 1, gsl_histogram_max_val(h));   
  hh.h = h;
  hh.binstart = binstart;
  hh.binend = binend;
  n = binend - binstart + 1;

  covar = gsl_matrix_alloc(p, p);

  f.f = Rayleigh_f;
  f.df = Rayleigh_df;
  f.fdf = Rayleigh_fdf;
  f.n = n;
  f.p = p;
  f.params = &hh;

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, n, p);
  gsl_multifit_fdfsolver_set(s, &f, x);

  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(s);
    if (status) break;
    status = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
  } while (status == GSL_CONTINUE);
  sigma = sqrt(gsl_vector_get(s->x, 0));
  height = gsl_vector_get(s->x, 1)*sigma*sigma;
  gsl_multifit_covar(s->J, 0.0, covar);
  chi2 = gsl_pow_2(gsl_blas_dnrm2(s->f));   /* not reduced chi-square */
  dof = n - p;
  errs = sqrt(chi2/dof*gsl_matrix_get(covar, 0, 0))/sigma/2;
  errh = sqrt(chi2/dof*gsl_matrix_get(covar, 1, 1));

  gsl_multifit_fdfsolver_free(s);
  gsl_vector_free(x);
  gsl_matrix_free(covar);
  return rb_ary_new3(6, rb_float_new(sigma), 
		     rb_float_new(height), rb_float_new(errs),
		     rb_float_new(errh),
		     rb_float_new(chi2), INT2FIX(dof));
}

/*
 * y(x) = Amp*exp(-a*x)
 */
static int xExponential_f(const gsl_vector *v, void *params, gsl_vector *f);
static int xExponential_df(const gsl_vector *v, void *params, gsl_matrix * J);
static int xExponential_f(const gsl_vector *v, void *params, gsl_vector *f)
{
  struct fit_histogram *hh;
  gsl_histogram *h = NULL;
  double amp, b, xl, xh, xi, yi, sqw;
  size_t i, binstart, binend;
  hh = (struct fit_histogram *) params;
  h = hh->h;
  binstart = hh->binstart;
  binend = hh->binend;
  b = gsl_vector_get(v, 0);
  amp = gsl_vector_get(v, 1);
  for (i = binstart; i <= binend; i++) {
    if (gsl_histogram_get_range(h, i, &xl, &xh))
      rb_raise(rb_eIndexError, "wrong index");
    xi = (xl + xh)/2.0;
    yi = h->bin[i];
    sqw = sqrt(yi);
    gsl_vector_set(f, i-binstart, (amp*xi*exp(-b*xi) - yi)*sqw);
  }
  return GSL_SUCCESS;
}

static int xExponential_df(const gsl_vector *v, void *params, gsl_matrix *J)
{
  struct fit_histogram *hh;
  gsl_histogram *h = NULL;
  double amp, b, xl, xh, xi, yi, y, sqw;
  size_t i, binstart, binend;
  hh = (struct fit_histogram *) params;
  h = hh->h;
  binstart = hh->binstart;
  binend = hh->binend;
  b = gsl_vector_get(v, 0);
  amp = gsl_vector_get(v, 1);
  for (i = binstart; i <= binend; i++) {
    if (gsl_histogram_get_range(h, i, &xl, &xh))
      rb_raise(rb_eIndexError, "wrong index");
    xi = (xl + xh)/2.0;
    yi = h->bin[i];
    sqw = sqrt(yi);
    y = xi*exp(-b*xi);
    gsl_matrix_set(J, i-binstart, 0, -amp*y*xi*sqw);
    gsl_matrix_set(J, i-binstart, 1, y*sqw);
  }
  return GSL_SUCCESS;
}

static int xExponential_fdf(const gsl_vector *v, void *params, gsl_vector *f, 
			gsl_matrix *J)
{
  xExponential_f(v, params, f);
  xExponential_df(v, params, J);
  return GSL_SUCCESS;
}


static VALUE rb_gsl_histogram_fit_xexponential(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h = NULL;
  struct fit_histogram hh;
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  size_t iter = 0, binstart, binend;
  size_t n, dof;      /* # of data points */
  size_t p = 2;  /* # of fitting parameters */
  gsl_multifit_function_fdf f;
  gsl_matrix *covar = NULL;
  gsl_vector *x = NULL;
  double b, height, errs, errh, chi2;
  Data_Get_Struct(obj, gsl_histogram, h);
  binstart = 0;
  binend = h->n - 1;
  switch (argc) {
  case 2:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    binstart = FIX2INT(argv[0]);
    binend = FIX2INT(argv[1]);
    if (binend >= h->n) binend = h->n - 1;
    break;
  case 0:
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments (%d for 0 or 2)", argc);
    break;
  }
  x = gsl_vector_alloc(p);
  gsl_vector_set(x, 0, gsl_histogram_sigma(h));   /* initial values, var = 1 */
  gsl_vector_set(x, 1, gsl_histogram_max_val(h));   
  hh.h = h;
  hh.binstart = binstart;
  hh.binend = binend;
  n = binend - binstart + 1;

  covar = gsl_matrix_alloc(p, p);

  f.f = xExponential_f;
  f.df = xExponential_df;
  f.fdf = xExponential_fdf;
  f.n = n;
  f.p = p;
  f.params = &hh;

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, n, p);
  gsl_multifit_fdfsolver_set(s, &f, x);

  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(s);
    if (status) break;
    status = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
  } while (status == GSL_CONTINUE);
  b = gsl_vector_get(s->x, 0);
  height = gsl_vector_get(s->x, 1);
  gsl_multifit_covar(s->J, 0.0, covar);
  chi2 = gsl_pow_2(gsl_blas_dnrm2(s->f));   /* not reduced chi-square */
  dof = n - p;
  errs = sqrt(chi2/dof*gsl_matrix_get(covar, 0, 0));
  errh = sqrt(chi2/dof*gsl_matrix_get(covar, 1, 1));

  gsl_multifit_fdfsolver_free(s);
  gsl_vector_free(x);
  gsl_matrix_free(covar);
  return rb_ary_new3(6, rb_float_new(b), 
		     rb_float_new(height), rb_float_new(errs),
		     rb_float_new(errh),
		     rb_float_new(chi2), INT2FIX(dof));
}

static VALUE rb_gsl_histogram_fit(int argc, VALUE *argv, VALUE obj)
{
  char fittype[32];
  if (argc < 1) rb_raise(rb_eArgError, "too few arguments");
  Check_Type(argv[0], T_STRING);
  strcpy(fittype, STR2CSTR(argv[0]));
  if (str_head_grep(fittype, "exp") == 0) {
    return rb_gsl_histogram_fit_exponential(argc-1, argv+1, obj);
  } else if (str_head_grep(fittype, "power") == 0) {
    return rb_gsl_histogram_fit_power(argc-1, argv+1, obj);
  } else if (str_head_grep(fittype, "gaus") == 0) {
    return rb_gsl_histogram_fit_gaussian(argc-1, argv+1, obj);
 } else if (str_head_grep(fittype, "rayleigh") == 0) {
    return rb_gsl_histogram_fit_rayleigh(argc-1, argv+1, obj);
  } else if (str_head_grep(fittype, "xexp") == 0) {
    return rb_gsl_histogram_fit_xexponential(argc-1, argv+1, obj);
  } else {
    rb_raise(rb_eRuntimeError, 
	     "unknown fitting type %s (exp, power, gaus expected)", fittype);
  }
  return Qnil;
}

/* Integrate histogram: the two histograms must have the same range and bins. */
void mygsl_histogram_integrate(const gsl_histogram *h, gsl_histogram *hi,
			       size_t istart, size_t iend)
{
  size_t i;
  if (iend >= istart) {
    //if (istart < 0) istart = 0;
    if (iend >= h->n) iend = h->n-1;
    hi->bin[istart] = h->bin[istart];
    for (i = istart+1; i <= iend; i++) hi->bin[i] = hi->bin[i-1] + h->bin[i];
  } else {
    if (istart >= h->n) istart = h->n-1;
    //if (iend < 0) iend = 0;
    hi->bin[istart] = h->bin[istart];
    for (i = istart-1; i >= iend; i--) {
      hi->bin[i] = hi->bin[i+1] + h->bin[i];
      if (i == 0) break;
    }
  }
}

void mygsl_histogram_differentiate(const gsl_histogram *hi, gsl_histogram *h)
{
  size_t i;
  h->bin[0] = hi->bin[0];
  for (i = 1; i < hi->n; i++) h->bin[i] = hi->bin[i] - hi->bin[i-1];
}

/* Create a histogram integrating the given histogram h */
gsl_histogram* mygsl_histogram_calloc_integrate(const gsl_histogram *h,
						size_t istart, size_t iend)
{
  gsl_histogram *hi = NULL;
  hi = gsl_histogram_calloc_range(h->n, h->range);
  mygsl_histogram_integrate(h, hi, istart, iend);
  return hi;
}

gsl_histogram* mygsl_histogram_calloc_differentiate(const gsl_histogram *hi)
{
  gsl_histogram *h = NULL;
  h = gsl_histogram_calloc_range(hi->n, hi->range);
  mygsl_histogram_differentiate(hi, h);
  return h;
}

static VALUE rb_gsl_histogram_integrate(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h, *hi;
  size_t istart, iend;
  int itmp;
  Data_Get_Struct(obj, gsl_histogram, h);
  switch (argc) {
  case 2:
    istart = FIX2INT(argv[0]);
    iend = FIX2INT(argv[1]);
    break;
  case 1:
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      istart = FIX2INT(rb_ary_entry(argv[0], 0));
      iend = FIX2INT(rb_ary_entry(argv[0], 1));
      break;
    case T_FIXNUM:
      itmp = FIX2INT(argv[0]);
      if (itmp == -1) {
	istart = h->n - 1;
	iend = 0;
      } else {
	istart = 0;
	iend = h->n - 1;
      }
      break;
    default:
      rb_raise(rb_eArgError, "wrong argument type %s (Arran or Fixnum expected)",
	       rb_class2name(CLASS_OF(argv[0])));
      break;
    }
    break;
  case 0:
    istart = 0;
    iend = h->n - 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0-2)", argc);
    break;
  }
  hi = mygsl_histogram_calloc_integrate(h, istart, iend);
  return Data_Wrap_Struct(cgsl_histogram_integ, 0, gsl_histogram_free, hi);
}

static VALUE rb_gsl_histogram_differentiate(VALUE obj)
{
  gsl_histogram *h, *hi;
  Data_Get_Struct(obj, gsl_histogram, hi);
  h = mygsl_histogram_calloc_differentiate(hi);
  return Data_Wrap_Struct(cgsl_histogram, 0, gsl_histogram_free, h);
}

static gsl_histogram* mygsl_histogram_rebin(const gsl_histogram *h, size_t m)
{
  gsl_histogram *hnew;
  double w;
  size_t n, i, j, k;
  if (m > h->n) m = h->n;
  n = (size_t) h->n/m;
  if (n*m != h->n) n += 1;
  w = (h->range[h->n] - h->range[0])/h->n;
  hnew = gsl_histogram_alloc(n);
  for (i = 0, j = 0; i <= n; i++) {
    if (i*m <= h->n) hnew->range[i] = h->range[i*m];
    else hnew->range[i] = w*m*i;
  }
  for (i = 0, j = 0; i < n; i++) {
    hnew->bin[i] = 0;
    for (k = 0; k < m && j < h->n; k++) hnew->bin[i] += h->bin[j++];
  }
  return hnew;
}

static VALUE rb_gsl_histogram_rebin(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram *h, *hnew;
  size_t m = 2;
  switch (argc) {
  case 1:
    CHECK_FIXNUM(argv[0]);
    m = (size_t) FIX2INT(argv[0]);
    break;
  case 0:
    m = 2;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    break;
  }
  Data_Get_Struct(obj, gsl_histogram, h);
  hnew = mygsl_histogram_rebin(h, m);
  return Data_Wrap_Struct(cgsl_histogram, 0, gsl_histogram_free, hnew);
}

static int mygsl_histogram_fread2(FILE * stream, gsl_histogram * h)
{  
  double min, max;
  int status;
  status = gsl_block_raw_fread(stream, &min, 1, 1);
  if (status)    return status;  
  status = gsl_block_raw_fread(stream, &max, 1, 1);
  if (status)    return status;  
  gsl_histogram_set_ranges_uniform(h, min, max);
  status = gsl_block_raw_fread (stream, h->bin, h->n, 1);  
  if (status)    return status;  
  return status;
}

static int mygsl_histogram_fwrite2(FILE * stream, const gsl_histogram * h)
{  
  int status;
  status = gsl_block_raw_fwrite (stream, h->range, 1, 1);
  if (status)    return status;  
  status = gsl_block_raw_fwrite (stream, h->range+h->n, 1, 1);
  if (status)    return status;  
  status = gsl_block_raw_fwrite (stream, h->bin, h->n, 1);  
  return status;
}

static VALUE rb_gsl_histogram_fwrite2(VALUE obj, VALUE io)
{
  gsl_histogram *h = NULL;
  FILE *f;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_histogram, h);
  f = rb_gsl_open_writefile(io, &flag);
  status = mygsl_histogram_fwrite2(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram_fread2(VALUE obj, VALUE io)
{
  gsl_histogram *h = NULL;
  FILE *f;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_histogram, h);
  f = rb_gsl_open_readfile(io, &flag);
  status = mygsl_histogram_fread2(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static gsl_histogram* mygsl_histogram_calloc_reverse(const gsl_histogram *h)
{
  gsl_histogram *hnew;
  size_t i, n;
  hnew = gsl_histogram_alloc(h->n);
  n = h->n;
  for (i = 0; i <= n; i++) hnew->range[i] = h->range[n-i];
  for (i = 0; i < n; i++) hnew->bin[i] = h->bin[n-1-i];
  return hnew;
}

static VALUE rb_gsl_histogram_reverse(VALUE obj)
{
  gsl_histogram *h, *hnew;
  Data_Get_Struct(obj, gsl_histogram, h);
  hnew = mygsl_histogram_calloc_reverse(h);
  return Data_Wrap_Struct(cgsl_histogram, 0, gsl_histogram_free, hnew);
}

/* The functions below are not included in GSL */
/*
 * Returns an x value at which the histogram integration
 * reaches the given percentile. The x value is calculated
 * by an interpolation between the ranges in which the percentile
 * is found.
 */
static double histogram_percentile(const gsl_histogram *h, double f)
{
  double sum = gsl_histogram_sum(h), sf;
  double val = 0, s = 0, x;
  double ri, ri1;
  size_t i;
  sf = sum * f;
  for (i = 0; i < h->n; i++) {
    val = gsl_histogram_get(h, i);
    if ((s+val) > sf) break;
    s += val;
  }
  ri = h->range[i];
  ri1 = h->range[i+1];
  x = (sf - s)*(ri1 - ri)/val + ri;
  return x;
}

static double histogram_median(const gsl_histogram *h)
{
  return histogram_percentile(h, 0.5);
}

static VALUE rb_gsl_histogram_percentile(VALUE obj, VALUE f)
{
  gsl_histogram *h;
  Data_Get_Struct(obj, gsl_histogram, h);
  return rb_float_new(histogram_percentile(h, NUM2DBL(f)));
}

static VALUE rb_gsl_histogram_median(VALUE obj)
{
  gsl_histogram *h;
  Data_Get_Struct(obj, gsl_histogram, h);
  return rb_float_new(histogram_median(h));
}

static double histogram_percentile_inv(const gsl_histogram *h, double x)
{
  double sum = gsl_histogram_sum(h);
  double val = 0, s = 0;
  double ri, ri1, q;
  size_t i;

  for (i = 0; i < h->n; i++) {
    val = gsl_histogram_get(h, i);
    if (h->range[i+1] > x) break;
    s += val;
  }
  ri = h->range[i];
  ri1 = h->range[i+1];
  q = s + val/(ri1 - ri)*(x - ri);
  return q/sum;
}

static VALUE rb_gsl_histogram_percentile_inv(VALUE obj, VALUE x)
{
  gsl_histogram *h;
  Data_Get_Struct(obj, gsl_histogram, h);
  return rb_float_new(histogram_percentile_inv(h, NUM2DBL(x)));
}

void Init_gsl_histogram(VALUE module)
{
  VALUE cgsl_histogram_pdf;

  cgsl_histogram = rb_define_class_under(module, "Histogram", cGSL_Object);
  cgsl_histogram_range = rb_define_class_under(cgsl_histogram, "Range", 
					       cgsl_vector_view_ro);
  cgsl_histogram_bin = rb_define_class_under(cgsl_histogram, "Bin", 
					     cgsl_vector_view);
  cgsl_histogram_integ = rb_define_class_under(cgsl_histogram, "Integral", 
					       cgsl_histogram);

#ifdef GSL_0_9_4_LATER
  rb_define_singleton_method(cgsl_histogram, "alloc", rb_gsl_histogram_alloc, -1);
  /*  rb_define_singleton_method(cgsl_histogram, "new", rb_gsl_histogram_alloc, -1);*/
  rb_define_singleton_method(cgsl_histogram, "[]", rb_gsl_histogram_alloc, -1);

  rb_define_singleton_method(cgsl_histogram, "alloc_uniform", 
			     rb_gsl_histogram_alloc_uniform, -1);
  rb_define_singleton_method(cgsl_histogram, "new_uniform", 
			     rb_gsl_histogram_alloc_uniform, -1);

  rb_define_singleton_method(cgsl_histogram, "alloc_with_min_max_step", 
			     rb_gsl_histogram_alloc_with_min_max_step, 3);
  rb_define_singleton_method(cgsl_histogram, "new_with_min_max_step", 
			     rb_gsl_histogram_alloc_with_min_max_step, 3);
#endif

  rb_define_singleton_method(cgsl_histogram, "calloc", 
			     rb_gsl_histogram_calloc, 1);
  rb_define_singleton_method(cgsl_histogram, "calloc_range", 
			     rb_gsl_histogram_calloc_range, -1);

  rb_define_method(cgsl_histogram, "bins", rb_gsl_histogram_bins, 0);
  rb_define_alias(cgsl_histogram, "n", "bins");
  rb_define_alias(cgsl_histogram, "size", "bins");
  rb_define_method(cgsl_histogram, "set_ranges", rb_gsl_histogram_set_ranges, -1);
  rb_define_method(cgsl_histogram, "range", rb_gsl_histogram_range, 0);
  rb_define_method(cgsl_histogram, "bin", rb_gsl_histogram_bin, 0);
  rb_define_method(cgsl_histogram, "set_ranges_uniform", 
		   rb_gsl_histogram_set_ranges_uniform, -1);
  rb_define_singleton_method(cgsl_histogram, "memcpy", rb_gsl_histogram_memcpy, 2);
  rb_define_method(cgsl_histogram, "clone", rb_gsl_histogram_clone, 0);
  rb_define_alias(cgsl_histogram, "duplicate", "clone");
  rb_define_method(cgsl_histogram, "increment", rb_gsl_histogram_accumulate, -1);
  rb_define_alias(cgsl_histogram, "fill", "increment");
  rb_define_alias(cgsl_histogram, "accumulate", "increment");
  rb_define_method(cgsl_histogram, "increment2", rb_gsl_histogram_accumulate2, -1);
  rb_define_alias(cgsl_histogram, "accumulate2", "increment2");
  rb_define_alias(cgsl_histogram, "fill2", "increment2");

  rb_define_method(cgsl_histogram, "get", rb_gsl_histogram_get, 1);
  rb_define_alias(cgsl_histogram, "[]", "get");
  rb_define_method(cgsl_histogram, "get_range", rb_gsl_histogram_get_range, 1);
  rb_define_method(cgsl_histogram, "max", rb_gsl_histogram_max, 0);
  rb_define_method(cgsl_histogram, "min", rb_gsl_histogram_min, 0);
  rb_define_method(cgsl_histogram, "reset", rb_gsl_histogram_reset, 0);
  rb_define_method(cgsl_histogram, "find", rb_gsl_histogram_find, 1);
  rb_define_method(cgsl_histogram, "max_val", rb_gsl_histogram_max_val, 0);
  rb_define_method(cgsl_histogram, "max_bin", rb_gsl_histogram_max_bin, 0);
  rb_define_method(cgsl_histogram, "min_val", rb_gsl_histogram_min_val, 0);
  rb_define_method(cgsl_histogram, "min_bin", rb_gsl_histogram_min_bin, 0);
  rb_define_method(cgsl_histogram, "mean", rb_gsl_histogram_mean, 0);
  rb_define_method(cgsl_histogram, "sigma", rb_gsl_histogram_sigma, 0);

  rb_define_method(cgsl_histogram, "sum", rb_gsl_histogram_integral, -1);
  rb_define_alias(cgsl_histogram, "integral", "sum");

  rb_define_method(cgsl_histogram, "equal_bins_p", 
		   rb_gsl_histogram_equal_bins_p, -1);
  rb_define_alias(cgsl_histogram, "equal_bins", "equal_bins_p");
  rb_define_singleton_method(cgsl_histogram, "equal_bins_p", 
			     rb_gsl_histogram_equal_bins_p, -1);
  rb_define_singleton_method(cgsl_histogram, "equal_bins", 
			     rb_gsl_histogram_equal_bins_p, -1);
  rb_define_method(cgsl_histogram, "equal_bins_p?", 
		   rb_gsl_histogram_equal_bins_p2, -1);
  rb_define_alias(cgsl_histogram, "equal_bins?", "equal_bins_p?");
  rb_define_singleton_method(cgsl_histogram, "equal_bins_p?", 
			     rb_gsl_histogram_equal_bins_p2, -1);
  rb_define_singleton_method(cgsl_histogram, "equal_bins?", 
			     rb_gsl_histogram_equal_bins_p2, -1);

  rb_define_method(cgsl_histogram, "add", rb_gsl_histogram_add, 1);
  rb_define_alias(cgsl_histogram, "+", "add");
  rb_define_method(cgsl_histogram, "sub", rb_gsl_histogram_sub, 1);
  rb_define_alias(cgsl_histogram, "-", "sub");
  rb_define_method(cgsl_histogram, "mul", rb_gsl_histogram_mul, 1);
  rb_define_alias(cgsl_histogram, "*", "mul");
  rb_define_method(cgsl_histogram, "div", rb_gsl_histogram_div, 1);
  rb_define_alias(cgsl_histogram, "/", "div");

  rb_define_method(cgsl_histogram, "add!", rb_gsl_histogram_add2, 1);
  rb_define_method(cgsl_histogram, "sub!", rb_gsl_histogram_sub2, 1);
  rb_define_method(cgsl_histogram, "mul!", rb_gsl_histogram_mul2, 1);
  rb_define_method(cgsl_histogram, "div!", rb_gsl_histogram_div2, 1);

  rb_define_method(cgsl_histogram, "scale!", rb_gsl_histogram_scale_bang, -1);
  rb_define_method(cgsl_histogram, "scale", rb_gsl_histogram_scale, -1);
  rb_define_method(cgsl_histogram, "shift!", rb_gsl_histogram_shift, 1);
  rb_define_method(cgsl_histogram, "shift", rb_gsl_histogram_shift2, 1);

  rb_define_method(cgsl_histogram, "fwrite", rb_gsl_histogram_fwrite, 1);
  rb_define_method(cgsl_histogram, "fread", rb_gsl_histogram_fread, 1);
  rb_define_method(cgsl_histogram, "fwrite2", rb_gsl_histogram_fwrite2, 1);
  rb_define_method(cgsl_histogram, "fread2", rb_gsl_histogram_fread2, 1);
  rb_define_method(cgsl_histogram, "fprintf", rb_gsl_histogram_fprintf, -1);
  rb_define_method(cgsl_histogram, "printf", rb_gsl_histogram_printf, -1);
  rb_define_method(cgsl_histogram, "fscanf", rb_gsl_histogram_fscanf, 1);
  rb_define_method(cgsl_histogram, "print", rb_gsl_histogram_print, 0);

  cgsl_histogram_pdf = rb_define_class_under(cgsl_histogram, "Pdf", cGSL_Object);
  rb_define_singleton_method(cgsl_histogram_pdf, "alloc", 
			     rb_gsl_histogram_pdf_alloc, 1);
  /*  rb_define_singleton_method(cgsl_histogram_pdf, "new", 
      rb_gsl_histogram_pdf_alloc, 1);*/
#ifdef GSL_0_9_4_LATER
  rb_define_method(cgsl_histogram_pdf, "init", rb_gsl_histogram_pdf_init, 1);
#endif
  rb_define_method(cgsl_histogram_pdf, "sample", rb_gsl_histogram_pdf_sample, 1);

  rb_define_method(cgsl_histogram_pdf, "range", rb_gsl_histogram_pdf_range, 0);
  rb_define_method(cgsl_histogram_pdf, "sum", rb_gsl_histogram_pdf_sum, 0);

  /*****/
  rb_define_method(cgsl_histogram, "graph", rb_gsl_histogram_graph, -1);
  rb_define_alias(cgsl_histogram, "draw", "graph");

  rb_define_method(cgsl_histogram, "plot", rb_gsl_histogram_plot, -1);

  rb_define_method(cgsl_histogram, "fit_gaussian", 
		   rb_gsl_histogram_fit_gaussian, -1);
  rb_define_method(cgsl_histogram, "fit_exponential", 
		   rb_gsl_histogram_fit_exponential, -1);
  rb_define_method(cgsl_histogram, "fit_xexponential", 
		   rb_gsl_histogram_fit_xexponential, -1);
  rb_define_method(cgsl_histogram, "fit_power", 
		   rb_gsl_histogram_fit_power, -1);
  rb_define_method(cgsl_histogram, "fit_rayleigh", 
		   rb_gsl_histogram_fit_rayleigh, -1);
  rb_define_method(cgsl_histogram, "fit", 
		   rb_gsl_histogram_fit, -1);

  rb_define_method(cgsl_histogram, "integrate", rb_gsl_histogram_integrate, -1);
  rb_undef_method(cgsl_histogram_integ, "integrate");
  rb_define_method(cgsl_histogram_integ, "differentiate", 
		   rb_gsl_histogram_differentiate, 0);
  rb_define_alias(cgsl_histogram_integ, "diff", "differentiate");

  rb_define_method(cgsl_histogram, "normalize", rb_gsl_histogram_normalize, 0);
  rb_define_method(cgsl_histogram, "normalize!", rb_gsl_histogram_normalize_bang, 0);

  rb_define_method(cgsl_histogram, "rebin", rb_gsl_histogram_rebin, -1);
  rb_define_alias(cgsl_histogram, "mergebin", "rebin");

  rb_define_method(cgsl_histogram, "reverse", rb_gsl_histogram_reverse, 0);

  rb_define_method(cgsl_histogram, "percentile", rb_gsl_histogram_percentile, 1);
  rb_define_method(cgsl_histogram, "median", rb_gsl_histogram_median, 0);
  rb_define_method(cgsl_histogram, "percentile_inv", rb_gsl_histogram_percentile_inv, 1);
}
