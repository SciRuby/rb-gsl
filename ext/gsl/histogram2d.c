/*
  histogram2d.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_histogram.h"
#include "include/rb_gsl_common.h"
#include "include/rb_gsl_array.h"

VALUE cgsl_histogram2d;
VALUE cgsl_histogram2d_view;
static VALUE cgsl_histogram2d_integ;

#ifdef GSL_0_9_4_LATER
static VALUE rb_gsl_histogram2d_alloc_uniform(int argc, VALUE *argv, VALUE klass);
static VALUE rb_gsl_histogram2d_alloc(int argc, VALUE *argv, VALUE klass)
{
  gsl_histogram2d *h = NULL;
  size_t xsize, ysize;
  gsl_vector *vx, *vy;
  switch (argc) {
  case 2:
    switch (TYPE(argv[0])) {
    case T_FIXNUM:
      CHECK_FIXNUM(argv[1]);
      h = gsl_histogram2d_calloc(FIX2INT(argv[0]), FIX2INT(argv[1]));
      break;
    case T_ARRAY:
      vx = make_cvector_from_rarray(argv[0]);
      vy = make_cvector_from_rarray(argv[1]);
      h = gsl_histogram2d_alloc(vx->size-1, vy->size-1);
      gsl_histogram2d_set_ranges(h, vx->data, vx->size, vy->data, vy->size);
      gsl_vector_free(vx);
      gsl_vector_free(vy);
      break;
    default:
      CHECK_VECTOR(argv[0]); CHECK_VECTOR(argv[1]);
      Data_Get_Struct(argv[0], gsl_vector, vx);
      Data_Get_Struct(argv[1], gsl_vector, vy);
      h = gsl_histogram2d_alloc(vx->size-1, vy->size-1);
      gsl_histogram2d_set_ranges(h, vx->data, vx->size, vy->data, vy->size);
      break;
    }
    return Data_Wrap_Struct(klass, 0, gsl_histogram2d_free, h);
    break;
  case 4:
    if (VECTOR_P(argv[0]) && VECTOR_P(argv[2])) {
      CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[3]);
      Data_Get_Struct(argv[0], gsl_vector, vx);
      Data_Get_Struct(argv[2], gsl_vector, vy);
      xsize = (size_t) FIX2UINT(argv[1]); ysize = (size_t) FIX2UINT(argv[3]);
      h = gsl_histogram2d_alloc(xsize-1, ysize-1);
      gsl_histogram2d_set_ranges(h, vx->data, xsize, vy->data, ysize);
      return Data_Wrap_Struct(klass, 0, gsl_histogram2d_free, h);
    } else {
      return rb_gsl_histogram2d_alloc_uniform(argc, argv, klass);
    }
    break;
  case 6:
    return rb_gsl_histogram2d_alloc_uniform(argc, argv, klass);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments %d", argc);
    break;
  }
  return Qnil;  /* never reach here */
}

static VALUE rb_gsl_histogram2d_alloc_uniform(int argc, VALUE *argv, VALUE klass)
{
  gsl_histogram2d *h = NULL;
  double xmin, xmax, ymin, ymax;
  switch (argc) {
  case 4:
    CHECK_FIXNUM(argv[0]);     CHECK_FIXNUM(argv[2]); 
    Check_Type(argv[1], T_ARRAY); Check_Type(argv[3], T_ARRAY);
    //    if (RARRAY(argv[1])->len != 2 || RARRAY(argv[3])->len != 2)
    if (RARRAY_LEN(argv[1]) != 2 || RARRAY_LEN(argv[3]) != 2)
      rb_raise(rb_eArgError, "array size must be 2");
    xmin = NUM2DBL(rb_ary_entry(argv[1], 0));
    xmax = NUM2DBL(rb_ary_entry(argv[1], 1));
    ymin = NUM2DBL(rb_ary_entry(argv[3], 0));
    ymax = NUM2DBL(rb_ary_entry(argv[3], 1));
    h = gsl_histogram2d_alloc(FIX2INT(argv[0]), FIX2INT(argv[2]));
    gsl_histogram2d_set_ranges_uniform(h, xmin, xmax, ymin, ymax);
    return Data_Wrap_Struct(klass, 0, gsl_histogram2d_free, h);
    break;
  case 6:
    CHECK_FIXNUM(argv[0]); 
    Need_Float(argv[1]); Need_Float(argv[2]);
    CHECK_FIXNUM(argv[3]);
    Need_Float(argv[4]); Need_Float(argv[5]);
    h = gsl_histogram2d_alloc(FIX2INT(argv[0]), FIX2INT(argv[3]));
    gsl_histogram2d_set_ranges_uniform(h, NUM2DBL(argv[1]), NUM2DBL(argv[2]),
               NUM2DBL(argv[4]), NUM2DBL(argv[5]));
    return Data_Wrap_Struct(klass, 0, gsl_histogram2d_free, h);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments %d", argc);
    break;
  }
  return Qnil;  /* never reach here */
}
#endif

static VALUE rb_gsl_histogram2d_set_ranges(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram2d *h = NULL;
  gsl_vector *vx, *vy;
  size_t xsize, ysize;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  switch (argc) {
  case 4:
    CHECK_VECTOR(argv[0]); CHECK_VECTOR(argv[2]);
    CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[3]);
    Data_Get_Struct(argv[0], gsl_vector, vx);
    Data_Get_Struct(argv[2], gsl_vector, vy);
    xsize = (size_t) FIX2UINT(argv[1]); ysize = (size_t) FIX2UINT(argv[3]);
    break;
  case 2:
    CHECK_VECTOR(argv[0]); CHECK_VECTOR(argv[1]);
    Data_Get_Struct(argv[0], gsl_vector, vx);
    Data_Get_Struct(argv[1], gsl_vector, vy);
    xsize = vx->size; ysize = vy->size;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 4)", argc);
  }
  gsl_histogram2d_set_ranges(h, vx->data, xsize, vy->data, ysize);
  return obj;
}

static VALUE rb_gsl_histogram2d_set_ranges_uniform(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram2d *h = NULL;
  double xmin, xmax, ymin, ymax;
  switch (argc) {
  case 2:
    Check_Type(argv[0], T_ARRAY); Check_Type(argv[1], T_ARRAY);
    xmin = NUM2DBL(rb_ary_entry(argv[0], 0));
    xmax = NUM2DBL(rb_ary_entry(argv[0], 1));
    ymin = NUM2DBL(rb_ary_entry(argv[1], 0));
    ymax = NUM2DBL(rb_ary_entry(argv[1], 1));
    break;
  case 4:
    xmin = NUM2DBL(argv[0]); xmax = NUM2DBL(argv[1]);
    ymin = NUM2DBL(argv[2]); ymax = NUM2DBL(argv[3]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 4)", argc);
    break;
  }
  Data_Get_Struct(obj, gsl_histogram2d, h);
  gsl_histogram2d_set_ranges_uniform(h, xmin, xmax, ymin, ymax);
  return obj;
}

static VALUE rb_gsl_histogram2d_clone(VALUE obj)
{
  gsl_histogram2d *h, *h2 = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  h2 = gsl_histogram2d_clone(h);
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_histogram2d_free, h2);
}

/* singleton */
static VALUE rb_gsl_histogram2d_memcpy(VALUE obj, VALUE vhdest, VALUE vhsrc)
{
  gsl_histogram2d *hdest, *hsrc;
  CHECK_HISTOGRAM2D(vhdest); CHECK_HISTOGRAM2D(vhsrc);
  Data_Get_Struct(vhdest, gsl_histogram2d, hdest);
  Data_Get_Struct(vhsrc, gsl_histogram2d, hsrc);
  gsl_histogram2d_memcpy(hdest, hsrc);
  return vhdest;
}

static VALUE rb_gsl_histogram2d_accumulate(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram2d *h = NULL;
  gsl_vector *vx, *vy;
  size_t n, i;
  double weight = 1;
  switch (argc) {
  case 3:
    Need_Float(argv[2]);
    weight = NUM2DBL(argv[2]);
    break;
  case 2:
    weight = 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  Data_Get_Struct(obj, gsl_histogram2d, h);
  if (VECTOR_P(argv[0]) && VECTOR_P(argv[1])) {
    Data_Get_Struct(argv[0], gsl_vector, vx);
    Data_Get_Struct(argv[1], gsl_vector, vy);
    n = (size_t) GSL_MIN_INT((int) vx->size, (int) vy->size);
    for (i = 0; i < n; i++)
      gsl_histogram2d_accumulate(h, gsl_vector_get(vx, i), gsl_vector_get(vy, i),
         weight);
  } else {
    gsl_histogram2d_accumulate(h, NUM2DBL(argv[0]), NUM2DBL(argv[1]), weight);
  }
  return obj;
}

static VALUE rb_gsl_histogram2d_accumulate2(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram2d *h = NULL;
  double x, y, weight=1;
  switch (argc) {
  case 3:
    Need_Float(argv[2]);
    weight = NUM2DBL(argv[2]);
    /* no break */
  case 2:
    Need_Float(argv[0]); Need_Float(argv[1]); 
    x = NUM2DBL(argv[0]);    y = NUM2DBL(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  Data_Get_Struct(obj, gsl_histogram2d, h);
  if (x < h->xrange[0]) x = h->xrange[0] + 4*GSL_DBL_EPSILON;
  if (x > h->xrange[h->nx]) x = h->xrange[h->nx] - 4*GSL_DBL_EPSILON;
  if (y < h->yrange[0]) y = h->yrange[0] + 4*GSL_DBL_EPSILON;
  if (y > h->yrange[h->ny]) y = h->yrange[h->ny] - 4*GSL_DBL_EPSILON;
  gsl_histogram2d_accumulate(h, x, y, weight);
  return obj;
}

static VALUE rb_gsl_histogram2d_get(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram2d *h2 = NULL;
  mygsl_histogram2d_view *h1 = NULL;
  size_t i;
  switch (argc) {
  case 2:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    Data_Get_Struct(obj, gsl_histogram2d, h2);
    return rb_float_new(gsl_histogram2d_get(h2, (size_t) FIX2INT(argv[0]), (size_t) FIX2INT(argv[1])));
    break;
  case 1:
    Data_Get_Struct(obj, gsl_histogram2d, h2);
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      return rb_float_new(gsl_histogram2d_get(h2, FIX2INT(rb_ary_entry(argv[0], 0)),
                FIX2INT(rb_ary_entry(argv[0], 1))));
      break;
    case T_FIXNUM:
      CHECK_FIXNUM(argv[0]); 
      i = (size_t) FIX2INT(argv[0]);
      if (i >= h2->ny)
  rb_raise(rb_eIndexError, "wrong index");
      h1 = ALLOC(mygsl_histogram2d_view);
      h1->h.n = h2->ny;
      h1->h.range = h2->yrange;
      h1->h.bin = h2->bin + i*h2->ny;
      return Data_Wrap_Struct(cgsl_histogram2d_view, 0, free, h1);
      break;
    default:
      rb_raise(rb_eTypeError, "wrong argument type %s (Array or Fixnum expected)",
         rb_class2name(CLASS_OF(argv[0])));
      break;
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 1)",argc);
    break;
  }
}

static VALUE rb_gsl_histogram2d_get_xrange(VALUE obj, VALUE i)
{
  gsl_histogram2d *h = NULL;
  double x1, x2;
  CHECK_FIXNUM(i);
  Data_Get_Struct(obj, gsl_histogram2d, h);
  gsl_histogram2d_get_xrange(h, FIX2INT(i), &x1, &x2);
  return rb_ary_new3(2, rb_float_new(x1), rb_float_new(x2));
}

static VALUE rb_gsl_histogram2d_get_yrange(VALUE obj, VALUE j)
{
  gsl_histogram2d *h = NULL;
  double y1, y2;
  CHECK_FIXNUM(j);
  Data_Get_Struct(obj, gsl_histogram2d, h);
  gsl_histogram2d_get_yrange(h, FIX2INT(j), &y1, &y2);
  return rb_ary_new3(2, rb_float_new(y1), rb_float_new(y2));
}

static VALUE rb_gsl_histogram2d_xmax(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return rb_float_new(gsl_histogram2d_xmax(h));
}

static VALUE rb_gsl_histogram2d_xmin(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return rb_float_new(gsl_histogram2d_xmin(h));
}

static VALUE rb_gsl_histogram2d_nx(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return INT2FIX(gsl_histogram2d_nx(h));
}

static VALUE rb_gsl_histogram2d_ymax(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return rb_float_new(gsl_histogram2d_ymax(h));
}

static VALUE rb_gsl_histogram2d_ymin(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return rb_float_new(gsl_histogram2d_ymin(h));
}

static VALUE rb_gsl_histogram2d_ny(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return INT2FIX(gsl_histogram2d_ny(h));
}

static VALUE rb_gsl_histogram2d_find(VALUE obj, VALUE x, VALUE y)
{
  gsl_histogram2d *h = NULL;
  size_t i, j;
  Need_Float(x);Need_Float(y);
  Data_Get_Struct(obj, gsl_histogram2d, h);
  gsl_histogram2d_find(h, NUM2DBL(x), NUM2DBL(y), &i, &j);
  return rb_ary_new3(2, INT2FIX(i), INT2FIX(j));
}

static VALUE rb_gsl_histogram2d_max_val(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return rb_float_new(gsl_histogram2d_max_val(h));
}

static VALUE rb_gsl_histogram2d_max_bin(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  size_t i, j;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  gsl_histogram2d_max_bin(h, &i, &j);
  return rb_ary_new3(2, INT2FIX(i), INT2FIX(j));
}

static VALUE rb_gsl_histogram2d_min_val(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return rb_float_new(gsl_histogram2d_min_val(h));
}

static VALUE rb_gsl_histogram2d_min_bin(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  size_t i, j;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  gsl_histogram2d_min_bin(h, &i, &j);
  return rb_ary_new3(2, INT2FIX(i), INT2FIX(j));
}

#ifdef GSL_1_1_LATER
static VALUE rb_gsl_histogram2d_xmean(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return rb_float_new(gsl_histogram2d_xmean(h));
}

static VALUE rb_gsl_histogram2d_ymean(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return rb_float_new(gsl_histogram2d_ymean(h));
}

static VALUE rb_gsl_histogram2d_xsigma(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return rb_float_new(gsl_histogram2d_xsigma(h));
}

static VALUE rb_gsl_histogram2d_ysigma(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return rb_float_new(gsl_histogram2d_ysigma(h));
}

static VALUE rb_gsl_histogram2d_cov(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return rb_float_new(gsl_histogram2d_cov(h));
}

static VALUE rb_gsl_histogram2d_sum(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  return rb_float_new(gsl_histogram2d_sum(h));
}
#endif

/* singleton */
static VALUE rb_gsl_histogram2d_equal_bins_p(VALUE obj, VALUE hh1, VALUE hh2)
{
  gsl_histogram2d *h1 = NULL, *h2 = NULL;
  CHECK_HISTOGRAM2D(hh1); CHECK_HISTOGRAM2D(hh2);
  Data_Get_Struct(hh1, gsl_histogram2d, h1);
  Data_Get_Struct(hh2, gsl_histogram2d, h2);
  return INT2FIX(gsl_histogram2d_equal_bins_p(h1, h2));
}

static VALUE rb_gsl_histogram2d_equal_bins_p2(VALUE obj, VALUE hh1, VALUE hh2)
{
  gsl_histogram2d *h1 = NULL, *h2 = NULL;
  CHECK_HISTOGRAM2D(hh1); CHECK_HISTOGRAM2D(hh2);
  Data_Get_Struct(hh1, gsl_histogram2d, h1);
  Data_Get_Struct(hh2, gsl_histogram2d, h2);
  if (gsl_histogram2d_equal_bins_p(h1, h2)) return Qtrue;
  else return Qfalse;
}

static VALUE rb_gsl_histogram2d_scale(VALUE obj, VALUE s)
{
  gsl_histogram2d *h = NULL;
  Need_Float(s);
  Data_Get_Struct(obj, gsl_histogram2d, h);
  gsl_histogram2d_scale(h, NUM2DBL(s));
  return obj;
}

static VALUE rb_gsl_histogram2d_shift(VALUE obj, VALUE s)
{
  gsl_histogram2d *h = NULL;
  Need_Float(s);
  Data_Get_Struct(obj, gsl_histogram2d, h);
  gsl_histogram2d_shift(h, NUM2DBL(s));
  return obj;
}

/*****/
static VALUE rb_gsl_histogram2d_add(VALUE obj, VALUE hh2)
{
  gsl_histogram2d *h1 = NULL, *h2 = NULL, *hnew = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h1);
  hnew = gsl_histogram2d_clone(h1);
  if (HISTOGRAM2D_P(hh2)) {
    Data_Get_Struct(hh2, gsl_histogram2d, h2);
    gsl_histogram2d_add(hnew, h2);
  } else {
    Need_Float(hh2);
    gsl_histogram2d_shift(hnew, NUM2DBL(hh2));
  }
  return Data_Wrap_Struct(CLASS_OF(h1), 0, gsl_histogram2d_free, hnew);
}

static VALUE rb_gsl_histogram2d_sub(VALUE obj, VALUE hh2)
{
  gsl_histogram2d *h1 = NULL, *h2 = NULL, *hnew = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h1);
  hnew = gsl_histogram2d_clone(h1);
  if (HISTOGRAM2D_P(hh2)) {
    Data_Get_Struct(hh2, gsl_histogram2d, h2);
    gsl_histogram2d_sub(hnew, h2);
  } else {
    Need_Float(hh2);
    gsl_histogram2d_shift(hnew, -NUM2DBL(hh2));
  }
  return Data_Wrap_Struct(CLASS_OF(h1), 0, gsl_histogram2d_free, hnew);
}

static VALUE rb_gsl_histogram2d_mul(VALUE obj, VALUE hh2)
{
  gsl_histogram2d *h1 = NULL, *h2 = NULL, *hnew = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h1);
  hnew = gsl_histogram2d_clone(h1);
  if (HISTOGRAM2D_P(hh2)) {
    Data_Get_Struct(hh2, gsl_histogram2d, h2);
    gsl_histogram2d_mul(hnew, h2);
  } else {
    Need_Float(hh2);
    gsl_histogram2d_scale(hnew, NUM2DBL(hh2));
  }
  return Data_Wrap_Struct(CLASS_OF(h1), 0, gsl_histogram2d_free, hnew);
}

static VALUE rb_gsl_histogram2d_div(VALUE obj, VALUE hh2)
{
  gsl_histogram2d *h1 = NULL, *h2 = NULL, *hnew = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h1);
  hnew = gsl_histogram2d_clone(h1);
  if (HISTOGRAM2D_P(hh2)) {
    Data_Get_Struct(hh2, gsl_histogram2d, h2);
    gsl_histogram2d_div(hnew, h2);
  } else {
    Need_Float(hh2);
    gsl_histogram2d_scale(hnew, 1.0/NUM2DBL(hh2));
  }
  return Data_Wrap_Struct(CLASS_OF(h1), 0, gsl_histogram2d_free, hnew);
}

static VALUE rb_gsl_histogram2d_scale2(VALUE obj, VALUE val)
{
  gsl_histogram2d *h1 = NULL, *hnew = NULL;
  Need_Float(val);
  Data_Get_Struct(obj, gsl_histogram2d, h1);
  hnew = gsl_histogram2d_clone(h1);
  gsl_histogram2d_scale(hnew, NUM2DBL(val));
  return Data_Wrap_Struct(CLASS_OF(h1), 0, gsl_histogram2d_free, hnew);
}

static VALUE rb_gsl_histogram2d_shift2(VALUE obj, VALUE val)
{
  gsl_histogram2d *h1 = NULL, *hnew = NULL;
  Need_Float(val);
  Data_Get_Struct(obj, gsl_histogram2d, h1);
  hnew = gsl_histogram2d_clone(h1);
  gsl_histogram2d_shift(hnew, NUM2DBL(val));
  return Data_Wrap_Struct(CLASS_OF(h1), 0, gsl_histogram2d_free, hnew);
}

static VALUE rb_gsl_histogram2d_fwrite(VALUE obj, VALUE io)
{
  gsl_histogram2d *h = NULL;
  FILE *f;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  f = rb_gsl_open_writefile(io, &flag);
  status = gsl_histogram2d_fwrite(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram2d_fread(VALUE obj, VALUE io)
{
  gsl_histogram2d *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  f = rb_gsl_open_readfile(io, &flag);
  status = gsl_histogram2d_fread(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram2d_fprintf(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram2d *h = NULL;
  FILE *fp;
  int status, flag = 0;
  if (argc != 1 && argc != 3) {
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 3)", argc);
  }
  Data_Get_Struct(obj, gsl_histogram2d, h);
  fp = rb_gsl_open_writefile(argv[0], &flag);
  if (argc == 3) {
    Check_Type(argv[1], T_STRING);
    Check_Type(argv[2], T_STRING);
    status = gsl_histogram2d_fprintf(fp, h, STR2CSTR(argv[1]), STR2CSTR(argv[2]));
  } else {
    status = gsl_histogram2d_fprintf(fp, h, "%g", "%g");
  }
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram2d_fscanf(VALUE obj, VALUE io)
{
  gsl_histogram2d *h = NULL;
  FILE *fp;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  fp = rb_gsl_open_readfile(io, &flag);
  status = gsl_histogram2d_fscanf(fp, h);
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram2d_reset(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  gsl_histogram2d_reset(h);
  return obj;
}

#ifdef GSL_0_9_4_LATER
static VALUE rb_gsl_histogram2d_pdf_alloc(VALUE klass, VALUE nx, VALUE ny)
{
  gsl_histogram2d_pdf *h = NULL;
  CHECK_FIXNUM(nx); CHECK_FIXNUM(ny);
  h = gsl_histogram2d_pdf_alloc(FIX2INT(nx), FIX2INT(ny));
  return Data_Wrap_Struct(klass, 0, gsl_histogram2d_pdf_free, h);
}

static VALUE rb_gsl_histogram2d_pdf_init(VALUE obj, VALUE hh)
{
  gsl_histogram2d_pdf *pdf = NULL;
  gsl_histogram2d *h = NULL;
  CHECK_HISTOGRAM2D(hh);
  Data_Get_Struct(obj, gsl_histogram2d_pdf, pdf);
  Data_Get_Struct(hh, gsl_histogram2d, h);
  gsl_histogram2d_pdf_init(pdf, h);
  return obj;
}
#else
static VALUE rb_gsl_histogram2d_pdf_alloc(VALUE klass, VALUE hhh)
{
  gsl_histogram2d_pdf *h = NULL;
  gsl_histogram2d *hh;
  Data_Get_Struct(hhh, gsl_histogram2d, hh);
  h = gsl_histogram2d_pdf_alloc(hh);
  return Data_Wrap_Struct(klass, 0, gsl_histogram2d_pdf_free, h);
}

#endif

static VALUE rb_gsl_histogram2d_pdf_sample(VALUE obj, VALUE r1, VALUE r2)
{
  gsl_histogram2d_pdf *pdf = NULL;
  double x, y;
  Need_Float(r1); Need_Float(r2);
  Data_Get_Struct(obj, gsl_histogram2d_pdf, pdf);
  gsl_histogram2d_pdf_sample(pdf, NUM2DBL(r1), NUM2DBL(r2), &x, &y);
  return rb_ary_new3(2, rb_float_new(x), rb_float_new(y));
}

static VALUE rb_gsl_histogram2d_xrange(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  v = gsl_vector_view_alloc(h->nx);
  v->vector.data = h->xrange;
  v->vector.size = h->nx + 1;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_histogram_range, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_histogram2d_yrange(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  v = gsl_vector_view_alloc(h->ny);
  v->vector.data = h->yrange;
  v->vector.size = h->ny + 1;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_histogram_range, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_histogram2d_bin(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  v = gsl_vector_view_alloc(h->nx*h->ny);
  v->vector.data = h->bin;
  v->vector.size = h->nx*h->ny;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_histogram_bin, 0, gsl_vector_view_free, v);
}

void mygsl_histogram2d_yproject(const gsl_histogram2d *h2, size_t istart,
        size_t iend, gsl_histogram *h)
{
  size_t i, j;
  double sum;
  for (j = 0; j < h2->ny; j++) {
    sum = 0.0;
    for (i = istart; i <= iend; i++) {
      if (i >= h2->nx) break;
      sum += gsl_histogram2d_get(h2, i, j);
    }
    h->bin[j] = sum;
  }
}

gsl_histogram* mygsl_histogram2d_calloc_yproject(const gsl_histogram2d *h2,
             size_t istart, size_t iend)
{
  gsl_histogram *h;
  h = gsl_histogram_calloc_range(h2->ny, h2->yrange);
  mygsl_histogram2d_yproject(h2, istart, iend, h);
  return h;
}

void mygsl_histogram2d_xproject(const gsl_histogram2d *h2, size_t jstart,
        size_t jend, gsl_histogram *h)
{
  size_t i, j;
  double sum;
  for (i = 0; i < h2->nx; i++) {
    sum = 0.0;
    for (j = jstart; j <= jend; j++) {
      if (j >= h2->ny) break;
      sum += gsl_histogram2d_get(h2, i, j);
    }
    h->bin[i] = sum;
  }
}

gsl_histogram* mygsl_histogram2d_calloc_xproject(const gsl_histogram2d *h2,
             size_t jstart, size_t jend)
{
  gsl_histogram *h;
  h = gsl_histogram_calloc_range(h2->nx, h2->xrange);
  mygsl_histogram2d_xproject(h2, jstart, jend, h);
  return h;
}

static VALUE rb_gsl_histogram2d_xproject(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram2d *h2 = NULL;
  gsl_histogram *h = NULL;
  size_t jstart = 0, jend;
  Data_Get_Struct(obj, gsl_histogram2d, h2);
  switch (argc) {
  case 2:
    jstart = (size_t) FIX2INT(argv[0]);
    jend = (size_t) FIX2INT(argv[1]);
    break;
  case 1:
    jstart = (size_t) FIX2INT(argv[0]);
    jend = h2->ny;
    break;
  case 0:
    jend = h2->ny;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0-2)", argc);
    break;
  }
  h = mygsl_histogram2d_calloc_xproject(h2, jstart, jend);
  return Data_Wrap_Struct(cgsl_histogram, 0, gsl_histogram_free, h);
}

static VALUE rb_gsl_histogram2d_yproject(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram2d *h2 = NULL;
  gsl_histogram *h = NULL;
  size_t istart = 0, iend;
  Data_Get_Struct(obj, gsl_histogram2d, h2);
  switch (argc) {
  case 2:
    istart = (size_t) FIX2INT(argv[0]);
    iend = (size_t) FIX2INT(argv[1]);
    break;
  case 1:
    istart = (size_t) FIX2INT(argv[0]);
    iend = h2->ny;
    break;
  case 0:
    iend = h2->ny;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0-2)", argc);
    break;
  }
  h = mygsl_histogram2d_calloc_yproject(h2, istart, iend);
  return Data_Wrap_Struct(cgsl_histogram, 0, gsl_histogram_free, h);
}

static int mygsl_histogram2d_fread2(FILE * stream, gsl_histogram2d * h)
{  
  double xmin, xmax, ymin, ymax;
  int status;
  status = gsl_block_raw_fread(stream, &xmin, 1, 1);
  if (status)    return status;  
  status = gsl_block_raw_fread(stream, &xmax, 1, 1);
  if (status)    return status;  
  status = gsl_block_raw_fread(stream, &ymin, 1, 1);
  if (status)    return status;  
  status = gsl_block_raw_fread(stream, &ymax, 1, 1);
  if (status)    return status;  
  gsl_histogram2d_set_ranges_uniform(h, xmin, xmax, ymin, ymax);
  status = gsl_block_raw_fread (stream, h->bin, h->nx*h->ny, 1);  
  if (status)    return status;  
  return status;
}

static int mygsl_histogram2d_fwrite2(FILE * stream, const gsl_histogram2d * h)
{  
  int status;
  status = gsl_block_raw_fwrite (stream, h->xrange, 1, 1);
  if (status)    return status;  
  status = gsl_block_raw_fwrite (stream, h->xrange+h->nx, 1, 1);
  if (status)    return status;  
  status = gsl_block_raw_fwrite (stream, h->yrange, 1, 1);
  if (status)    return status;  
  status = gsl_block_raw_fwrite (stream, h->yrange+h->ny, 1, 1);
  if (status)    return status;  
  status = gsl_block_raw_fwrite (stream, h->bin, h->nx*h->ny, 1);  
  return status;
}

static VALUE rb_gsl_histogram2d_fwrite2(VALUE obj, VALUE io)
{
  gsl_histogram2d *h = NULL;
  FILE *f;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  f = rb_gsl_open_writefile(io, &flag);
  status = mygsl_histogram2d_fwrite2(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram2d_fread2(VALUE obj, VALUE io)
{
  gsl_histogram2d *h = NULL;
  FILE *f;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  f = rb_gsl_open_readfile(io, &flag);
  status = mygsl_histogram2d_fread2(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static gsl_histogram2d* mygsl_histogram2d_calloc_integrate(const gsl_histogram2d *h,
                 int flag)
{
  gsl_histogram2d *hi;
  size_t i, j, k;
  size_t nx = h->nx, ny = h->ny, n = h->nx*h->ny;
  hi = gsl_histogram2d_calloc(nx, ny);
  gsl_histogram2d_set_ranges(hi, h->xrange, nx+1, h->yrange, ny+1);
  if (flag == -1) {
    hi->bin[n-1] = h->bin[n-1];
    i = nx - 1;
    for (j = ny-2, k = 0;; j--, k++) {
      hi->bin[n-1-k] = gsl_histogram2d_get(hi, i, j+1) + gsl_histogram2d_get(h, i, j);
      if (j == 0) break;
    }
    j = ny - 1;
    for (i = nx-2;; i--) {
      hi->bin[i*ny + j] = gsl_histogram2d_get(hi, i+1, j) + gsl_histogram2d_get(h, i, j);
      if (i == 0) break;
    }
    for (i = nx-2;; i--) {
      for (j = ny-2;; j--) {
  hi->bin[i*ny+j] = ((gsl_histogram2d_get(hi, i+1, j) 
          + gsl_histogram2d_get(hi, i, j+1))
         - gsl_histogram2d_get(hi, i+1, j+1))
    + gsl_histogram2d_get(h, i, j);
  if (j == 0) break;
      }
      if (i == 0) break;
    }
  } else {
    hi->bin[0] = h->bin[0];
    for (j = 1; j < ny; j++) hi->bin[j] = gsl_histogram2d_get(hi, 0, j-1) 
             + gsl_histogram2d_get(h, 0, j);
    for (i = 1; i < nx; i++) hi->bin[i*ny] = gsl_histogram2d_get(hi, i-1, 0)
             + gsl_histogram2d_get(h, i, 0);
    for (i = 1; i < nx; i++) {
      for (j = 1; j < ny; j++) {
  hi->bin[i*ny+j] = ((gsl_histogram2d_get(hi, i-1, j) 
          + gsl_histogram2d_get(hi, i, j-1))
         - gsl_histogram2d_get(hi, i-1, j-1))
    + gsl_histogram2d_get(h, i, j);
      }
    }
  }
  return hi;
}

static VALUE rb_gsl_histogram2d_integrate(int argc, VALUE *argv, VALUE obj)
{
  gsl_histogram2d *h = NULL, *hi = NULL;
  int flag;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  switch (argc) {
  case 0:
    flag = 1;
    break;
  case 1:
    flag = FIX2INT(argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    break;
  }
  hi = mygsl_histogram2d_calloc_integrate(h, flag);
  return Data_Wrap_Struct(cgsl_histogram2d_integ, 0, gsl_histogram2d_free, hi);
}

static VALUE rb_gsl_histogram2d_normalize_bang(VALUE obj)
{
  gsl_histogram2d *h = NULL;
  double scale;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  if (CLASS_OF(obj) == cgsl_histogram2d_integ)
    scale = 1.0/h->bin[h->nx*h->ny-1];
  else
    scale = 1.0/gsl_histogram2d_sum(h);
  gsl_histogram2d_scale(h, scale);
  return obj;
}

static VALUE rb_gsl_histogram2d_normalize(VALUE obj)
{
  gsl_histogram2d *h = NULL, *hnew = NULL;
  Data_Get_Struct(obj, gsl_histogram2d, h);
  hnew = gsl_histogram2d_clone(h);
  return rb_gsl_histogram2d_normalize_bang(Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_histogram2d_free, hnew));
}

void Init_gsl_histogram2d(VALUE module)
{
  VALUE cgsl_histogram2d_pdf;

  cgsl_histogram2d = rb_define_class_under(module, "Histogram2d", cGSL_Object);
  cgsl_histogram2d_view = rb_define_class_under(cgsl_histogram2d, "View",
            cgsl_histogram);

  cgsl_histogram2d_integ = rb_define_class_under(cgsl_histogram2d, "Integral", 
                 cgsl_histogram2d);
#ifdef GSL_0_9_4_LATER
  /*  rb_define_singleton_method(cgsl_histogram2d, "new", rb_gsl_histogram2d_alloc, -1);*/
  rb_define_singleton_method(cgsl_histogram2d, "alloc", rb_gsl_histogram2d_alloc, -1);
  rb_define_singleton_method(cgsl_histogram2d, "new_uniform", 
           rb_gsl_histogram2d_alloc_uniform, -1);
  rb_define_singleton_method(cgsl_histogram2d, "alloc_uniform", 
           rb_gsl_histogram2d_alloc_uniform, -1);
#endif

  rb_define_singleton_method(cgsl_histogram2d, "equal_bins_p", 
           rb_gsl_histogram2d_equal_bins_p, 2);
  rb_define_singleton_method(cgsl_histogram2d, "equal_bins_p?", 
           rb_gsl_histogram2d_equal_bins_p2, 2);

  rb_define_method(cgsl_histogram2d, "set_ranges",  
       rb_gsl_histogram2d_set_ranges, -1);
  rb_define_method(cgsl_histogram2d, "set_ranges_uniform", 
       rb_gsl_histogram2d_set_ranges_uniform, -1);

  rb_define_singleton_method(cgsl_histogram2d, "memcpy",  
           rb_gsl_histogram2d_memcpy, 2);
  rb_define_method(cgsl_histogram2d, "clone",  rb_gsl_histogram2d_clone, 0);
  rb_define_alias(cgsl_histogram2d, "duplicate", "clone");
  rb_define_method(cgsl_histogram2d, "increment",  
       rb_gsl_histogram2d_accumulate, -1);
  rb_define_alias(cgsl_histogram2d, "fill", "increment");
  rb_define_alias(cgsl_histogram2d, "accumulate", "increment");

  rb_define_method(cgsl_histogram2d, "increment2",  
       rb_gsl_histogram2d_accumulate2, -1);
  rb_define_alias(cgsl_histogram2d, "accumulate2", "increment2");
  rb_define_alias(cgsl_histogram2d, "fill2", "increment2");

  rb_define_method(cgsl_histogram2d, "get",  rb_gsl_histogram2d_get, -1);
  rb_define_alias(cgsl_histogram2d, "[]", "get");

  rb_define_method(cgsl_histogram2d, "get_xrange", rb_gsl_histogram2d_get_xrange, 1);
  rb_define_method(cgsl_histogram2d, "get_yrange", rb_gsl_histogram2d_get_yrange, 1);
  rb_define_method(cgsl_histogram2d, "xmax",  rb_gsl_histogram2d_xmax, 0);
  rb_define_method(cgsl_histogram2d, "xmin",  rb_gsl_histogram2d_xmin, 0);
  rb_define_method(cgsl_histogram2d, "ymax",  rb_gsl_histogram2d_ymax, 0);
  rb_define_method(cgsl_histogram2d, "ymin",  rb_gsl_histogram2d_ymin, 0);
  rb_define_method(cgsl_histogram2d, "nx",  rb_gsl_histogram2d_nx, 0);
  rb_define_method(cgsl_histogram2d, "ny",  rb_gsl_histogram2d_ny, 0);

  rb_define_method(cgsl_histogram2d, "find",  rb_gsl_histogram2d_find, 2);

  rb_define_method(cgsl_histogram2d, "max_val",  rb_gsl_histogram2d_max_val, 0);
  rb_define_method(cgsl_histogram2d, "max_bin",  rb_gsl_histogram2d_max_bin, 0);
  rb_define_method(cgsl_histogram2d, "min_val",  rb_gsl_histogram2d_min_val, 0);
  rb_define_method(cgsl_histogram2d, "min_bin",  rb_gsl_histogram2d_min_bin, 0);

#ifdef GSL_1_1_LATER
  rb_define_method(cgsl_histogram2d, "xmean",  rb_gsl_histogram2d_xmean, 0);
  rb_define_method(cgsl_histogram2d, "ymean",  rb_gsl_histogram2d_ymean, 0);
  rb_define_method(cgsl_histogram2d, "xsigma",  rb_gsl_histogram2d_xsigma, 0);
  rb_define_method(cgsl_histogram2d, "ysigma",  rb_gsl_histogram2d_ysigma, 0);
  rb_define_method(cgsl_histogram2d, "cov",  rb_gsl_histogram2d_cov, 0);
  rb_define_method(cgsl_histogram2d, "sum",  rb_gsl_histogram2d_sum, 0);
  rb_define_alias(cgsl_histogram2d, "integral", "sum");
#endif

  rb_define_method(cgsl_histogram2d, "add",  rb_gsl_histogram2d_add, 1);
  rb_define_alias(cgsl_histogram2d, "+", "add");
  rb_define_method(cgsl_histogram2d, "sub",  rb_gsl_histogram2d_sub, 1);
  rb_define_alias(cgsl_histogram2d, "-", "sub");
  rb_define_method(cgsl_histogram2d, "mul",  rb_gsl_histogram2d_mul, 1);
  rb_define_alias(cgsl_histogram2d, "*", "mul");
  rb_define_method(cgsl_histogram2d, "div",  rb_gsl_histogram2d_div, 1);
  rb_define_alias(cgsl_histogram2d, "/", "div");
  rb_define_method(cgsl_histogram2d, "scale",  rb_gsl_histogram2d_scale2, 1);
  rb_define_method(cgsl_histogram2d, "shift",  rb_gsl_histogram2d_shift2, 1);

  rb_define_method(cgsl_histogram2d, "scale!",  rb_gsl_histogram2d_scale, 1);
  rb_define_method(cgsl_histogram2d, "shift!",  rb_gsl_histogram2d_shift, 1);

  rb_define_method(cgsl_histogram2d, "fwrite",  rb_gsl_histogram2d_fwrite, 1);
  rb_define_method(cgsl_histogram2d, "fread",  rb_gsl_histogram2d_fread, 1);
  rb_define_method(cgsl_histogram2d, "fwrite2",  rb_gsl_histogram2d_fwrite2, 1);
  rb_define_method(cgsl_histogram2d, "fread2",  rb_gsl_histogram2d_fread2, 1);
  rb_define_method(cgsl_histogram2d, "fprintf",  rb_gsl_histogram2d_fprintf, -1);
  rb_define_method(cgsl_histogram2d, "fscanf",  rb_gsl_histogram2d_fscanf, 3);

  cgsl_histogram2d_pdf = rb_define_class_under(cgsl_histogram2d, "Pdf", cGSL_Object);
#ifdef GSL_0_9_4_LATER
  /*  rb_define_singleton_method(cgsl_histogram2d_pdf, "new", 
      rb_gsl_histogram2d_pdf_alloc, 2);*/
  rb_define_singleton_method(cgsl_histogram2d_pdf, "alloc", 
           rb_gsl_histogram2d_pdf_alloc, 2);
  rb_define_method(cgsl_histogram2d_pdf, "init",  rb_gsl_histogram2d_pdf_init, 1);
#else
  /*  rb_define_singleton_method(cgsl_histogram2d_pdf, "new", 
      rb_gsl_histogram2d_pdf_alloc, 1);*/
  rb_define_singleton_method(cgsl_histogram2d_pdf, "alloc", 
           rb_gsl_histogram2d_pdf_alloc, 1);
#endif

  rb_define_method(cgsl_histogram2d_pdf, "sample",  rb_gsl_histogram2d_pdf_sample, 2);

  rb_define_method(cgsl_histogram2d, "xrange",  rb_gsl_histogram2d_xrange, 0);
  rb_define_method(cgsl_histogram2d, "yrange",  rb_gsl_histogram2d_yrange, 0);
  rb_define_method(cgsl_histogram2d, "bin",  rb_gsl_histogram2d_bin, 0);

  rb_define_method(cgsl_histogram2d, "reset",  rb_gsl_histogram2d_reset, 0);
  rb_define_method(cgsl_histogram2d, "xproject",  rb_gsl_histogram2d_xproject, -1);
  rb_define_method(cgsl_histogram2d, "yproject",  rb_gsl_histogram2d_yproject, -1);

  rb_define_method(cgsl_histogram2d, "integrate",  rb_gsl_histogram2d_integrate, -1);
  rb_undef_method(cgsl_histogram2d_integ, "integrate");

  rb_define_method(cgsl_histogram2d, "normalize",  rb_gsl_histogram2d_normalize, 0);
  rb_define_method(cgsl_histogram2d, "normalize!",  rb_gsl_histogram2d_normalize_bang, 0);
}

#ifdef HISTOGRAM2D_P
#undef HISTOGRAM2D_P
#endif
#ifdef CHECK_HISTOGRAM2D
#undef CHECK_HISTOGRAM2D
#endif

