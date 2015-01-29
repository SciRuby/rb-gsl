/*
  histogram3d.c
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

static VALUE cgsl_histogram3d;
static VALUE cgsl_histogram3d_view;

static VALUE rb_gsl_histogram3d_new(int argc, VALUE *argv, VALUE klass)
{
  mygsl_histogram3d *h = NULL;
  gsl_vector *xrange, *yrange, *zrange;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  size_t nx, ny, nz;
  switch (argc) {
  case 3:
    if (TYPE(argv[0]) == T_FIXNUM 
  && TYPE(argv[1]) == T_FIXNUM && TYPE(argv[2]) == T_FIXNUM) {
      h = mygsl_histogram3d_alloc(FIX2INT(argv[0]), 
          FIX2INT(argv[1]), FIX2INT(argv[2])); 
    } else if (VECTOR_P(argv[0]) && VECTOR_P(argv[1]) && VECTOR_P(argv[2])) {
      Data_Get_Struct(argv[0], gsl_vector, xrange);
      Data_Get_Struct(argv[1], gsl_vector, yrange);
      Data_Get_Struct(argv[2], gsl_vector, zrange);
      h = mygsl_histogram3d_alloc(xrange->size-1, yrange->size-1, zrange->size-1);
      mygsl_histogram3d_set_ranges(h, xrange->data, xrange->size,
           yrange->data, yrange->size,
           zrange->data, zrange->size);
    } else if (TYPE(argv[0]) == T_ARRAY 
         && TYPE(argv[1]) == T_ARRAY && TYPE(argv[2]) == T_ARRAY) {
      xrange = make_cvector_from_rarray(argv[0]);
      yrange = make_cvector_from_rarray(argv[1]);
      zrange = make_cvector_from_rarray(argv[2]);
      h = mygsl_histogram3d_alloc(xrange->size-1, yrange->size-1, zrange->size-1);
      mygsl_histogram3d_set_ranges(h, xrange->data, xrange->size,
           yrange->data, yrange->size,
           zrange->data, zrange->size);
      gsl_vector_free(zrange);
      gsl_vector_free(yrange);
      gsl_vector_free(xrange);
    } else {
      rb_raise(rb_eTypeError, "wrong argument types");
    }
    break;
  case 6:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[2]); CHECK_FIXNUM(argv[4]);
    Check_Type(argv[1], T_ARRAY); Check_Type(argv[3], T_ARRAY); 
    Check_Type(argv[5], T_ARRAY);
    nx = FIX2INT(argv[0]); ny = FIX2INT(argv[2]); nz = FIX2INT(argv[4]);
    xmin = NUM2DBL(rb_ary_entry(argv[1], 0));
    xmax = NUM2DBL(rb_ary_entry(argv[1], 1));
    ymin = NUM2DBL(rb_ary_entry(argv[3], 0));
    ymax = NUM2DBL(rb_ary_entry(argv[3], 1));
    zmin = NUM2DBL(rb_ary_entry(argv[5], 0));
    zmax = NUM2DBL(rb_ary_entry(argv[5], 1));
    h = mygsl_histogram3d_calloc_uniform(nx, ny, nz, xmin, xmax, ymin, ymax,
           zmin, zmax);
    break;
  case 9:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[3]); CHECK_FIXNUM(argv[6]);
    nx = FIX2INT(argv[0]); ny = FIX2INT(argv[3]); nz = FIX2INT(argv[6]);
    xmin = NUM2DBL(argv[1]);    xmax = NUM2DBL(argv[2]);
    ymin = NUM2DBL(argv[4]);    ymax = NUM2DBL(argv[5]);
    zmin = NUM2DBL(argv[7]);    zmax = NUM2DBL(argv[8]);
    h = mygsl_histogram3d_calloc_uniform(nx, ny, nz, xmin, xmax, ymin, ymax,
           zmin, zmax);
    break;
  default:
    break;
  }
  return Data_Wrap_Struct(cgsl_histogram3d, 0, mygsl_histogram3d_free, h);
}

static VALUE rb_gsl_histogram3d_nx(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return INT2FIX(h->nx);
}

static VALUE rb_gsl_histogram3d_ny(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return INT2FIX(h->ny);
}

static VALUE rb_gsl_histogram3d_nz(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return INT2FIX(h->nz);
}

static VALUE rb_gsl_histogram3d_size(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return INT2NUM(h->nx*h->ny*h->nz);
}

static VALUE rb_gsl_histogram3d_shape(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_ary_new3(3, INT2FIX(h->nx), INT2FIX(h->ny), INT2FIX(h->nz));
}

static VALUE rb_gsl_histogram3d_xrange(VALUE obj)
{
  mygsl_histogram3d *h = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  v = gsl_vector_view_alloc(h->nx);
  v->vector.data = h->xrange;
  v->vector.size = h->nx + 1;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_histogram_range, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_histogram3d_yrange(VALUE obj)
{
  mygsl_histogram3d *h = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  v = gsl_vector_view_alloc(h->ny);
  v->vector.data = h->yrange;
  v->vector.size = h->ny + 1;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_histogram_range, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_histogram3d_zrange(VALUE obj)
{
  mygsl_histogram3d *h = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  v = gsl_vector_view_alloc(h->nz);
  v->vector.data = h->zrange;
  v->vector.size = h->nz + 1;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_histogram_range, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_histogram3d_bin(VALUE obj)
{
  mygsl_histogram3d *h = NULL;
  size_t n = h->nx*h->ny*h->nz;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  v = gsl_vector_view_alloc(n);
  v->vector.data = h->bin;
  v->vector.size = n;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_vector_view, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_histogram3d_get(int argc, VALUE *argv, VALUE obj)
{
  mygsl_histogram3d *h = NULL;
  mygsl_histogram3d_view *h2 = NULL;
  mygsl_histogram2d_view *h1 = NULL;
  size_t i, j, k;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  switch (argc) {
  case 1:
    switch (TYPE(argv[0])) {
    case T_FIXNUM:
      i = FIX2INT(argv[0]);
      h2 = ALLOC(mygsl_histogram3d_view);
      h2->h.nx = h->ny;
      h2->h.ny = h->nz;
      h2->h.xrange = h->yrange;
      h2->h.yrange = h->zrange;
      h2->h.bin = h->bin + i*h->ny*h->nz;
      return Data_Wrap_Struct(cgsl_histogram3d_view, 0, free, h2);
      break;
    case T_ARRAY:
      //      switch (RARRAY(argv[0])->len) {
      switch (RARRAY_LEN(argv[0])) {
      case 1:
  i = FIX2INT(rb_ary_entry(argv[0], 0));
  h2 = ALLOC(mygsl_histogram3d_view);
  h2->h.nx = h->ny;
  h2->h.ny = h->nz;
  h2->h.xrange = h->yrange;
  h2->h.yrange = h->zrange;
  h2->h.bin = h->bin + i*h->ny*h->nz;
  return Data_Wrap_Struct(cgsl_histogram3d_view, 0, free, h2);
  break;
      case 2:
  i = FIX2INT(rb_ary_entry(argv[0], 0));
  j = FIX2INT(rb_ary_entry(argv[0], 1));
  h1 = ALLOC(mygsl_histogram2d_view);
  h1->h.n = h->nz;
  h1->h.range = h->zrange;
  h1->h.bin = h->bin + i*h->ny*h->nz + j*h->nz;
  return Data_Wrap_Struct(cgsl_histogram2d_view, 0, free, h1);    
  break;
      case 3:
  i = FIX2INT(rb_ary_entry(argv[0], 0));
  j = FIX2INT(rb_ary_entry(argv[0], 1));
  k = FIX2INT(rb_ary_entry(argv[0], 2));
  /* do the last line of this function */
  break;
      default:
  rb_raise(rb_eRuntimeError, "wrong array size");
      }
      break;
    default:
      rb_raise(rb_eTypeError, "wrong argument type %s (Fixnum or Array expected)",
         rb_class2name(CLASS_OF(argv[0])));
      break;
    }
    break;
  case 2:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    i = FIX2INT(argv[0]); j = FIX2INT(argv[1]);
    h1 = ALLOC(mygsl_histogram2d_view);
    h1->h.n = h->nz;
    h1->h.range = h->zrange;
    h1->h.bin = h->bin + i*h->ny*h->nz + j*h->nz;
    return Data_Wrap_Struct(cgsl_histogram2d_view, 0, free, h1);    
    break;
  case 3:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
    i = FIX2INT(argv[0]); j = FIX2INT(argv[1]); k = FIX2INT(argv[2]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arugments");
    break;
  }
  return rb_float_new(mygsl_histogram3d_get(h, i, j, k));
}

static VALUE rb_gsl_histogram3d_increment(int argc, VALUE *argv, VALUE obj)
{
  mygsl_histogram3d *h = NULL;
  double x, y, z, weight = 1;
  switch (argc) {
  case 4:
    Need_Float(argv[3]);
    weight = NUM2DBL(argv[3]);
    /* no break */
  case 3:
    Need_Float(argv[0]); Need_Float(argv[1]); Need_Float(argv[2]);
    x = NUM2DBL(argv[0]); y = NUM2DBL(argv[1]); z = NUM2DBL(argv[2]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arugments (%d for 3 or 4", argc);
    break;
  }
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  mygsl_histogram3d_accumulate(h, x, y, z, weight);
  return obj;
}

static VALUE rb_gsl_histogram3d_increment2(int argc, VALUE *argv, VALUE obj)
{
  mygsl_histogram3d *h = NULL;
  double x, y, z, weight = 1;
  switch (argc) {
  case 4:
    Need_Float(argv[3]);
    weight = NUM2DBL(argv[3]);
    /* no break */
  case 3:
    Need_Float(argv[0]); Need_Float(argv[1]); Need_Float(argv[2]);
    x = NUM2DBL(argv[0]); y = NUM2DBL(argv[1]); z = NUM2DBL(argv[2]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arugments (%d for 3 or 4", argc);
    break;
  }
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  mygsl_histogram3d_accumulate2(h, x, y, z, weight);
  return obj;
}

static VALUE rb_gsl_histogram3d_get_xrange(VALUE obj, VALUE ii)
{
  mygsl_histogram3d *h = NULL;
  size_t i;
  double xlower, xupper;
  CHECK_FIXNUM(ii);
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  i = (size_t) FIX2INT(ii);
  mygsl_histogram3d_get_xrange(h, i, &xlower, &xupper);
  return rb_ary_new3(2, rb_float_new(xlower), rb_float_new(xupper));
}

static VALUE rb_gsl_histogram3d_get_yrange(VALUE obj, VALUE jj)
{
  mygsl_histogram3d *h = NULL;
  size_t j;
  double ylower, yupper;
  CHECK_FIXNUM(jj);
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  j = (size_t) FIX2INT(jj);
  mygsl_histogram3d_get_yrange(h, j, &ylower, &yupper);
  return rb_ary_new3(2, rb_float_new(ylower), rb_float_new(yupper));
}

static VALUE rb_gsl_histogram3d_get_zrange(VALUE obj, VALUE kk)
{
  mygsl_histogram3d *h = NULL;
  size_t k;
  double zlower, zupper;
  CHECK_FIXNUM(kk);
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  k = (size_t) FIX2INT(kk);
  mygsl_histogram3d_get_zrange(h, k, &zlower, &zupper);
  return rb_ary_new3(2, rb_float_new(zlower), rb_float_new(zupper));
}

static VALUE rb_gsl_histogram3d_find(VALUE obj, VALUE x, VALUE y, VALUE z)
{
  mygsl_histogram3d *h = NULL;
  size_t i, j, k;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  mygsl_histogram3d_find(h, NUM2DBL(x), NUM2DBL(y), NUM2DBL(z),
       &i, &j, &k);
  return rb_ary_new3(3, INT2FIX(i), INT2FIX(j), INT2FIX(k));
}

static VALUE rb_gsl_histogram3d_set_ranges(VALUE obj, VALUE xx, VALUE yy, VALUE zz)
{
  mygsl_histogram3d *h = NULL;
  gsl_vector *xrange, *yrange, *zrange;
  int flagx = 0, flagy = 0, flagz = 0;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  if (TYPE(xx) == T_ARRAY) {
    xrange = make_cvector_from_rarray(xx);
    flagx = 1;
  } else if (VECTOR_P(xx)) {
    Data_Get_Struct(xx, gsl_vector, xrange);
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (Array or Vector expected)",
       rb_class2name(CLASS_OF(xx)));
  }
  if (xrange->size != h->nx+1)
    rb_raise(rb_eIndexError, "xrange length is different");
  if (TYPE(yy) == T_ARRAY) {
    yrange = make_cvector_from_rarray(yy);
    flagy = 1;
  } else if (VECTOR_P(yy)) {
    Data_Get_Struct(yy, gsl_vector, yrange);
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (Array or Vector expected)",
       rb_class2name(CLASS_OF(yy)));
  }
  if (yrange->size != h->ny+1)
    rb_raise(rb_eIndexError, "yrange length is different");
  if (TYPE(zz) == T_ARRAY) {
    zrange = make_cvector_from_rarray(zz);
    flagz = 1;
  } else if (VECTOR_P(zz)) {
    Data_Get_Struct(zz, gsl_vector, zrange);
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (Array or Vector expected)",
       rb_class2name(CLASS_OF(zz)));
  }
  if (zrange->size != h->nz+1)
    rb_raise(rb_eIndexError, "zrange length is different");
  mygsl_histogram3d_set_ranges(h, xrange->data, xrange->size,
             yrange->data, yrange->size, 
             zrange->data, zrange->size);
  if (flagz == 1) gsl_vector_free(zrange);
  if (flagy == 1) gsl_vector_free(yrange);
  if (flagx == 1) gsl_vector_free(xrange);
  return obj;
}

static VALUE rb_gsl_histogram3d_set_ranges_uniform(int argc, VALUE *argv, VALUE obj)
{
  mygsl_histogram3d *h = NULL;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  switch (argc) {
  case 3:
    Check_Type(argv[0], T_ARRAY); Check_Type(argv[1], T_ARRAY);
    Check_Type(argv[2], T_ARRAY);
    xmin = NUM2DBL(rb_ary_entry(argv[0], 0));
    xmax = NUM2DBL(rb_ary_entry(argv[0], 1));
    ymin = NUM2DBL(rb_ary_entry(argv[1], 0));
    ymax = NUM2DBL(rb_ary_entry(argv[1], 1));
    zmin = NUM2DBL(rb_ary_entry(argv[2], 0));
    zmax = NUM2DBL(rb_ary_entry(argv[2], 1));
    break;
  case 6:
    xmin = NUM2DBL(argv[0]); xmax = NUM2DBL(argv[1]);
    ymin = NUM2DBL(argv[2]); ymax = NUM2DBL(argv[3]);
    zmin = NUM2DBL(argv[4]); zmax = NUM2DBL(argv[5]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 3 or 6)", argc);
    break;
  }
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  mygsl_histogram3d_set_ranges_uniform(h, xmin, xmax, ymin, ymax, zmin, zmax);
  return obj;
}

static VALUE rb_gsl_histogram3d_memcpy(VALUE obj, VALUE a, VALUE b)
{
  mygsl_histogram3d *dst, *src;
  Data_Get_Struct(a, mygsl_histogram3d, dst);
  Data_Get_Struct(b, mygsl_histogram3d, src);
  mygsl_histogram3d_memcpy(dst, src);
  return a;
}

static VALUE rb_gsl_histogram3d_clone(VALUE obj)
{
  mygsl_histogram3d *h, *hnew;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  hnew = mygsl_histogram3d_clone(h);
  return Data_Wrap_Struct(cgsl_histogram3d, 0, mygsl_histogram3d_free, hnew);
}

static VALUE rb_gsl_histogram3d_xyproject(int argc, VALUE *argv, VALUE obj)
{
  mygsl_histogram3d *h3;
  gsl_histogram2d *h2;
  size_t kstart = 0, kend;
  Data_Get_Struct(obj, mygsl_histogram3d, h3);
  switch (argc) {
  case 2:
    kstart = (size_t) FIX2INT(argv[0]);
    kend = (size_t) FIX2INT(argv[1]);
    break;
  case 1:
    kstart = (size_t) FIX2INT(argv[0]);
    kend = h3->nz;
    break;
  case 0:
    kend = h3->nz;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0-2)", argc);
    break;
  }
  h2 = mygsl_histogram3d_xyproject(h3, kstart, kend);
  return Data_Wrap_Struct(cgsl_histogram2d, 0, gsl_histogram2d_free, h2);
}

static VALUE rb_gsl_histogram3d_xzproject(int argc, VALUE *argv, VALUE obj)
{
  mygsl_histogram3d *h3;
  gsl_histogram2d *h2;
  size_t jstart = 0, jend;
  Data_Get_Struct(obj, mygsl_histogram3d, h3);
  switch (argc) {
  case 2:
    jstart = (size_t) FIX2INT(argv[0]);
    jend = (size_t) FIX2INT(argv[1]);
    break;
  case 1:
    jstart = (size_t) FIX2INT(argv[0]);
    jend = h3->ny;
    break;
  case 0:
    jend = h3->ny;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0-2)", argc);
    break;
  }
  h2 = mygsl_histogram3d_xzproject(h3, jstart, jend);
  return Data_Wrap_Struct(cgsl_histogram2d, 0, gsl_histogram2d_free, h2);
}

static VALUE rb_gsl_histogram3d_yzproject(int argc, VALUE *argv, VALUE obj)
{
  mygsl_histogram3d *h3;
  gsl_histogram2d *h2;
  size_t istart = 0, iend;
  Data_Get_Struct(obj, mygsl_histogram3d, h3);
  switch (argc) {
  case 2:
    istart = (size_t) FIX2INT(argv[0]);
    iend = (size_t) FIX2INT(argv[1]);
    break;
  case 1:
    istart = (size_t) FIX2INT(argv[0]);
    iend = h3->nx;
    break;
  case 0:
    iend = h3->nx;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0-2)", argc);
    break;
  }
  h2 = mygsl_histogram3d_yzproject(h3, istart, iend);
  return Data_Wrap_Struct(cgsl_histogram2d, 0, gsl_histogram2d_free, h2);
}

static VALUE rb_gsl_histogram3d_scale_bang(VALUE obj, VALUE s)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  mygsl_histogram3d_scale(h, NUM2DBL(s));
  return obj;
}

static VALUE rb_gsl_histogram3d_scale(VALUE obj, VALUE s)
{
  mygsl_histogram3d *h, *hnew;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  hnew = mygsl_histogram3d_clone(h);
  mygsl_histogram3d_scale(hnew, NUM2DBL(s));
  return Data_Wrap_Struct(cgsl_histogram3d, 0, mygsl_histogram3d_free, hnew);
}

static VALUE rb_gsl_histogram3d_shift_bang(VALUE obj, VALUE s)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  mygsl_histogram3d_shift(h, NUM2DBL(s));
  return obj;
}

static VALUE rb_gsl_histogram3d_shift(VALUE obj, VALUE s)
{
  mygsl_histogram3d *h, *hnew;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  hnew = mygsl_histogram3d_clone(h);
  mygsl_histogram3d_shift(hnew, NUM2DBL(s));
  return Data_Wrap_Struct(cgsl_histogram3d, 0, mygsl_histogram3d_free, hnew);
}

static VALUE rb_gsl_histogram3d_xmax(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_xmax(h));
}

static VALUE rb_gsl_histogram3d_xmin(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_xmin(h));
}

static VALUE rb_gsl_histogram3d_ymax(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_ymax(h));
}

static VALUE rb_gsl_histogram3d_ymin(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_ymin(h));
}

static VALUE rb_gsl_histogram3d_zmax(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_zmax(h));
}

static VALUE rb_gsl_histogram3d_zmin(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_zmin(h));
}

static VALUE rb_gsl_histogram3d_fwrite(VALUE obj, VALUE io)
{
  mygsl_histogram3d *h = NULL;
  FILE *f;
  int status, flag = 0;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  f = rb_gsl_open_writefile(io, &flag);
  status = mygsl_histogram3d_fwrite(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram3d_fread(VALUE obj, VALUE io)
{
  mygsl_histogram3d *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  f = rb_gsl_open_readfile(io, &flag);
  status = mygsl_histogram3d_fread(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_histogram3d_max_val(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_max_val(h));
}

static VALUE rb_gsl_histogram3d_min_val(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_min_val(h));
}

static VALUE rb_gsl_histogram3d_max_bin(VALUE obj)
{
  mygsl_histogram3d *h;
  size_t i, j, k;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  mygsl_histogram3d_max_bin(h, &i, &j, &k);
  return rb_ary_new3(3, INT2FIX(i), INT2FIX(j), INT2FIX(k));
}

static VALUE rb_gsl_histogram3d_min_bin(VALUE obj)
{
  mygsl_histogram3d *h;
  size_t i, j, k;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  mygsl_histogram3d_min_bin(h, &i, &j, &k);
  return rb_ary_new3(3, INT2FIX(i), INT2FIX(j), INT2FIX(k));
}

static VALUE rb_gsl_histogram3d_sum(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_sum(h));
}

static VALUE rb_gsl_histogram3d_xmean(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_xmean(h));
}

static VALUE rb_gsl_histogram3d_ymean(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_ymean(h));
}

static VALUE rb_gsl_histogram3d_zmean(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_zmean(h));
}

static VALUE rb_gsl_histogram3d_xsigma(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_xsigma(h));
}

static VALUE rb_gsl_histogram3d_ysigma(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_ysigma(h));
}

static VALUE rb_gsl_histogram3d_zsigma(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  return rb_float_new(mygsl_histogram3d_zsigma(h));
}

static VALUE rb_gsl_histogram3d_reset(VALUE obj)
{
  mygsl_histogram3d *h;
  Data_Get_Struct(obj, mygsl_histogram3d, h);
  mygsl_histogram3d_reset(h);
  return obj;
}

static VALUE rb_gsl_histogram3d_oper(VALUE obj, VALUE hh,
             int (*func)(mygsl_histogram3d *, const mygsl_histogram3d*))
{
  mygsl_histogram3d *h1, *h2, *hnew;
  CHECK_HISTOGRAM3D(hh);
  Data_Get_Struct(obj, mygsl_histogram3d, h1);
  Data_Get_Struct(hh, mygsl_histogram3d, h2);
  hnew = mygsl_histogram3d_clone(h1);
  (*func)(hnew, h2);
  return Data_Wrap_Struct(cgsl_histogram, 0, mygsl_histogram3d_free, hnew);
}

static VALUE rb_gsl_histogram3d_add(VALUE obj, VALUE hh)
{
  return rb_gsl_histogram3d_oper(obj, hh, mygsl_histogram3d_add);
}

static VALUE rb_gsl_histogram3d_sub(VALUE obj, VALUE hh)
{
  return rb_gsl_histogram3d_oper(obj, hh, mygsl_histogram3d_sub);
}

static VALUE rb_gsl_histogram3d_mul(VALUE obj, VALUE hh)
{
  return rb_gsl_histogram3d_oper(obj, hh, mygsl_histogram3d_mul);
}

static VALUE rb_gsl_histogram3d_div(VALUE obj, VALUE hh)
{
  return rb_gsl_histogram3d_oper(obj, hh, mygsl_histogram3d_div);
}

static VALUE rb_gsl_histogram3d_add_shift(VALUE obj, VALUE hh)
{
  switch (TYPE(hh)) {
  case T_FIXNUM: case T_BIGNUM: case T_FLOAT:
    return rb_gsl_histogram3d_shift(obj, hh);
    break;
  default:
    CHECK_HISTOGRAM3D(hh);
    return rb_gsl_histogram3d_oper(obj, hh, mygsl_histogram3d_add);
    break;
  }
}

static VALUE rb_gsl_histogram3d_sub_shift(VALUE obj, VALUE hh)
{
  switch (TYPE(hh)) {
  case T_FIXNUM: case T_BIGNUM: case T_FLOAT:
    return rb_gsl_histogram3d_shift(obj, rb_float_new(-NUM2DBL(hh)));
    break;
  default:
    CHECK_HISTOGRAM3D(hh);
    return rb_gsl_histogram3d_oper(obj, hh, mygsl_histogram3d_sub);
    break;
  }
}

static VALUE rb_gsl_histogram3d_mul_scale(VALUE obj, VALUE hh)
{
  switch (TYPE(hh)) {
  case T_FIXNUM: case T_BIGNUM: case T_FLOAT:
    return rb_gsl_histogram3d_scale(obj, hh);
    break;
  default:
    CHECK_HISTOGRAM3D(hh);
    return rb_gsl_histogram3d_oper(obj, hh, mygsl_histogram3d_mul);
    break;
  }
}

static VALUE rb_gsl_histogram3d_div_scale(VALUE obj, VALUE hh)
{
  switch (TYPE(hh)) {
  case T_FIXNUM: case T_BIGNUM: case T_FLOAT:
    return rb_gsl_histogram3d_scale(obj, rb_float_new(1.0/NUM2DBL(hh)));
    break;
  default:
    CHECK_HISTOGRAM3D(hh);
    return rb_gsl_histogram3d_oper(obj, hh, mygsl_histogram3d_div);
    break;
  }
}

void Init_gsl_histogram3d(VALUE module)
{
  cgsl_histogram3d = rb_define_class_under(module, "Histogram3d", cGSL_Object);
  cgsl_histogram3d_view = rb_define_class_under(cgsl_histogram3d, "View",
            cgsl_histogram2d);

  /*  rb_define_singleton_method(cgsl_histogram3d, "new", rb_gsl_histogram3d_new, -1);*/
  rb_define_singleton_method(cgsl_histogram3d, "alloc", 
           rb_gsl_histogram3d_new, -1);
  rb_define_singleton_method(cgsl_histogram3d, "memcpy", 
           rb_gsl_histogram3d_memcpy, 2);
  /******/
  rb_define_method(cgsl_histogram3d, "set_ranges", 
       rb_gsl_histogram3d_set_ranges, 3);
  rb_define_method(cgsl_histogram3d, "set_ranges_uniform", 
       rb_gsl_histogram3d_set_ranges_uniform, -1);

  rb_define_method(cgsl_histogram3d, "nx", rb_gsl_histogram3d_nx, 0);
  rb_define_method(cgsl_histogram3d, "ny", rb_gsl_histogram3d_ny, 0);
  rb_define_method(cgsl_histogram3d, "nz", rb_gsl_histogram3d_nz, 0);
  rb_define_method(cgsl_histogram3d, "size", rb_gsl_histogram3d_size, 0);
  rb_define_method(cgsl_histogram3d, "shape", rb_gsl_histogram3d_shape, 0);

  rb_define_method(cgsl_histogram3d, "xrange", rb_gsl_histogram3d_xrange, 0);
  rb_define_method(cgsl_histogram3d, "yrange", rb_gsl_histogram3d_yrange, 0);
  rb_define_method(cgsl_histogram3d, "zrange", rb_gsl_histogram3d_zrange, 0);
  rb_define_method(cgsl_histogram3d, "bin", rb_gsl_histogram3d_bin, 0);

  rb_define_method(cgsl_histogram3d, "get", rb_gsl_histogram3d_get, -1);
  rb_define_alias(cgsl_histogram3d, "[]", "get");

  rb_define_method(cgsl_histogram3d, "increment", 
       rb_gsl_histogram3d_increment, -1);
  rb_define_alias(cgsl_histogram3d, "fill", "increment");
  rb_define_alias(cgsl_histogram3d, "accumulate", "increment");

  rb_define_method(cgsl_histogram3d, "increment2", 
       rb_gsl_histogram3d_increment2, -1);
  rb_define_alias(cgsl_histogram3d, "fill2", "increment2");
  rb_define_alias(cgsl_histogram3d, "accumulate2", "increment2");

  rb_define_method(cgsl_histogram3d, "get_xrange", 
       rb_gsl_histogram3d_get_xrange, 1);
  rb_define_method(cgsl_histogram3d, "get_yrange", 
       rb_gsl_histogram3d_get_yrange, 1);
  rb_define_method(cgsl_histogram3d, "get_zrange", 
       rb_gsl_histogram3d_get_zrange, 1);

  rb_define_method(cgsl_histogram3d, "find", rb_gsl_histogram3d_find, 3);

  rb_define_method(cgsl_histogram3d, "clone", rb_gsl_histogram3d_clone, 0);
  rb_define_alias(cgsl_histogram3d, "duplicate", "clone");

  rb_define_method(cgsl_histogram3d, "reset", rb_gsl_histogram3d_reset, 0);

  /*****/

  rb_define_method(cgsl_histogram3d, "xyproject", 
       rb_gsl_histogram3d_xyproject, -1);
  rb_define_method(cgsl_histogram3d, "xzproject", 
       rb_gsl_histogram3d_xzproject, -1);
  rb_define_method(cgsl_histogram3d, "yzproject", 
       rb_gsl_histogram3d_yzproject, -1);

  rb_define_method(cgsl_histogram3d, "scale", rb_gsl_histogram3d_scale, 1);
  rb_define_method(cgsl_histogram3d, "scale!", rb_gsl_histogram3d_scale_bang, 1);
  rb_define_method(cgsl_histogram3d, "shift", rb_gsl_histogram3d_shift, 1);
  rb_define_method(cgsl_histogram3d, "shift!", rb_gsl_histogram3d_shift_bang, 1);

  rb_define_method(cgsl_histogram3d, "xmax", rb_gsl_histogram3d_xmax, 0);
  rb_define_method(cgsl_histogram3d, "xmin", rb_gsl_histogram3d_xmin, 0);
  rb_define_method(cgsl_histogram3d, "ymax", rb_gsl_histogram3d_ymax, 0);
  rb_define_method(cgsl_histogram3d, "ymin", rb_gsl_histogram3d_ymin, 0);
  rb_define_method(cgsl_histogram3d, "zmax", rb_gsl_histogram3d_zmax, 0);
  rb_define_method(cgsl_histogram3d, "zmin", rb_gsl_histogram3d_zmin, 0);
  /*****/
  rb_define_method(cgsl_histogram3d, "fwrite", rb_gsl_histogram3d_fwrite, 1);
  rb_define_method(cgsl_histogram3d, "fread", rb_gsl_histogram3d_fread, 1);

  rb_define_method(cgsl_histogram3d, "max_val", rb_gsl_histogram3d_max_val, 0);
  rb_define_method(cgsl_histogram3d, "min_val", rb_gsl_histogram3d_min_val, 0);
  rb_define_method(cgsl_histogram3d, "max_bin", rb_gsl_histogram3d_max_bin, 0);
  rb_define_method(cgsl_histogram3d, "min_bin", rb_gsl_histogram3d_min_bin, 0);

  rb_define_method(cgsl_histogram3d, "sum", rb_gsl_histogram3d_sum, 0);
  rb_define_alias(cgsl_histogram3d, "integral", "sum");
  rb_define_method(cgsl_histogram3d, "xmean", rb_gsl_histogram3d_xmean, 0);
  rb_define_method(cgsl_histogram3d, "ymean", rb_gsl_histogram3d_ymean, 0);
  rb_define_method(cgsl_histogram3d, "zmean", rb_gsl_histogram3d_zmean, 0);
  rb_define_method(cgsl_histogram3d, "xsigma", rb_gsl_histogram3d_xsigma, 0);
  rb_define_method(cgsl_histogram3d, "ysigma", rb_gsl_histogram3d_ysigma, 0);
  rb_define_method(cgsl_histogram3d, "zsigma", rb_gsl_histogram3d_zsigma, 0);

  /*****/
  rb_define_method(cgsl_histogram3d, "add", rb_gsl_histogram3d_add, 1);
  rb_define_method(cgsl_histogram3d, "sub", rb_gsl_histogram3d_sub, 1);
  rb_define_method(cgsl_histogram3d, "mul", rb_gsl_histogram3d_mul, 1);
  rb_define_method(cgsl_histogram3d, "div", rb_gsl_histogram3d_div, 1);
  rb_define_method(cgsl_histogram3d, "+", rb_gsl_histogram3d_add_shift, 1);
  rb_define_method(cgsl_histogram3d, "-", rb_gsl_histogram3d_sub_shift, 1);
  rb_define_method(cgsl_histogram3d, "*", rb_gsl_histogram3d_mul_scale, 1);
  rb_define_method(cgsl_histogram3d, "/", rb_gsl_histogram3d_div_scale, 1);

}
