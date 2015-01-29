/*
  vector.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada
        Modified by Seiya Nishizawa        14/Apr/2004

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_array.h"
#include "include/rb_gsl_histogram.h"
#include "include/rb_gsl_complex.h"
#include "include/rb_gsl_poly.h"
#ifdef HAVE_NARRAY_H
#include "include/rb_gsl_with_narray.h"
#endif

VALUE rb_gsl_vector_inner_product(int argc, VALUE *argv, VALUE obj);
static VALUE rb_gsl_vector_product_to_m(int argc, VALUE *argv, VALUE obj);
VALUE rb_gsl_vector_int_to_f(VALUE obj);
VALUE rb_ary_to_gv(VALUE klass, VALUE ary);

static VALUE rb_gsl_vector_Xspace(double min, double max, int i,
          gsl_vector* (*f)(const double, const double, const size_t))
{
  gsl_vector *v = NULL;
  size_t n;
  if (i <= 0) rb_raise(rb_eArgError, "npoints must be greater than 0");
  n = (size_t) i;
  if (n == 1 && min != max) rb_raise(rb_eArgError, "npoints is 1, but x1 != x2");
  v = (*f)(min, max, n);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

static gsl_vector* gsl_vector_linspace(const double min, 
               const double max, const size_t n)
{
  gsl_vector *v = NULL;
  double dx;
  size_t i;
  v = gsl_vector_alloc(n);
  if (n > 1) {
    dx = (max - min)/(n-1);
    gsl_vector_set(v, 0, min);
    for (i = 1; i < n-1; i++) gsl_vector_set(v, i, i*dx + min);
    gsl_vector_set(v, n-1, max);
  } else {
    gsl_vector_set(v, 0, min);
  }
  return v;
}

static gsl_vector* gsl_vector_logspace(const double min, 
               const double max, const size_t n)
{
  gsl_vector *v = NULL;
  double dx;
  size_t i;
  v = gsl_vector_alloc(n);
  if (n > 1) {
    dx = (max - min)/(n-1);
    gsl_vector_set(v, 0, pow(10.0, min));
    for (i = 1; i < n-1; i++) gsl_vector_set(v, i, pow(10.0, dx*i + min));
    gsl_vector_set(v, n-1, pow(10.0, max));
  } else {
    gsl_vector_set(v, 0, pow(10.0, min));
  }
  return v;
}

static VALUE rb_gsl_vector_linspace(int argc, VALUE *argv, VALUE klass)
{
  size_t n = 10;
  switch (argc) {
  case 3:
    CHECK_FIXNUM(argv[2]);
    n = (size_t) FIX2UINT(argv[2]);
    break;
  case 2:
    /* do nothing */
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  Need_Float(argv[0]);   Need_Float(argv[1]); 
  return rb_gsl_vector_Xspace(NUM2DBL(argv[0]), NUM2DBL(argv[1]), n,
            gsl_vector_linspace);
}

static VALUE rb_gsl_vector_logspace(int argc, VALUE *argv, VALUE klass)
{
  size_t n = 10;
  switch (argc) {
  case 3:
    CHECK_FIXNUM(argv[2]);
    n = (size_t) FIX2INT(argv[2]);
    break;
  case 2:
    /* do nothing */
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  Need_Float(argv[0]);   Need_Float(argv[1]); 
  return rb_gsl_vector_Xspace(NUM2DBL(argv[0]), NUM2DBL(argv[1]), n,
            gsl_vector_logspace);
}

static VALUE rb_gsl_vector_logspace2(int argc, VALUE *argv, VALUE klass)
{
  size_t n = 10;
  switch (argc) {
  case 3:
    CHECK_FIXNUM(argv[2]);
    n = (size_t) FIX2INT(argv[2]);
    break;
  case 2:
    /* do nothing */
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  Need_Float(argv[0]);   Need_Float(argv[1]); 
  return rb_gsl_vector_Xspace(log10(NUM2DBL(argv[0])), 
            log10(NUM2DBL(argv[1])), n,
            gsl_vector_logspace);
}

/********************************************************/
VALUE rb_gsl_vector_do_something(VALUE obj, void (*func)(gsl_vector*));

/* singleton */

enum {
  GSL_VECTOR_ADD,
  GSL_VECTOR_SUB,
  GSL_VECTOR_MUL,
  GSL_VECTOR_DIV,
};

static VALUE rb_gsl_vector_arithmetics(int flag, VALUE obj, VALUE bb) 
{
  gsl_vector *v = NULL, *vnew = NULL, *b = NULL;
  gsl_vector_complex *cvnew = NULL, *cb = NULL;
  gsl_complex *c = NULL;
  Data_Get_Struct(obj, gsl_vector, v);
  switch (TYPE(bb)) {
  case T_FLOAT:
  case T_FIXNUM:
    switch (flag) {
    case GSL_VECTOR_ADD:
      vnew = make_vector_clone(v);
      gsl_vector_add_constant(vnew, NUM2DBL(bb));
      break;
    case GSL_VECTOR_SUB:
      vnew = make_vector_clone(v);
      gsl_vector_add_constant(vnew, -NUM2DBL(bb));
      break;
    case GSL_VECTOR_MUL:
      vnew = make_vector_clone(v);
      gsl_vector_scale(vnew, NUM2DBL(bb));
      break;
    case GSL_VECTOR_DIV:
      vnew = make_vector_clone(v);
      gsl_vector_scale(vnew, 1.0/NUM2DBL(bb));
      break;
    }
    if (!VECTOR_VIEW_P(obj)) 
      return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_vector_free, vnew);
    else 
      return Data_Wrap_Struct(VECTOR_ROW_COL(obj), 0, gsl_vector_free, vnew);
    break;
  default:
    if (VECTOR_INT_P(bb)) bb = rb_gsl_vector_int_to_f(bb);
    if (VECTOR_P(bb)) {
      Data_Get_Struct(bb, gsl_vector, b);
      switch (flag) {
      case GSL_VECTOR_ADD:
  vnew = make_vector_clone(v);
  gsl_vector_add(vnew, b);
  break;
      case GSL_VECTOR_SUB:
  vnew = make_vector_clone(v);
  gsl_vector_sub(vnew, b);
  break;
      case GSL_VECTOR_MUL:
  vnew = make_vector_clone(v);
  gsl_vector_mul(vnew, b);
  break;
      case GSL_VECTOR_DIV:
  vnew = make_vector_clone(v);
  gsl_vector_div(vnew, b);
  break;
      }
      if (!VECTOR_VIEW_P(obj)) 
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_vector_free, vnew);
      else 
  return Data_Wrap_Struct(VECTOR_ROW_COL(obj), 0, gsl_vector_free, vnew);
    } else if (VECTOR_COMPLEX_P(bb)) {
      Data_Get_Struct(bb, gsl_vector_complex, cb);
      cvnew = vector_to_complex(v);
      switch (flag) {
      case GSL_VECTOR_ADD:
  gsl_vector_complex_add(cvnew, cb);
  break;
      case GSL_VECTOR_SUB:
  gsl_vector_complex_sub(cvnew, cb);
  break;
      case GSL_VECTOR_MUL:
  gsl_vector_complex_mul(cvnew, cb);
  break;
      case GSL_VECTOR_DIV:
  gsl_vector_complex_div(cvnew, cb);
  break;
      }
      if (VECTOR_COL_P(obj))
  return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, cvnew);
      else
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cvnew);
    } else if (COMPLEX_P(bb)) {
      Data_Get_Struct(bb, gsl_complex, c);
      cvnew = vector_to_complex(v);
      switch (flag) {
      case GSL_VECTOR_ADD:
  gsl_vector_complex_add_constant(cvnew, *c);
  break;
      case GSL_VECTOR_SUB:
  gsl_vector_complex_add_constant(cvnew, gsl_complex_negative(*c));
  break;
      case GSL_VECTOR_MUL:
  gsl_vector_complex_scale(cvnew, *c);
  break;
      case GSL_VECTOR_DIV:
  gsl_vector_complex_scale(cvnew, gsl_complex_inverse(*c));
  break;
      }
      if (VECTOR_COL_P(obj)) 
  return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, cvnew);
      else
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cvnew);
    } else {
      rb_raise(rb_eTypeError, "wrong type argument %s", rb_class2name(CLASS_OF(bb)));
    }
    break;
  }
  /* never reach here */
  return Qnil;
}

VALUE rb_gsl_vector_add(VALUE obj, VALUE b)
{
  return rb_gsl_vector_arithmetics(GSL_VECTOR_ADD, obj, b);
}

VALUE rb_gsl_vector_sub(VALUE obj, VALUE b)
{
  return rb_gsl_vector_arithmetics(GSL_VECTOR_SUB, obj, b);
}

gsl_vector* mygsl_vector_mul_matrix(gsl_vector *v, gsl_matrix *m);
VALUE rb_gsl_vector_mul(VALUE obj, VALUE b)
{
  VALUE argv[2];
  gsl_vector *v, *vnew;
  gsl_matrix *m;
  if (VECTOR_ROW_P(obj) && VECTOR_COL_P(b)) {
    argv[0] = obj;
    argv[1] = b;
    return rb_gsl_vector_inner_product(2, argv, CLASS_OF(obj));
  }
  if (VECTOR_ROW_P(obj) && MATRIX_P(b)) {
    Data_Get_Struct(obj, gsl_vector, v);
    Data_Get_Struct(b, gsl_matrix, m);
    vnew = mygsl_vector_mul_matrix(v, m);
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
  }
  if (VECTOR_COL_P(obj) && VECTOR_ROW_P(b)) {
    argv[0] = obj;
    argv[1] = b;
    return rb_gsl_vector_product_to_m(2, argv, CLASS_OF(obj));
  }
  return rb_gsl_vector_arithmetics(GSL_VECTOR_MUL, obj, b);
}

VALUE rb_gsl_vector_div(VALUE obj, VALUE b)
{
  return rb_gsl_vector_arithmetics(GSL_VECTOR_DIV, obj, b);
}

VALUE rb_ary_to_gv0(VALUE ary)
{
  gsl_vector *v = NULL;
  size_t i, size;
  size = RARRAY_LEN(ary);
  v = gsl_vector_alloc(size);
  if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
  for (i = 0; i < size; i++) {
    gsl_vector_set(v, i, NUM2DBL(rb_ary_entry(ary, i)));
  }
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

VALUE rb_ary_to_gv(VALUE klass, VALUE ary)
{
  gsl_vector *v = NULL;
  size_t i, size;
  size = RARRAY_LEN(ary);
  v = gsl_vector_alloc(size);
  if (v == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
  for (i = 0; i < size; i++) {
    gsl_vector_set(v, i, NUM2DBL(rb_ary_entry(ary, i)));
  }
  return Data_Wrap_Struct(klass, 0, gsl_vector_free, v);
}

VALUE rb_gsl_range_to_gv(VALUE obj)
{
  int beg, en;
  size_t n, i;
  gsl_vector *v = NULL;
  beg = NUM2INT(rb_funcall3(obj, rb_gsl_id_beg, 0, NULL));
  en = NUM2INT(rb_funcall3(obj, rb_gsl_id_end, 0, NULL));
  if (RTEST(rb_funcall3(obj, rb_gsl_id_excl, 0, NULL))) n = en - beg;
  else n = en - beg + 1;
  v = gsl_vector_alloc(n);
  for (i = 0; i < n; i++) gsl_vector_set(v, i, beg + (int)i);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}
gsl_vector_view* gsl_vector_view_alloc()
{
  gsl_vector_view *vv = NULL;
  vv = ALLOC(gsl_vector_view);
  if (vv == NULL) rb_raise(rb_eNoMemError, "malloc failed");
  vv->vector.owner = 0;
  return vv;
}

void gsl_vector_view_free(gsl_vector_view * vv)
{
  free((gsl_vector_view *) vv);
}

static VALUE rb_gsl_vector_to_complex(VALUE obj)
{
  gsl_vector *v = NULL;
  gsl_vector_complex *cv = NULL;
  gsl_complex z;
  size_t i;
  double x;
  Data_Get_Struct(obj, gsl_vector, v);
  cv = gsl_vector_complex_alloc(v->size);
  for (i = 0; i < v->size; i++) {
    x = gsl_vector_get(v, i);
    z.dat[0] = x;
    z.dat[1] = 0;
    gsl_vector_complex_set(cv, i, z);
  }
  if (VECTOR_COL_P(obj))
    return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, cv);
  else
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cv);
}

static VALUE rb_gsl_vector_to_complex2(VALUE obj)
{
  gsl_vector *v = NULL;
  gsl_vector_complex *cv = NULL;
  gsl_complex z;
  size_t i;
  double re, im;
  Data_Get_Struct(obj, gsl_vector, v);
  cv = gsl_vector_complex_alloc(ceil((double)v->size/2));
  for (i = 0; i < v->size; i += 2) {
    re = gsl_vector_get(v, i);
    if (i+1 == v->size) im = 0.0;
    else im = gsl_vector_get(v, i+1);
    z.dat[0] = re;
    z.dat[1] = im;
    gsl_vector_complex_set(cv, i/2, z);
  }
  if (VECTOR_COL_P(obj))
    return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, cv);
  else
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cv);
}

static VALUE rb_gsl_vector_coerce(VALUE obj, VALUE other)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_vector_complex *cv = NULL;
  gsl_complex *c = NULL;
  VALUE vv;
  Data_Get_Struct(obj, gsl_vector, v);
  switch (TYPE(other)) {
  case T_FLOAT:
  case T_FIXNUM:
    vnew = gsl_vector_alloc(v->size);
    if (vnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
    gsl_vector_set_all(vnew, NUM2DBL(other));
    vv = Data_Wrap_Struct(VECTOR_ROW_COL(obj), 0, gsl_vector_free, vnew);
    return rb_ary_new3(2, vv, obj);
    break;
  default:
    if (COMPLEX_P(other)) {
      Data_Get_Struct(other, gsl_complex, c);
      cv = gsl_vector_complex_alloc(v->size);
      if (cv == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
      gsl_vector_complex_set_all(cv, *c);
      if (VECTOR_ROW_P(obj)) 
  vv = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cv);
      else
  vv = Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, cv);
      return rb_ary_new3(2, vv, obj);
    } else if (VECTOR_COMPLEX_P(other)) {
      cv = vector_to_complex(v);
      if (VECTOR_ROW_P(obj)) 
  vv = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cv);
      else
  vv = Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, cv);
      return rb_ary_new3(2, other, vv);
    } else {
      rb_raise(rb_eTypeError, "cannot coerced");
    }
    break;
  }
  return Qnil; /* never reach here */
}

static VALUE rb_gsl_vector_product_to_m(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v = NULL, *v2 = NULL;
  gsl_matrix *m = NULL;
  size_t i, j;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    if (!VECTOR_COL_P(argv[0]))
      rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Vector::Col expected)",
         rb_class2name(CLASS_OF(argv[0])));
    if (!VECTOR_ROW_P(argv[1]))
      rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Vector expected)",
         rb_class2name(CLASS_OF(argv[1])));
    Data_Get_Struct(argv[0], gsl_vector, v);
    Data_Get_Struct(argv[1], gsl_vector, v2);
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    if (!VECTOR_COL_P(obj))
      rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Vector::Col expected)",
         rb_class2name(CLASS_OF(obj)));
    if (!VECTOR_ROW_P(argv[0]))
      rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Vector expected)",
         rb_class2name(CLASS_OF(argv[0])));
    Data_Get_Struct(obj, gsl_vector, v);
    Data_Get_Struct(argv[0], gsl_vector, v2);
    break;
  }
  m = gsl_matrix_alloc(v->size, v2->size);
  for (i = 0; i < v->size; i++) {
    for (j = 0; j < v2->size; j++) {
      gsl_matrix_set(m, i, j, gsl_vector_get(v, i)*gsl_vector_get(v2, j));
    }
  }
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, m);
}

static VALUE rb_gsl_vector_to_f(VALUE obj)
{
  return obj;
}

VALUE rb_gsl_vector_to_i(VALUE obj)
{
  gsl_vector *v = NULL;
  gsl_vector_int *vi = NULL;
  size_t i;
  int val;
  Data_Get_Struct(obj, gsl_vector, v);
  vi = gsl_vector_int_alloc(v->size);
  for (i = 0; i < v->size; i++) {
    val = (int) gsl_vector_get(v, i);
    gsl_vector_int_set(vi, i, val);
  }
  if (VECTOR_COL_P(obj)) 
    return Data_Wrap_Struct(cgsl_vector_int_col, 0, gsl_vector_int_free, vi);
  else
    return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, vi);
}

/* singleton */
#ifdef HAVE_GNU_GRAPH
static void draw_hist(VALUE obj, FILE *fp);
static void draw_vector(VALUE obj, FILE *fp);
static void draw_vector2(VALUE xx, VALUE yy, FILE *fp);
static void draw_vector_array(VALUE ary, FILE *fp);
#ifdef HAVE_NARRAY_H
static void draw_narray(VALUE obj, FILE *fp);
#endif // HAVE_NARRAY_H
#endif // HAVE_GNU_GRAPH

static VALUE rb_gsl_vector_graph2(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_GNU_GRAPH
  size_t i, iend, j, n = 0;
  gsl_vector *x = NULL, *y = NULL;
  gsl_histogram *h = NULL;
  VALUE vx = (VALUE) NULL;
  char command[1024];
  int flag = 0;
  FILE *fp = NULL;
  if (argc < 1)
    rb_raise(rb_eArgError, "two few arguments");

  if (TYPE(argv[argc-1]) == T_STRING) {
    sprintf(command, "graph -T X %s", STR2CSTR(argv[argc-1]));
    iend = argc-1;
  } else {
    strcpy(command, "graph -T X -C -g 3");
    iend = argc;
  }

  if (iend == 1) {
    fp = popen(command, "w");
    if (fp == NULL) rb_raise(rb_eIOError, "GNU graph not found.");
    if (VECTOR_P(argv[0])) {
      draw_vector(argv[0], fp);
    } else if (HISTOGRAM_P(argv[0])) {
      draw_hist(argv[0], fp);
    } else if (TYPE(argv[0]) == T_ARRAY) {
      draw_vector_array(argv[0], fp);
#ifdef HAVE_NARRAY_H
    } else if (NA_IsNArray(argv[0])) {
      draw_narray(argv[0], fp);
#endif
    } else {
      if (fp) pclose(fp);
      rb_raise(rb_eTypeError, "wrong argument type %s", 
         rb_class2name(CLASS_OF(argv[0])));
    }
    if (fp) pclose(fp);
    return Qtrue;
  } else {
    fp = popen(command, "w");
    if (fp == NULL) rb_raise(rb_eIOError, "GNU graph not found.");
    if (VECTOR_P(argv[0])) {
      Data_Get_Struct(argv[0], gsl_vector, x);
      vx = argv[0];
      n = x->size;
      flag = 0;
    } else if (HISTOGRAM_P(argv[0])) {
      Data_Get_Struct(argv[0], gsl_histogram, h);
      x = gsl_vector_alloc(h->n);
      n = x->size;
      for (j = 0; j < x->size; j++)
  gsl_vector_set(x, j, h->range[j]);
      flag = 1;
      draw_hist(argv[0], fp);
      fprintf(fp, "\n");
    } else if (TYPE(argv[0]) == T_ARRAY) {
      draw_vector_array(argv[0], fp);
      fprintf(fp, "\n");
#ifdef HAVE_NARRAY_H
    } else if (NA_IsNArray(argv[0])) {
      vx = argv[0];
      n = NA_TOTAL(argv[0]);
      flag = 0;
#endif
    } else if (NIL_P(argv[0])) {
      if (argc < 2) rb_raise(rb_eArgError, "too few arguments");
      if (VECTOR_P(argv[1])) {
  Data_Get_Struct(argv[1], gsl_vector, y);
  n = y->size;
      } else if (HISTOGRAM_P(argv[1])) {
  Data_Get_Struct(argv[1], gsl_histogram, h);
  n = h->n;
#ifdef HAVE_NARRAY_H
      } else if (NA_IsNArray(argv[1])) {
  n = NA_TOTAL(argv[1]);
#endif
      } else {
  rb_raise(rb_eTypeError, "wrong argument type %s", 
     rb_class2name(CLASS_OF(argv[0])));
      }
      x = gsl_vector_alloc(n);
      for (j = 0; j < n; j++) gsl_vector_set(x, j, (double) j);
      flag = 1;
    } else {
      if (fp) pclose(fp);
      rb_raise(rb_eTypeError, "wrong argument type %s", 
         rb_class2name(CLASS_OF(argv[0])));
    }
    if (flag == 1) vx = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, x);
    for (i = 1; i < iend; i++) {
      if (VECTOR_P(argv[i])) draw_vector2(vx, argv[i], fp);
      else if (HISTOGRAM_P(argv[i])) draw_hist(argv[i], fp);
      else if (TYPE(argv[i]) == T_ARRAY) draw_vector_array(argv[i], fp);
#ifdef HAVE_NARRAY_H
      else if (NA_IsNArray(argv[i])) draw_vector2(vx, argv[i], fp);
#endif
      else 
  rb_raise(rb_eTypeError, "wrong argument type %s", 
     rb_class2name(CLASS_OF(argv[i])));
      fprintf(fp, "\n");
      fflush(fp);
    }
    fclose(fp);
    return Qtrue;
  }
  return Qtrue;
#else
  rb_raise(rb_eNoMethodError, "GNU plotutils required");
  return Qfalse;
#endif
}

#ifdef HAVE_GNU_GRAPH
static void draw_vector(VALUE obj, FILE *fp)
{
  gsl_vector *x = NULL;
  size_t j;
  Data_Get_Vector(obj, x);
  for (j = 0; j < x->size; j++)
    fprintf(fp, "%d %g\n", (int) j, gsl_vector_get(x, j));
  fflush(fp);
}

static void draw_vector2(VALUE xx, VALUE yy, FILE *fp)
{
#ifdef HAVE_NARRAY_H
  struct NARRAY *nax, *nay;
#endif // HAVE_NARRAY_H
  double *ptr1 = NULL, *ptr2 = NULL;
  gsl_vector *vx, *vy;
  size_t j, n, stridex = 1, stridey = 1;
  if (VECTOR_P(xx)) {
    Data_Get_Struct(xx, gsl_vector, vx);
    ptr1 = vx->data;
    n = vx->size;
    stridex = vx->stride;
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(xx)) {
    GetNArray(xx, nax);
    ptr1 = (double *) nax->ptr;
    n = nax->total;
    stridex = 1;
#endif // HAVE_NARRAY_H
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (Vector expected)",
       rb_class2name(CLASS_OF(xx)));
  }
  if (VECTOR_P(yy)) {
    Data_Get_Struct(yy, gsl_vector, vy);
    ptr2 = vy->data;
    n = vy->size;
    stridey = vy->stride;
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(yy)) {
    GetNArray(yy, nay);
    ptr2 = (double *) nay->ptr;
    stridey = 1;
#endif // HAVE_NARRAY_H
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (Vector expected)",
       rb_class2name(CLASS_OF(yy)));
  }
  for (j = 0; j < n; j++)
    fprintf(fp, "%g %g\n", ptr1[j*stridex], ptr2[j*stridey]);
  fflush(fp);
}

#ifdef HAVE_NARRAY_H
static void draw_narray(VALUE obj, FILE *fp)
{
  struct NARRAY *na;
  double *ptr;
  size_t j;
  GetNArray(obj, na);
  ptr = (double *) na->ptr;
  for (j = 0; j < na->total; j++)
    fprintf(fp, "%d %g\n", (int) j, ptr[j]);
  fflush(fp);
}
#endif // HAVE_NARRAY_H

static void draw_hist(VALUE obj, FILE *fp)
{
  gsl_histogram *h = NULL;
  size_t j;
  Data_Get_Struct(obj, gsl_histogram, h);
  for (j = 0; j < h->n; j++) {
    fprintf(fp, "%g %g\n%g %g\n", 
      h->range[j], h->bin[j], h->range[j+1], h->bin[j]);
  }
  fflush(fp);
}

static void draw_vector_array(VALUE ary, FILE *fp)
{
  double *ptrx = NULL, *ptry = NULL, *ptrz = NULL;
  VALUE vx;
  size_t j, n, stridex, stridey, stridez;
  int flag = 0;
  switch (RARRAY_LEN(ary)) {
  case 1:
    flag = 1;
    ptry = get_vector_ptr(rb_ary_entry(ary, 0), &stridey, &n);
    break;
  case 2:
    ptry = get_vector_ptr(rb_ary_entry(ary, 1), &stridey, &n);
    vx = rb_ary_entry(ary, 0);
    if (NIL_P(vx)) {flag = 1;}
    else {
      ptrx = get_vector_ptr(vx, &stridex, &n);
    }
    break;
  case 3:
    ptrz = get_vector_ptr(rb_ary_entry(ary, 2), &stridez, &n);
    ptry = get_vector_ptr(rb_ary_entry(ary, 1), &stridey, &n);
    vx = rb_ary_entry(ary, 0);
    if (NIL_P(vx)) {flag = 2;}
    else {
      ptrx = get_vector_ptr(vx, &stridex, &n);
      flag = 3;
    }
    break;
  default:
    rb_raise(rb_eRuntimeError, "wrong array length (%d for 1 or 2)", 
       (int) RARRAY_LEN(ary));
    break;
  }
  switch (flag) {
  case 0:
    for (j = 0; j < n; j++) 
      fprintf(fp, "%g %g\n", ptrx[j*stridex], ptry[j*stridey]);
    break;
  case 1:
    for (j = 0; j < n; j++) 
      fprintf(fp, "%d %g\n", (int) j, ptry[j*stridey]);
    break;
  case 2:
    for (j = 0; j < n; j++) 
      fprintf(fp, "%d %g %g\n", (int) j, ptry[j*stridey], ptrz[j*stridez]);
    break;
  case 3:
    for (j = 0; j < n; j++) 
      fprintf(fp, "%g %g %g\n", ptrx[j*stridex], ptry[j*stridey], ptrz[j*stridez]);
    break;
  default:
    break;
  }
  fflush(fp);
}
#endif // HAVE_GNU_GRAPH

/* singleton */
static VALUE rb_gsl_vector_plot2(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *x = NULL, *y = NULL, *xerr = NULL, *yerr = NULL;
  FILE *fp = NULL;
  size_t i, n;
  char command[1024];
  fp = popen("gnuplot -persist", "w");
  if (fp == NULL) rb_raise(rb_eIOError, "GNU graph not found.");
  strcpy(command, "plot '-'");
  switch (argc) {
  case 5:
    if (TYPE(argv[4]) == T_STRING)
      sprintf(command, "%s %s", command, STR2CSTR(argv[4]));
    /* no break */
  case 4:
    if (TYPE(argv[3]) == T_STRING) {
      sprintf(command, "%s %s", command, STR2CSTR(argv[3]));
    } else if (VECTOR_P(argv[3])) {
      Data_Get_Struct(argv[3], gsl_vector, yerr);
    } else {
      rb_raise(rb_eTypeError, "argv[3] wrong type %s (String or Vector expected)",
         rb_class2name(CLASS_OF(argv[3])));
    }
    /* no break */
  case 3:
    if (TYPE(argv[2]) == T_STRING) {
      sprintf(command, "%s %s", command, STR2CSTR(argv[2]));
    } else if (VECTOR_P(argv[2])) {
      Data_Get_Struct(argv[2], gsl_vector, xerr);
    } else {
      rb_raise(rb_eTypeError, "argv[2] wrong type %s (String or Vector expected)",
         rb_class2name(CLASS_OF(argv[2])));
    }
    /* no break */
  case 2:
    if (TYPE(argv[1]) == T_STRING) {
      sprintf(command, "%s %s", command, STR2CSTR(argv[1]));
    } else if (VECTOR_P(argv[1])) {
      Data_Get_Struct(argv[1], gsl_vector, y);
    } else {
      rb_raise(rb_eTypeError, "argv[1] wrong type %s (String or Vector expected)",
         rb_class2name(CLASS_OF(argv[1])));
    }
    /* no break */
  case 1:
    if (TYPE(argv[0]) == T_STRING) {
      sprintf(command, "%s %s", command, STR2CSTR(argv[0]));
    } else if (VECTOR_P(argv[0])) {
      Data_Get_Struct(argv[0], gsl_vector, x);
    } else {
      rb_raise(rb_eTypeError, "argv[0] wrong type %s (String or Vector expected)",
         rb_class2name(CLASS_OF(argv[0])));
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of argumens (%d for 1 - 5)", argc);
    break;
  }
  if (x == NULL) rb_raise(rb_eRuntimeError, "x data is not given");
  n = x->size;
  fprintf(fp, "%s\n", command);
  for (i = 0; i < n; i++) {
    if (y == NULL) 
      fprintf(fp, "%d %g\n", (int) i, gsl_vector_get(x, i));
    else if (yerr == NULL)
      fprintf(fp, "%g %g\n", gsl_vector_get(x, i), gsl_vector_get(y, i));
    else if (xerr) 
      fprintf(fp, "%g %g %g %g\n", gsl_vector_get(x, i), gsl_vector_get(y, i),
        gsl_vector_get(xerr, i), gsl_vector_get(yerr, i));
    else
     fprintf(fp, "%g %g %g\n", gsl_vector_get(x, i), gsl_vector_get(y, i),
       gsl_vector_get(yerr, i));
  }
  fprintf(fp, "e\n");
  fflush(fp);
  pclose(fp);
  fp = NULL;
  return Qtrue;
}

static VALUE rb_gsl_vector_normalize(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v = NULL, *vnew = NULL;
  //  double mean;
  double nrm;
  switch (argc) {
  case 0:
    nrm = 1.0;
    break;
  case 1:
    Need_Float(argv[0]);
    nrm = NUM2DBL(argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  Data_Get_Vector(obj, v);
  vnew = make_vector_clone(v);
  /*  mean = gsl_stats_mean(v->data, v->stride, v->size);
  gsl_vector_add_constant(vnew, -mean);
  sd = gsl_stats_sd(vnew->data, vnew->stride, vnew->size);  
  gsl_vector_scale(vnew, sqrt(nrm)/sd);*/
  gsl_vector_scale(vnew, nrm/gsl_blas_dnrm2(v));
  return Data_Wrap_Struct(VECTOR_ROW_COL(obj), 0, gsl_vector_free, vnew);  
}

static VALUE rb_gsl_vector_normalize_bang(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v = NULL;
  //  double mean;
  double nrm;
  double factor;
  switch (argc) {
  case 0:
    nrm = 1.0;
    break;
  case 1:
    Need_Float(argv[0]);
    nrm = NUM2DBL(argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  Data_Get_Vector(obj, v);
  /*  mean = gsl_stats_mean(v->data, v->stride, v->size);
  gsl_vector_add_constant(v, -mean);
  sd = gsl_stats_sd(v->data, v->stride, v->size);  
  gsl_vector_scale(v, sqrt(nrm)/sd);*/
  factor = nrm/gsl_blas_dnrm2(v);
  gsl_vector_scale(v, factor);
  return obj;
}

#ifdef HAVE_NARRAY_H
static VALUE rb_gsl_vector_filescan_na(VALUE klass, VALUE file)
{
  FILE *fp = NULL;
  int nn;
  char buf[1024], filename[1024], *p;
  size_t n, lines, i, j;
  double **ptr;
  double val;
  int shape[1];
  VALUE ary, na;
  Check_Type(file, T_STRING);
  strcpy(filename, STR2CSTR(file));
  sprintf(buf, "wc %s", filename);
  fp = popen(buf, "r");
  if (fgets(buf, 1024, fp) == NULL)
    rb_sys_fail(0);
  pclose(fp);
  sscanf(buf, "%d", &nn);
  lines = (size_t) nn;      /* vector length */
  shape[0] = lines;
  fp = fopen(filename, "r");
  if (fgets(buf, 1024, fp) == NULL)
    rb_sys_fail(0);
  n = count_columns(buf);   /* number of vectors created */
  ptr = (double**) xmalloc(sizeof(double**)*n);
  ary = rb_ary_new2(n);
  p = buf;
  for (j = 0; j < n; j++) {
    na = na_make_object(NA_DFLOAT, 1, shape, cNArray);
    rb_ary_store(ary, j, na);
    ptr[j] = NA_PTR_TYPE(na, double*);
    p = str_scan_double(p, &val);
    if (p) ptr[j][0] = val;
    else break;
  }
  for (i = 1; i < lines; i++) {
    if (fgets(buf, 1024, fp) == NULL)
      rb_sys_fail(0);
    p = buf;
    for (j = 0; j < n; j++) {
      p = str_scan_double(p, &val);
      if (p) ptr[j][i] = val;
      else break;      
    }
  }
  fclose(fp);
  free(ptr);
  return ary;
}
#endif

static VALUE rb_gsl_vector_decimate(VALUE obj, VALUE nn)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_vector_view vv;
  size_t i, n, n2, n3;
  CHECK_FIXNUM(nn);
  Data_Get_Vector(obj, v);
  n = (size_t) FIX2INT(nn);
  if (n > v->size) 
    rb_raise(rb_eArgError, 
       "decimation factor must be smaller than the vector length.");
  if (n == 0) rb_raise(rb_eArgError, "decimation factor must be greater than 1");
  n2 = (size_t) ceil((double)v->size/n);
  vnew = gsl_vector_alloc(n2);
  n3 = n - (n2*n - v->size);
  for (i = 0; i < n2; i++) {
    if (i == n2-1) vv = gsl_vector_subvector(v, i*n, n3);
    else vv = gsl_vector_subvector(v, i*n, n);
    gsl_vector_set(vnew, i, gsl_stats_mean(vv.vector.data, vv.vector.stride,
             vv.vector.size));
  }
  return Data_Wrap_Struct(VECTOR_ROW_COL(obj), 0, gsl_vector_free, vnew);
}

static VALUE rb_gsl_vector_xxx(VALUE obj, double (*f)(double))
{
  gsl_vector *v;
  gsl_vector_int *vnew;
  size_t i;
  Data_Get_Struct(obj, gsl_vector, v);
  vnew = gsl_vector_int_alloc(v->size);
  for (i = 0; i < v->size; i++) {
    gsl_vector_int_set(vnew, i, (int) (*f)(gsl_vector_get(v, i)));
  }
  return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, vnew);
}

static VALUE rb_gsl_vector_floor(VALUE obj)
{
  return rb_gsl_vector_xxx(obj, floor);
}

static VALUE rb_gsl_vector_ceil(VALUE obj)
{
  return rb_gsl_vector_xxx(obj, ceil);
}

#ifdef HAVE_ROUND
double round(double x);
#define rb_gsl_round_native round
#else
static double rb_gsl_round_native(double x)
{
  if(!gsl_finite(x)) {
    return x; // nan, +inf, -inf
  } else if(x > 0.0) {
    return floor(x+0.5);
  } else if(x < 0.0) {
    return ceil(x-0.5);
  }
  return x; // +0 or -0
}
#endif

static VALUE rb_gsl_vector_round(VALUE obj)
{
  return rb_gsl_vector_xxx(obj, rb_gsl_round_native);
}

static VALUE rb_gsl_vector_dB(VALUE obj)
{
  gsl_vector *v, *vnew;
  double x;
  size_t i;
  Data_Get_Struct(obj, gsl_vector, v);
  vnew = gsl_vector_alloc(v->size);
  for (i = 0; i < v->size; i++) {
    x = gsl_vector_get(v, i);
    if (x <= 0.0) rb_raise(rb_eRuntimeError, "negative value found.\n");
    gsl_vector_set(vnew, i, 20*log(x));
  }
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
}

static VALUE rb_gsl_vector_sin(VALUE obj)
{
  return rb_gsl_sf_eval1(sin, obj);
}

static VALUE rb_gsl_vector_cos(VALUE obj)
{
  return rb_gsl_sf_eval1(cos, obj);
}

static VALUE rb_gsl_vector_tan(VALUE obj)
{
  return rb_gsl_sf_eval1(tan, obj);
}

static VALUE rb_gsl_vector_exp(VALUE obj)
{
  return rb_gsl_sf_eval1(exp, obj);
}

static VALUE rb_gsl_vector_log(VALUE obj)
{
  return rb_gsl_sf_eval1(log, obj);
}

static VALUE rb_gsl_vector_log10(VALUE obj)
{
  return rb_gsl_sf_eval1(log10, obj);
}

static VALUE rb_gsl_vector_rotate_bang(int argc, VALUE *argv, VALUE klass)
{
  double rad;
  double x, y, c, s;
  gsl_vector *vx, *vy;
  VALUE v0, v1, retval;
  size_t i, n;
  switch (argc) {
  case 2:
    if (TYPE(argv[0]) == T_ARRAY) {
      v0 = rb_ary_entry(argv[0], 0);
      v1 = rb_ary_entry(argv[0], 1);
      if (VECTOR_P(v0) && VECTOR_P(v1)) {
  Data_Get_Struct(v0, gsl_vector, vx);
  Data_Get_Struct(v1, gsl_vector, vy);
  n = (size_t) GSL_MIN(vx->size, vy->size);
  rad = NUM2DBL(argv[1]);
  retval = argv[0];
      } else {
  x = NUM2DBL(rb_ary_entry(argv[0], 0));
  y = NUM2DBL(rb_ary_entry(argv[0], 1));
  rad = NUM2DBL(argv[1]);
  c = cos(rad); s = sin(rad);
  return rb_ary_new3(2, rb_float_new(c*x - s*y), rb_float_new(s*x + c*y));
      }
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (Array expected)",
         rb_class2name(CLASS_OF(argv[0])));
    }
    break;
  case 3:
    if (VECTOR_P(argv[0]) && VECTOR_P(argv[1])) {
      Data_Get_Struct(argv[0], gsl_vector, vx);
      Data_Get_Struct(argv[1], gsl_vector, vy);
      n = (size_t) GSL_MIN(vx->size, vy->size);
      rad = NUM2DBL(argv[1]);
      retval = rb_ary_new3(2, argv[0], argv[1]);
    } else {
      x = NUM2DBL(argv[0]);
      y = NUM2DBL(argv[1]);
      rad = NUM2DBL(argv[2]);
      c = cos(rad); s = sin(rad);
      return rb_ary_new3(2, rb_float_new(c*x - s*y), rb_float_new(s*x + c*y));
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  c = cos(rad); s = sin(rad);
  for (i = 0; i < n; i++) {
    x = gsl_vector_get(vx, i);
    y = gsl_vector_get(vy, i);
    gsl_vector_set(vx, i, c*x - s*y);
    gsl_vector_set(vy, i, s*x + c*y);
  }
  return retval;
}

static VALUE rb_gsl_vector_rotate(int argc, VALUE *argv, VALUE klass)
{
  double rad;
  double x, y, c, s;
  gsl_vector *vx, *vy, *vxnew, *vynew;
  VALUE v0, v1;
  size_t i, n;
  switch (argc) {
  case 2:
    if (TYPE(argv[0]) == T_ARRAY) {
      v0 = rb_ary_entry(argv[0], 0);
      v1 = rb_ary_entry(argv[0], 1);
      if (VECTOR_P(v0) && VECTOR_P(v1)) {
  Data_Get_Struct(v0, gsl_vector, vx);
  Data_Get_Struct(v1, gsl_vector, vy);
  rad = NUM2DBL(argv[1]);
      } else {
  x = NUM2DBL(rb_ary_entry(argv[0], 0));
  y = NUM2DBL(rb_ary_entry(argv[0], 1));
  rad = NUM2DBL(argv[1]);
  c = cos(rad); s = sin(rad);
  return rb_ary_new3(2, rb_float_new(c*x - s*y), rb_float_new(s*x + c*y));
      }
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (Array expected)",
         rb_class2name(CLASS_OF(argv[0])));
    }
    break;
  case 3:
    if (VECTOR_P(argv[0]) && VECTOR_P(argv[1])) {
      Data_Get_Struct(argv[0], gsl_vector, vx);
      Data_Get_Struct(argv[1], gsl_vector, vy);
      rad = NUM2DBL(argv[1]);
    } else {
      x = NUM2DBL(argv[0]);
      y = NUM2DBL(argv[1]);
      rad = NUM2DBL(argv[2]);
      c = cos(rad); s = sin(rad);
      return rb_ary_new3(2, rb_float_new(c*x - s*y), rb_float_new(s*x + c*y));
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  n = (size_t) GSL_MIN(vx->size, vy->size);
  vxnew = gsl_vector_alloc(n);
  vynew = gsl_vector_alloc(n);
  c = cos(rad); s = sin(rad);
  for (i = 0; i < n; i++) {
    x = gsl_vector_get(vx, i);
    y = gsl_vector_get(vy, i);
    gsl_vector_set(vxnew, i, c*x - s*y);
    gsl_vector_set(vynew, i, s*x + c*y);
  }
  return rb_ary_new3(2, Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vxnew),
         Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vynew));
}

#include "gsl/gsl_fit.h"

static VALUE rb_gsl_vector_linearfit(int argc, VALUE *argv, VALUE klass)
{
  gsl_vector *x, *y, *w = NULL;
  double c0, c1, c00, c01, c11, sumsq;
  switch (argc) {
  case 3:
    CHECK_VECTOR(argv[0]); CHECK_VECTOR(argv[1]); CHECK_VECTOR(argv[2]);
    Data_Get_Struct(argv[0], gsl_vector, x);
    Data_Get_Struct(argv[1], gsl_vector, w);
    Data_Get_Struct(argv[2], gsl_vector, y);
    gsl_fit_wlinear(x->data, x->stride, w->data, w->stride,
        y->data, y->stride, y->size, &c0, &c1, &c00, &c01, &c11,
        &sumsq);
    break;
  case 2:
    CHECK_VECTOR(argv[0]); CHECK_VECTOR(argv[1]); 
    Data_Get_Struct(argv[0], gsl_vector, x);
    Data_Get_Struct(argv[1], gsl_vector, y);
    gsl_fit_linear(x->data, x->stride, y->data,y->stride, y->size, 
       &c0, &c1, &c00, &c01, &c11, &sumsq);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 2 or 3).\n", argc);
  }
  return rb_ary_new3(6, rb_float_new(c0), rb_float_new(c1), rb_float_new(c00),
        rb_float_new(c01), rb_float_new(c11), rb_float_new(sumsq));
}

static VALUE rb_gsl_vector_center(VALUE obj)
{
  gsl_vector *v, *vnew;
  double mean;
  Data_Get_Struct(obj, gsl_vector, v);
  mean = gsl_stats_mean(v->data, v->stride, v->size);
  vnew = gsl_vector_alloc(v->size);
  gsl_vector_memcpy(vnew, v);
  gsl_vector_add_constant(vnew, -mean);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
}

static VALUE rb_gsl_vector_clip(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v, *vnew;  
  double hi = 1.0, lo = 0.0;
  double x;
  size_t i;
  Data_Get_Struct(obj, gsl_vector, v);
  switch (argc) {
  case 0:
    break;
  case 1:
    if (TYPE(argv[0]) == T_ARRAY) {
      lo = NUM2DBL(rb_ary_entry(argv[0], 0));
      hi = NUM2DBL(rb_ary_entry(argv[0], 1));
    } else {
      hi = NUM2DBL(argv[0]);
    }
    break;
  case 2:
    lo = NUM2DBL(argv[0]);
    hi = NUM2DBL(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 0-2).\n", argc);
  }
  vnew = gsl_vector_alloc(v->size);
  for (i = 0; i < v->size; i++) {
    x = gsl_vector_get(v, i);
    if (x > hi) x = hi;
    else if (x < lo) x = lo;
    else {};
    gsl_vector_set(vnew, i, x);
  }
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
}

static VALUE rb_gsl_vector_amp_phase(VALUE obj)
{
  gsl_vector *v;
  gsl_vector *amp, *phase;
  double re, im;
  VALUE vamp, vphase;
  size_t i;
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


static VALUE rb_gsl_vector_clean(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v = NULL, *vnew = NULL;
  double eps = 1e-10;
  size_t n, i;
  switch (argc) {
  case 0:
    /* do nothing */
    break;
  case 1:
    Need_Float(argv[0]);
    eps = NUM2DBL(argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments (%d for 0 or 1)", argc);
    break;
  }
  Data_Get_Struct(obj, gsl_vector, v);
  vnew = make_vector_clone(v);
  n = v->size;
  for (i = 0; i < n; i++) if (fabs(vnew->data[i]) < eps) vnew->data[i] = 0.0;
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);;
}

static VALUE rb_gsl_vector_clean_bang(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v = NULL;
  double eps = 1e-10;
  size_t n, i;
  switch (argc) {
  case 0:
    /* do nothing */
    break;
  case 1:
    Need_Float(argv[0]);
    eps = NUM2DBL(argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments (%d for 0 or 1)", argc);
    break;
  }
  Data_Get_Struct(obj, gsl_vector, v);
  n = v->size;
  for (i = 0; i < n; i++) if (fabs(v->data[i]) < eps) v->data[i] = 0.0;
  return obj;
}

VALUE rb_gsl_pow(VALUE obj, VALUE xx, VALUE nn);
VALUE rb_gsl_vector_pow(VALUE obj, VALUE p)
{
  return rb_gsl_pow(Qnil, obj, p);
}

VALUE rb_gsl_vector_pow_bang(VALUE obj, VALUE pp)
{
  gsl_vector *v;
  double p;
  size_t i;
  Data_Get_Struct(obj, gsl_vector, v);
  p = NUM2DBL(pp);
  for (i = 0; i < v->size; i++) {
    gsl_vector_set(v, i, pow(gsl_vector_get(v, i), p));
  }
  return obj;
}

void Init_gsl_vector(VALUE module)
{
  rb_define_singleton_method(cgsl_vector, "linspace", rb_gsl_vector_linspace, -1);
  rb_define_module_function(module, "linspace", rb_gsl_vector_linspace, -1);
  rb_define_singleton_method(cgsl_vector, "logspace", rb_gsl_vector_logspace, -1);
  rb_define_module_function(module, "logspace", rb_gsl_vector_logspace, -1);
  rb_define_singleton_method(cgsl_vector, "logspace2", rb_gsl_vector_logspace2, -1);
  rb_define_module_function(module, "logspace2", rb_gsl_vector_logspace2, -1);

  rb_define_method(cgsl_vector, "add", rb_gsl_vector_add, 1);
  rb_define_alias(cgsl_vector, "+", "add");
  rb_define_method(cgsl_vector, "sub", rb_gsl_vector_sub, 1);
  rb_define_alias(cgsl_vector, "-", "sub");
  rb_define_method(cgsl_vector, "mul", rb_gsl_vector_mul, 1);
  rb_define_alias(cgsl_vector, "*", "mul");
  rb_define_method(cgsl_vector, "div", rb_gsl_vector_div, 1);
  rb_define_alias(cgsl_vector, "/", "div");

  rb_define_method(cgsl_vector, "to_complex", rb_gsl_vector_to_complex, 0);
  rb_define_method(cgsl_vector, "to_complex2", rb_gsl_vector_to_complex2, 0);

  /*****/

  rb_define_method(cgsl_vector, "coerce", rb_gsl_vector_coerce, 1);

  /*****/
  rb_define_method(rb_cArray, "to_gv", rb_ary_to_gv0, 0);
  rb_define_alias(rb_cArray, "to_gslv", "to_gv");
  rb_define_alias(rb_cArray, "to_gsl_vector", "to_gv");
  rb_define_method(rb_cRange, "to_gv", rb_gsl_range_to_gv, 0);
  rb_define_alias(rb_cRange, "to_gslv", "to_gv");
  rb_define_alias(rb_cRange, "to_gsl_vector", "to_gv");
  rb_define_singleton_method(cgsl_vector, "ary_to_gv", rb_ary_to_gv, 1);

  /*****/
  rb_define_method(cgsl_vector, "to_i", rb_gsl_vector_to_i, 0);
  rb_define_method(cgsl_vector, "to_f", rb_gsl_vector_to_f, 0);

  rb_define_singleton_method(cgsl_vector, "graph", rb_gsl_vector_graph2, -1);
  rb_define_module_function(module, "graph", rb_gsl_vector_graph2, -1);
  rb_define_singleton_method(cgsl_vector, "plot", rb_gsl_vector_plot2, -1);

  /*****/

  rb_define_method(cgsl_vector, "normalize", rb_gsl_vector_normalize, -1);
  rb_define_method(cgsl_vector, "normalize!", rb_gsl_vector_normalize_bang, -1);

#ifdef HAVE_NARRAY_H
  rb_define_singleton_method(cgsl_vector, "filescan_na", rb_gsl_vector_filescan_na, 1);
  rb_define_singleton_method(cNArray, "filescan", rb_gsl_vector_filescan_na, 1);
#endif

  rb_define_method(cgsl_vector, "decimate", rb_gsl_vector_decimate, 1);

  rb_define_method(cgsl_vector, "floor", rb_gsl_vector_floor, 0);
  rb_define_method(cgsl_vector, "ceil", rb_gsl_vector_ceil, 0);
  rb_define_method(cgsl_vector, "round", rb_gsl_vector_round, 0);

  rb_define_method(cgsl_vector, "dB", rb_gsl_vector_dB, 0);

  rb_define_method(cgsl_vector, "sin", rb_gsl_vector_sin, 0);
  rb_define_method(cgsl_vector, "cos", rb_gsl_vector_cos, 0);
  rb_define_method(cgsl_vector, "tan", rb_gsl_vector_tan, 0);
  rb_define_method(cgsl_vector, "exp", rb_gsl_vector_exp, 0);
  rb_define_method(cgsl_vector, "log", rb_gsl_vector_log, 0);
  rb_define_method(cgsl_vector, "log10", rb_gsl_vector_log10, 0);

  rb_define_singleton_method(cgsl_vector, "rotate", rb_gsl_vector_rotate, -1);
  rb_define_singleton_method(cgsl_vector, "rotate!", rb_gsl_vector_rotate_bang, -1);

  rb_define_singleton_method(cgsl_vector, "linearfit", rb_gsl_vector_linearfit, -1);

  rb_define_method(cgsl_vector, "center", rb_gsl_vector_center, 0);
  rb_define_method(cgsl_vector, "clip", rb_gsl_vector_clip, -1);

  rb_define_method(cgsl_vector, "amp_phase", rb_gsl_vector_amp_phase, 0);
  rb_define_method(cgsl_vector, "clean", rb_gsl_vector_clean, -1);
  rb_define_method(cgsl_vector, "clean!", rb_gsl_vector_clean_bang, -1);
  
  rb_define_method(cgsl_vector, "pow", rb_gsl_vector_pow, 1);
  rb_define_alias(cgsl_vector, "**", "pow");
  
  rb_define_method(cgsl_vector, "pow!", rb_gsl_vector_pow_bang, 1);  
  
  Init_gsl_vector_init(module);
}
