/*
  spline2d.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/
#ifdef GSL_2_0_LATER
#include "include/rb_gsl_interp2d.h"

EXTERN VALUE cgsl_interp2d_accel;  /* defined in interp2d.c */
static void rb_gsl_spline2d_free(rb_gsl_spline2d *fr);

static VALUE rb_gsl_spline2d_alloc(int argc, VALUE *argv, VALUE self)
{
  rb_gsl_spline2d *sp = NULL;
  const gsl_interp2d_type *T = NULL;
  double *xptr = NULL, *yptr = NULL, *zptr = NULL;
  size_t sizex = 0, sizey = 0, sizez = 0, stride = 1;

  T = get_interp2d_type(argv[0]);
  if (argc == 3) { // type, sizex, sizey
    sizex = FIX2INT(argv[1]);
    sizey = FIX2INT(argv[2]);
  }
  else if (argc == 4) { // type, xarr, yarr, zarr
    xptr = get_vector_ptr(argv[1], &stride, &sizex);
    yptr = get_vector_ptr(argv[2], &stride, &sizey);
    zptr = get_vector_ptr(argv[3], &stride, &sizez);
  }
  else {
    rb_raise(rb_eRuntimeError, "Expected args: (interp2d type, xn, yn) \
      or (interp2d type, xarr, yarr, zarr");
  }

  if (sizex == 0 || sizey == 0) rb_raise(rb_eRuntimeError, "Spline2d size not given.");
  sp = ALLOC(rb_gsl_spline2d);
  sp->s = gsl_spline2d_alloc(T, sizex, sizey);
  sp->xacc = gsl_interp_accel_alloc();
  sp->yacc = gsl_interp_accel_alloc();

  if (xptr && yptr && zptr) gsl_spline2d_init(sp->s, xptr, yptr, zptr, sizex, sizey);

  return Data_Wrap_Struct(self, 0, rb_gsl_spline2d_free, sp);
}

static void rb_gsl_spline2d_free(rb_gsl_spline2d *fr)
{
  gsl_spline2d_free(fr->s);
  gsl_interp_accel_free(fr->xacc);
  gsl_interp_accel_free(fr->yacc);
  free((rb_gsl_interp2d*) fr);
}

static VALUE rb_gsl_spline2d_init(VALUE self, VALUE xarr, VALUE yarr, VALUE zarr)
{
  rb_gsl_spline2d *rgs = NULL;
  double *xptr = NULL, *yptr = NULL, *zptr = NULL;
  size_t xsize, ysize, zsize, stride;
  xptr = get_vector_ptr(xarr, &stride, &xsize);
  yptr = get_vector_ptr(yarr, &stride, &ysize);
  zptr = get_vector_ptr(zarr, &stride, &zsize);

  Data_Get_Struct(self, rb_gsl_spline2d, rgs);
  gsl_spline2d_init(rgs->s, xptr, yptr, zptr, xsize, ysize);

  return self;
}

static VALUE rb_gsl_spline2d_accel(VALUE self)
{
  rb_gsl_spline2d *rgs = NULL;
  VALUE ary = rb_ary_new();
  Data_Get_Struct(self, rb_gsl_spline2d, rgs);

  rb_ary_push(ary, Data_Wrap_Struct(cgsl_interp2d_accel, 0, NULL, rgs->xacc));
  rb_ary_push(ary, Data_Wrap_Struct(cgsl_interp2d_accel, 0, NULL, rgs->yacc));

  return ary;
}

static VALUE rb_gsl_spline2d_evaluate(VALUE self, VALUE xx, VALUE yy,
  double (*eval)(const gsl_spline2d *, double, double, gsl_interp_accel *, 
    gsl_interp_accel *))
{ 
  VALUE is_swapped = rb_cvar_get(CLASS_OF(self), rb_intern("@@swapped"));
  VALUE temp;

  if (is_swapped != Qnil && is_swapped == Qtrue) {
    temp = xx;
    xx = yy;
    yy = temp;
  }
  
  rb_gsl_spline2d *rgs = NULL;
  gsl_vector *vx = NULL, *vy = NULL, *vnew = NULL;
  gsl_matrix *mx = NULL, *my = NULL, *mnew = NULL;
  VALUE ary, x, y;
  double val;
  size_t i, j;

  Data_Get_Struct(self, rb_gsl_spline2d, rgs);

  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  if (CLASS_OF(yy) == rb_cRange) yy = rb_gsl_range2ary(yy);

  if (TYPE(xx) != TYPE(yy)) {
    rb_raise(rb_eTypeError,"xx and yy must be same type. xx = %d yy = %d.",
      TYPE(xx), TYPE(yy));
  }

  switch (TYPE(xx)) {
  case T_FIXNUM:  case T_BIGNUM:  case T_FLOAT:
    Need_Float(xx);
    return rb_float_new((*eval)(rgs->s, NUM2DBL(xx), NUM2DBL(yy),
      rgs->xacc, rgs->yacc));
  case T_ARRAY:
    if (RARRAY_LEN(xx) != RARRAY_LEN(yy)) {
      rb_raise(rb_eRuntimeError, "xx and yy must be same sized Array.");
    }
    ary = rb_ary_new2(RARRAY_LEN(xx));

    for (i = 0; i < (unsigned)RARRAY_LEN(xx); i++) {
      x = rb_ary_entry(xx, i);
      y = rb_ary_entry(yy, i);
      Need_Float(x);
      Need_Float(y);

      val = (*eval)(rgs->s, NUM2DBL(x), NUM2DBL(y), rgs->xacc, rgs->yacc);
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
  default:
    if (VECTOR_P(xx)) {
      Data_Get_Struct(xx, gsl_vector, vx);
      Data_Get_Struct(yy, gsl_vector, vy);
      if (vx->size != vy->size) {
        rb_raise(rb_eRuntimeError, "xx and yy must be same size Vectors.");
      }
      vnew = gsl_vector_alloc(vx->size);

      for (i = 0; i < vx->size; i++) {
        val = (*eval)(rgs->s, gsl_vector_get(vx, i), gsl_vector_get(vy, i), 
          rgs->xacc, rgs->yacc);
        gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MATRIX_P(xx)) {
      Data_Get_Struct(xx, gsl_matrix, mx);
      Data_Get_Struct(xx, gsl_matrix, my);

      if ((mx->size1 != my->size1) || (mx->size2 != my->size2)) {
        rb_raise(rb_eRuntimeError, "xx and yy must be same sized Matrices.");
      }
      mnew = gsl_matrix_alloc(mx->size1, mx->size2);
      for (i = 0; i < mx->size1; i++) {
        for (j = 0; j < mx->size2; j++) {
          val = (*eval)(rgs->s, gsl_matrix_get(mx, i, j),
            gsl_matrix_get(my, i, j), rgs->xacc, rgs->yacc);
          gsl_matrix_set(mnew, i, j, val);
        }
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s", rb_class2name(CLASS_OF(xx)));
    }
  }
}

static VALUE rb_gsl_spline2d_eval(VALUE self, VALUE xx, VALUE yy)
{ 
  return rb_gsl_spline2d_evaluate(self, xx, yy, gsl_spline2d_eval);
}

static VALUE rb_gsl_spline2d_info(VALUE self)
{
  rb_gsl_spline2d *p = NULL;
  char buf[256];
  Data_Get_Struct(self, rb_gsl_spline2d, p);
  sprintf(buf, "Class:      %s\n", rb_class2name(CLASS_OF(self)));
  sprintf(buf, "%sSuperClass: %s\n", buf, rb_class2name(RCLASS_SUPER(CLASS_OF(self))));
  sprintf(buf, "%sType:       %s\n", buf, gsl_interp2d_name(&p->s->interp_object));
  sprintf(buf, "%sxmin:       %f\n", buf, p->s->interp_object.xmin);
  sprintf(buf, "%sxmax:       %f\n", buf, p->s->interp_object.xmax);
  sprintf(buf, "%symin:       %f\n", buf, p->s->interp_object.ymin);
  sprintf(buf, "%symax:       %f\n", buf, p->s->interp_object.ymax);
  sprintf(buf, "%sxSize:       %d\n", buf, (int) p->s->interp_object.xsize);
  sprintf(buf, "%sySize:       %d\n", buf, (int) p->s->interp_object.ysize);
  return rb_str_new2(buf);
}

void Init_gsl_spline2d(VALUE module)
{
  VALUE cgsl_spline2d;

  cgsl_spline2d = rb_define_class_under(module, "Spline2d", cGSL_Object);

  rb_define_singleton_method(cgsl_spline2d, "alloc", rb_gsl_spline2d_alloc, -1);

  /*****/

  rb_define_method(cgsl_spline2d, "init", rb_gsl_spline2d_init, 3);
  rb_define_method(cgsl_spline2d, "accel", rb_gsl_spline2d_accel, 0);
  rb_define_method(cgsl_spline2d, "eval", rb_gsl_spline2d_eval, 2);
  rb_define_alias(cgsl_spline2d, "[]", "eval");
  rb_define_method(cgsl_spline2d, "info", rb_gsl_spline2d_info, 0);
}
#endif
