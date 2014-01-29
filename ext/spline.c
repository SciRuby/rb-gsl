/*
  spline.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_config.h"
#include "rb_gsl_interp.h"
EXTERN VALUE cgsl_interp_accel;  /* defined in interp.c */

static void rb_gsl_spline_free(rb_gsl_spline *sp);

static VALUE rb_gsl_spline_new(int argc, VALUE *argv, VALUE klass)
{
  rb_gsl_spline *sp = NULL;
  const gsl_interp_type *T = NULL;
  double *ptrx = NULL, *ptry = NULL;
  size_t sizex = 0, sizey = 0, size = 0, stride = 1;
  int i;
  for (i = 0; i < argc; i++) {
    switch (TYPE(argv[i])) {
    case T_STRING:
      T = get_interp_type(argv[i]);
      break;
    case T_FIXNUM:
      if (T) size = FIX2INT(argv[i]);
      else T = get_interp_type(argv[i]);
      break;
    default:
      if (ptrx == NULL) {
	ptrx = get_vector_ptr(argv[i], &stride, &sizex);
      } else {
	ptry = get_vector_ptr(argv[i], &stride, &sizey);
	size = GSL_MIN_INT(sizex, sizey);
      }
      break;
    }
  }
  if (size == 0) rb_raise(rb_eRuntimeError, "spline size is not given.");
  sp = ALLOC(rb_gsl_spline);
  if (T == NULL) T = gsl_interp_cspline;
  sp->s = gsl_spline_alloc(T, size);
  sp->a = gsl_interp_accel_alloc();
  if (ptrx && ptry) gsl_spline_init(sp->s, ptrx, ptry, size);
  return Data_Wrap_Struct(klass, 0, rb_gsl_spline_free, sp);
}

static void rb_gsl_spline_free(rb_gsl_spline *sp)
{
  gsl_spline_free(sp->s);
  gsl_interp_accel_free(sp->a);
  free((rb_gsl_spline *) sp);
}

static VALUE rb_gsl_spline_init(VALUE obj, VALUE xxa, VALUE yya)
{
  rb_gsl_spline *sp = NULL;
  gsl_spline *p = NULL;
  gsl_vector *xa = NULL, *ya = NULL;
  size_t i, size;
  int flagx = 0, flagy = 0;
  double *ptr1 = NULL, *ptr2 = NULL;
#ifdef HAVE_NARRAY_H
  struct NARRAY *nax = NULL, *nay = NULL;
#endif
  Data_Get_Struct(obj, rb_gsl_spline, sp);
  p = sp->s;
  if (TYPE(xxa) == T_ARRAY) {
    //    size = RARRAY(xxa)->len;
    size = RARRAY_LEN(xxa);
    xa = gsl_vector_alloc(size);
    for (i = 0; i < size; i++) gsl_vector_set(xa, i, NUM2DBL(rb_ary_entry(xxa, i)));
    ptr1 = xa->data;
    flagx = 1;
  } else if (VECTOR_P(xxa)) {
    Data_Get_Struct(xxa, gsl_vector, xa);
    size = xa->size;
    ptr1 = xa->data;
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(xxa)) {
      GetNArray(xxa, nax);
      size = nax->total;
      ptr1 = (double *) nax->ptr;
#endif
  } else {
    rb_raise(rb_eTypeError, "not a vector");
  }
  if (TYPE(yya) == T_ARRAY) {
    ya = gsl_vector_alloc(size);
    for (i = 0; i < size; i++) gsl_vector_set(ya, i, NUM2DBL(rb_ary_entry(yya, i)));
    ptr2 = ya->data;
    flagy = 1;
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(yya)) {
      GetNArray(yya, nay);
      ptr2 = (double *) nay->ptr;
#endif
  } else if (VECTOR_P(yya)) {
    Data_Get_Struct(yya, gsl_vector, ya);
    ptr2 = ya->data;
  } else {
    rb_raise(rb_eTypeError, "not a vector");
  }
  gsl_spline_init(p, ptr1, ptr2, size);
  if (flagx == 1) gsl_vector_free(xa);
  if (flagy == 1) gsl_vector_free(ya);
  return obj;
}

static VALUE rb_gsl_spline_accel(VALUE obj)
{
  rb_gsl_spline *rgi = NULL;
  Data_Get_Struct(obj, rb_gsl_spline, rgi);
  return Data_Wrap_Struct(cgsl_interp_accel, 0, NULL, rgi->a);
}

static VALUE rb_gsl_spline_evaluate(VALUE obj, VALUE xx,
				    double (*eval)(const gsl_spline *, double, 
						   gsl_interp_accel *))
{
  rb_gsl_spline *rgs = NULL;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  VALUE ary, x;
  double val;
  size_t n, i, j;
#ifdef HAVE_NARRAY_H
  double *ptr1 = NULL, *ptr2 = NULL;
  struct NARRAY *na = NULL;
#endif
  Data_Get_Struct(obj, rb_gsl_spline, rgs);
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:  case T_BIGNUM:  case T_FLOAT:
    Need_Float(xx);
    return rb_float_new((*eval)(rgs->s, NUM2DBL(xx), rgs->a));
    break;
  case T_ARRAY:
    //    n = RARRAY(xx)->len;
    n = RARRAY_LEN(xx);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(x);
      val = (*eval)(rgs->s, NUM2DBL(x), rgs->a);
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      GetNArray(xx, na);
      ptr1 = (double *) na->ptr;
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) 
	ptr2[i] = (*eval)(rgs->s, ptr1[i], rgs->a);
      return ary;
    }
#endif
    if (VECTOR_P(xx)) {
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	val = (*eval)(rgs->s, gsl_vector_get(v, i), rgs->a);
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MATRIX_P(xx)) {
      Data_Get_Struct(xx, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  val = (*eval)(rgs->s, gsl_matrix_get(m, i, j), rgs->a);
	  gsl_matrix_set(mnew, i, j, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s", rb_class2name(CLASS_OF(xx)));
    }
    break;
  }

  /* never reach here */
  return Qnil;
}

static VALUE rb_gsl_spline_eval(VALUE obj, VALUE xx)
{
  return rb_gsl_spline_evaluate(obj, xx, gsl_spline_eval);
}

static VALUE rb_gsl_spline_eval_deriv(VALUE obj, VALUE xx)
{
  return rb_gsl_spline_evaluate(obj, xx, gsl_spline_eval_deriv);
}

static VALUE rb_gsl_spline_eval_deriv2(VALUE obj, VALUE xx)
{
  return rb_gsl_spline_evaluate(obj, xx, gsl_spline_eval_deriv2);
}

static VALUE rb_gsl_spline_eval_integ(VALUE obj, VALUE aa, VALUE bb)
{
  rb_gsl_spline *sp = NULL;
  gsl_spline *s = NULL;
  gsl_interp_accel *acc = NULL;
  double a, b;
  Need_Float(aa);
  Need_Float(bb);
  Data_Get_Struct(obj, rb_gsl_spline, sp);
  s = sp->s;
  acc = sp->a;
  a = NUM2DBL(aa);
  b = NUM2DBL(bb);
  return rb_float_new(gsl_spline_eval_integ(s, a, b, acc));
}

static VALUE rb_gsl_spline_find(VALUE obj, VALUE vv, VALUE xx)
{
  rb_gsl_spline *sp = NULL;
  double *ptr = NULL, x;
  size_t size, stride;
  Data_Get_Struct(obj, rb_gsl_spline, sp);
  ptr = get_vector_ptr(vv, &stride, &size);
  //  x = RFLOAT(xx)->value;
  x = NUM2DBL(xx);
  return INT2FIX(gsl_interp_accel_find(sp->a, ptr, size, x));
}

static VALUE rb_gsl_spline_eval_e(VALUE obj, VALUE xx)
{
  rb_gsl_spline *rgs = NULL;
  double val;
  int status;
  Data_Get_Struct(obj, rb_gsl_spline, rgs);
  Need_Float(xx);
  status = gsl_spline_eval_e(rgs->s, NUM2DBL(xx), rgs->a, &val);
  switch (status) {
  case GSL_EDOM:
    rb_gsl_error_handler("gsl_spline_eval_e error", __FILE__, __LINE__, status);
    break;
  default:
    return rb_float_new(val);
    break;
  }
  return Qnil;
}

static VALUE rb_gsl_spline_eval_deriv_e(VALUE obj, VALUE xx)
{
  rb_gsl_spline *rgs = NULL;
  double val;
  int status;
  Data_Get_Struct(obj, rb_gsl_spline, rgs);
  Need_Float(xx);
  status = gsl_spline_eval_deriv_e(rgs->s, NUM2DBL(xx), rgs->a, &val);
  switch (status) {
  case GSL_EDOM:
    rb_gsl_error_handler("gsl_spline_eval_deriv_e error", __FILE__, __LINE__, status);
    break;
  default:
    return rb_float_new(val);
    break;
  } 
  return Qnil;
}

static VALUE rb_gsl_spline_eval_deriv2_e(VALUE obj, VALUE xx)
{
  rb_gsl_spline *rgs = NULL;
  double val;
  int status;
  Data_Get_Struct(obj, rb_gsl_spline, rgs);
  Need_Float(xx);
  status = gsl_spline_eval_deriv2_e(rgs->s, NUM2DBL(xx), rgs->a, &val);
  switch (status) {
  case GSL_EDOM:
    rb_gsl_error_handler("gsl_spline_eval_deriv2_e error", __FILE__, __LINE__, status);
    break;
  default:
    return rb_float_new(val);
    break;
  } 
  return Qnil;
}

static VALUE rb_gsl_spline_eval_integ_e(VALUE obj, VALUE a, VALUE b)
{
  rb_gsl_spline *rgs = NULL;
  double val;
  int status;
  Data_Get_Struct(obj, rb_gsl_spline, rgs);
  Need_Float(a); Need_Float(b);
  status = gsl_spline_eval_integ_e(rgs->s, NUM2DBL(a), NUM2DBL(b), rgs->a, &val);
  switch (status) {
  case GSL_EDOM:
    rb_gsl_error_handler("gsl_spline_eval_integ_e error", __FILE__, __LINE__, status);
    break;
  default:
    return rb_float_new(val);
    break;
  } 
  return Qnil;
}

static VALUE rb_gsl_spline_info(VALUE obj)
{
  rb_gsl_spline *p = NULL;
  char buf[256];
  Data_Get_Struct(obj, rb_gsl_spline, p);
  sprintf(buf, "Class:      %s\n", rb_class2name(CLASS_OF(obj)));
  sprintf(buf, "%sSuperClass: %s\n", buf, rb_class2name(RCLASS_SUPER(CLASS_OF(obj))));
  sprintf(buf, "%sType:       %s\n", buf, gsl_interp_name(p->s->interp));
  sprintf(buf, "%sxmin:       %f\n", buf, p->s->interp->xmin);
  sprintf(buf, "%sxmax:       %f\n", buf, p->s->interp->xmax);
  sprintf(buf, "%sSize:       %d\n", buf, (int) p->s->size);
  return rb_str_new2(buf);
}

#ifdef GSL_1_8_LATER
static VALUE rb_gsl_spline_name(VALUE obj)
{
  rb_gsl_spline *p = NULL;
  Data_Get_Struct(obj, rb_gsl_spline, p);
  return rb_str_new2(gsl_spline_name(p->s));
}
static VALUE rb_gsl_spline_min_size(VALUE obj)
{
  rb_gsl_spline *sp = NULL;
  Data_Get_Struct(obj, rb_gsl_spline, sp);
  return UINT2NUM(gsl_spline_min_size(sp->s));
}

#else

static VALUE rb_gsl_spline_name(VALUE obj)
{
  rb_gsl_spline *sp = NULL;
  Data_Get_Struct(obj, rb_gsl_spline, sp);
  return rb_str_new2(gsl_interp_name(sp->s->interp));
}

#endif

void Init_gsl_spline(VALUE module)
{
  VALUE cgsl_spline;

  cgsl_spline = rb_define_class_under(module, "Spline", cGSL_Object);

  rb_define_singleton_method(cgsl_spline, "alloc", rb_gsl_spline_new, -1);

  /*****/

  rb_define_method(cgsl_spline, "init", rb_gsl_spline_init, 2);
  rb_define_method(cgsl_spline, "accel", rb_gsl_spline_accel, 0);
  rb_define_method(cgsl_spline, "eval", rb_gsl_spline_eval, 1);
  rb_define_alias(cgsl_spline, "[]", "eval");
  rb_define_method(cgsl_spline, "eval_deriv", rb_gsl_spline_eval_deriv, 1);
  rb_define_alias(cgsl_spline, "deriv", "eval_deriv");
  rb_define_method(cgsl_spline, "eval_deriv2", rb_gsl_spline_eval_deriv2, 1);
  rb_define_alias(cgsl_spline, "deriv2", "eval_deriv2");
  rb_define_method(cgsl_spline, "eval_integ", rb_gsl_spline_eval_integ, 2);
  rb_define_alias(cgsl_spline, "integ", "eval_integ");

  rb_define_method(cgsl_spline, "name", rb_gsl_spline_name, 0);
  rb_define_alias(cgsl_spline, "type", "name");

  rb_define_method(cgsl_spline, "find", rb_gsl_spline_find, 2);
  rb_define_alias(cgsl_spline, "accel_find", "find");

  rb_define_method(cgsl_spline, "eval_e", rb_gsl_spline_eval_e, 1);
  rb_define_method(cgsl_spline, "eval_deriv_e", rb_gsl_spline_eval_deriv_e, 1);
  rb_define_alias(cgsl_spline, "deriv_e", "eval_deriv_e");
  rb_define_method(cgsl_spline, "eval_deriv2_e", rb_gsl_spline_eval_deriv2_e, 1);
  rb_define_alias(cgsl_spline, "deri2v_e", "eval_deriv2_e");
  rb_define_method(cgsl_spline, "eval_integ_e", rb_gsl_spline_eval_integ_e, 1);
  rb_define_alias(cgsl_spline, "integ_e", "eval_integ_e");

  rb_define_method(cgsl_spline, "info", rb_gsl_spline_info, 0);

#ifdef GSL_1_8_LATER
  rb_define_method(cgsl_spline, "min_size", rb_gsl_spline_min_size, 0);
#endif

}
