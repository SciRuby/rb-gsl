/*
  interp.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_config.h"
#include "rb_gsl_interp.h"

VALUE cgsl_interp_accel; /* this is used also in spline.c */
EXTERN VALUE cgsl_vector, cgsl_matrix;

static void rb_gsl_interp_free(rb_gsl_interp *sp);

static VALUE rb_gsl_interp_new(int argc, VALUE *argv, VALUE klass)
{
  rb_gsl_interp *sp = NULL;
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
  if (size == 0) rb_raise(rb_eRuntimeError, "interp size is not given.");
  sp = ALLOC(rb_gsl_interp);
  if (T == NULL) T = gsl_interp_cspline;
  sp->p = gsl_interp_alloc(T, size);
  sp->a = gsl_interp_accel_alloc();
  if (ptrx && ptry) gsl_interp_init(sp->p, ptrx, ptry, size);
  return Data_Wrap_Struct(klass, 0, rb_gsl_interp_free, sp);
}

static void rb_gsl_interp_free(rb_gsl_interp *sp)
{
  gsl_interp_free(sp->p);
  gsl_interp_accel_free(sp->a);
  free((rb_gsl_interp *) sp);
}

static VALUE rb_gsl_interp_init(VALUE obj, VALUE xxa, VALUE yya)
{
  rb_gsl_interp *rgi = NULL;
  double *ptrx = NULL, *ptry = NULL;
  size_t size, stride;
  ptrx = get_vector_ptr(xxa, &stride, &size);
  ptry = get_vector_ptr(yya, &stride, &size);
  Data_Get_Struct(obj, rb_gsl_interp, rgi);
  gsl_interp_init(rgi->p, ptrx, ptry, size);
  return obj;
}

static VALUE rb_gsl_interp_name(VALUE obj)
{
  rb_gsl_interp *rgi = NULL;
  Data_Get_Struct(obj, rb_gsl_interp, rgi);
  return rb_str_new2(gsl_interp_name(rgi->p));
}

static VALUE rb_gsl_interp_min_size(VALUE obj)
{
  rb_gsl_interp *rgi = NULL;
  Data_Get_Struct(obj, rb_gsl_interp, rgi);
  return INT2FIX(gsl_interp_min_size(rgi->p));
}

static VALUE rb_gsl_interp_bsearch(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v = NULL;
  double x;
  size_t indexl, indexh;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    switch (argc) {
    case 2:
      CHECK_VECTOR(argv[0]);
      Need_Float(argv[1]);
      Data_Get_Struct(argv[0], gsl_vector, v);
      x = NUM2DBL(argv[1]);
      indexl = gsl_vector_get(v, 0);
      indexh = gsl_vector_get(v, v->size-1);
      break;
    case 4:
      CHECK_VECTOR(argv[0]);
      Need_Float(argv[1]); Need_Float(argv[2]); Need_Float(argv[3]);
      Data_Get_Struct(argv[0], gsl_vector, v);
      x = NUM2DBL(argv[1]);
      indexl = NUM2DBL(argv[2]);
      indexh = NUM2DBL(argv[3]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 4)", argc);
      break;
    }
    break;
  default:
    Data_Get_Struct(obj, gsl_vector, v);
    switch (argc) {
    case 1:
      Need_Float(argv[0]);
      x = NUM2DBL(argv[0]);
      indexl = gsl_vector_get(v, 0);
      indexh = gsl_vector_get(v, v->size-1);
      break;
    case 3:
      Need_Float(argv[0]); Need_Float(argv[1]); Need_Float(argv[2]);
      x = NUM2DBL(argv[0]);
      indexl = NUM2DBL(argv[1]);
      indexh = NUM2DBL(argv[2]);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 3)", argc);
      break;
    }
    break; 
  }
  return INT2FIX(gsl_interp_bsearch(v->data, x, indexl, indexh));
}

static VALUE rb_gsl_interp_accel(VALUE obj)
{
  rb_gsl_interp *rgi = NULL;
  Data_Get_Struct(obj, rb_gsl_interp, rgi);
  return Data_Wrap_Struct(cgsl_interp_accel, 0, NULL, rgi->a);
}

static VALUE rb_gsl_interp_find(VALUE obj, VALUE vv, VALUE xx)
{
  rb_gsl_interp *rgi = NULL;
  double *ptr = NULL, x;
  size_t size, stride;
  Need_Float(xx);
  Data_Get_Struct(obj, rb_gsl_interp, rgi);
  ptr = get_vector_ptr(vv, &stride, &size);
  x = NUM2DBL(xx);
  return INT2FIX(gsl_interp_accel_find(rgi->a, ptr, size, x));
}

static VALUE rb_gsl_interp_accel_find(VALUE obj, VALUE vv, VALUE xx)
{
  gsl_interp_accel *a = NULL;
  double x, *ptr = NULL;
  size_t size, stride;
  Need_Float(xx);
  Data_Get_Struct(obj, gsl_interp_accel, a);
  ptr = get_vector_ptr(vv, &stride, &size);
  Need_Float(xx);
  x = NUM2DBL(xx);
  return INT2FIX(gsl_interp_accel_find(a, ptr, size, x));
}

static VALUE rb_gsl_interp_evaluate(VALUE obj, VALUE xxa, VALUE yya, VALUE xx,
				    double (*eval)(const gsl_interp *, const double [], 
						   const double [], double, 
						   gsl_interp_accel *))
{
  rb_gsl_interp *rgi = NULL;
  double *ptrx = NULL, *ptry = NULL;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  VALUE ary, x;
  double val;
  size_t n, i, j, size, stridex, stridey;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na = NULL;
  double *ptrz = NULL, *ptr = NULL;
#endif
  Data_Get_Struct(obj, rb_gsl_interp, rgi);
  ptrx = get_vector_ptr(xxa, &stridex, &size);
  if (size != rgi->p->size ){
    rb_raise(rb_eTypeError, "size mismatch (xa:%d != %d)",  (int) size, (int) rgi->p->size);
  }
  ptry = get_vector_ptr(yya, &stridey, &size);
  if (size != rgi->p->size ){
    rb_raise(rb_eTypeError, "size mismatch (ya:%d != %d)", (int) size, (int) rgi->p->size);
  }
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:  case T_BIGNUM:  case T_FLOAT:
    Need_Float(xx);
    return rb_float_new((*eval)(rgi->p, ptrx, ptry, NUM2DBL(xx), rgi->a));
    break;
  case T_ARRAY:
    //    n = RARRAY(xx)->len;
    n = RARRAY_LEN(xx);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(x);
      val = (*eval)(rgi->p, ptrx, ptry, NUM2DBL(x), rgi->a);
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      GetNArray(xx, na);
      ptrz = (double*) na->ptr;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr = NA_PTR_TYPE(ary, double*);
      for (i = 0; (int) i < na->total; i++)
	ptr[i] = (*eval)(rgi->p, ptrx, ptry, ptrz[i], rgi->a);
      return ary;
    }
#endif
    if (VECTOR_P(xx)) {
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	val = (*eval)(rgi->p, ptrx, ptry, gsl_vector_get(v, i), rgi->a);
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MATRIX_P(xx)) {
      Data_Get_Struct(xx, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  val = (*eval)(rgi->p, ptrx, ptry, gsl_matrix_get(m, i, j), rgi->a);
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

static VALUE rb_gsl_interp_eval(VALUE obj, VALUE xxa, VALUE yya, VALUE xx)
{
  return rb_gsl_interp_evaluate(obj, xxa, yya, xx, gsl_interp_eval);
}

static VALUE rb_gsl_interp_eval_e(VALUE obj, VALUE xxa, VALUE yya, VALUE xx)
{
  rb_gsl_interp *rgi = NULL;
  double *ptr1 = NULL, *ptr2 = NULL;
  size_t size, stridex, stridey;
  double x, y;
  int status;
  Need_Float(xx);
  Data_Get_Struct(obj, rb_gsl_interp, rgi);
  ptr1 = get_vector_ptr(xxa, &stridex, &size);
  ptr2 = get_vector_ptr(yya, &stridey, &size);
  x = NUM2DBL(xx);
  status = gsl_interp_eval_e(rgi->p, ptr1, ptr2, x, rgi->a, &y);
  switch (status) {
  case GSL_EDOM:
    rb_gsl_error_handler("gsl_interp_eval_e error", __FILE__, __LINE__, status);
    break;
  default:
    return rb_float_new(y);
    break;
  }
  return Qnil;
}

static VALUE rb_gsl_interp_eval_deriv(VALUE obj, VALUE xxa, VALUE yya, VALUE xx)
{
  return rb_gsl_interp_evaluate(obj, xxa, yya, xx, gsl_interp_eval_deriv);
}

static VALUE rb_gsl_interp_eval_deriv_e(VALUE obj, VALUE xxa, VALUE yya, VALUE xx)
{
  rb_gsl_interp *rgi = NULL;
  double *ptr1 = NULL, *ptr2 = NULL;
  size_t size, stridex, stridey;
  double x, y;
  int status;
  Need_Float(xx);
  Data_Get_Struct(obj, rb_gsl_interp, rgi);
  ptr1 = get_vector_ptr(xxa, &stridex, &size);
  ptr2 = get_vector_ptr(yya, &stridey, &size);
  x = NUM2DBL(xx);
  status = gsl_interp_eval_deriv_e(rgi->p, ptr1, ptr2, x, rgi->a, &y);
  switch (status) {
  case GSL_EDOM:
    rb_gsl_error_handler("gsl_interp_eval_deriv_e error", __FILE__, __LINE__, status);
    break;
  default:
    return rb_float_new(y);
    break;
  } 
  return Qnil;
}

static VALUE rb_gsl_interp_eval_deriv2(VALUE obj, VALUE xxa, VALUE yya, VALUE xx)
{
  return rb_gsl_interp_evaluate(obj, xxa, yya, xx, gsl_interp_eval_deriv2);
}

static VALUE rb_gsl_interp_eval_deriv2_e(VALUE obj, VALUE xxa, VALUE yya, VALUE xx)
{
  rb_gsl_interp *rgi = NULL;
  double *ptr1 = NULL, *ptr2 = NULL, x, y;
  size_t size, stridex, stridey;
  int status;
  Need_Float(xx);
  Data_Get_Struct(obj, rb_gsl_interp, rgi);
  ptr1 = get_vector_ptr(xxa, &stridex, &size);
  ptr2 = get_vector_ptr(yya, &stridey, &size);
  x = NUM2DBL(xx);
  status = gsl_interp_eval_deriv2_e(rgi->p, ptr1, ptr2, x, rgi->a, &y);
  switch (status) {
  case GSL_EDOM:
    rb_gsl_error_handler("gsl_interp_eval_deriv2_e error", __FILE__, __LINE__, status);
    break;
  default:
    return rb_float_new(y);
    break;
  }
  return Qnil;
}

static VALUE rb_gsl_interp_eval_integ(VALUE obj, VALUE xxa, VALUE yya, 
				      VALUE aa, VALUE bb)
{
  rb_gsl_interp *rgi = NULL;
  double *ptr1 = NULL, *ptr2 = NULL;
  size_t size, stridex, stridey;
  double a, b;
  Need_Float(aa); Need_Float(bb);
  Data_Get_Struct(obj, rb_gsl_interp, rgi);
  ptr1 = get_vector_ptr(xxa, &stridex, &size);
  ptr2 = get_vector_ptr(yya, &stridey, &size);
  a = NUM2DBL(aa);
  b = NUM2DBL(bb);
  return rb_float_new(gsl_interp_eval_integ(rgi->p, ptr1, ptr2, a, b, rgi->a));
}

static VALUE rb_gsl_interp_eval_integ_e(VALUE obj, VALUE xxa, VALUE yya, 
					VALUE aa, VALUE bb)
{
  rb_gsl_interp *rgi = NULL;
  double *ptr1 = NULL, *ptr2 = NULL;
  size_t size, stridex, stridey;
  double y, a, b;
  int status;
  Need_Float(aa);
  Need_Float(bb);
  Data_Get_Struct(obj, rb_gsl_interp, rgi);
  ptr1 = get_vector_ptr(xxa, &stridex, &size);
  ptr2 = get_vector_ptr(yya, &stridey, &size);
  a = NUM2DBL(aa);
  b = NUM2DBL(bb);
  status = gsl_interp_eval_integ_e(rgi->p, ptr1, ptr2, a, b, rgi->a, &y);
  switch (status) {
  case GSL_EDOM:
    rb_gsl_error_handler("gsl_interp_eval_integ_e error", __FILE__, __LINE__, status);
    break;
  default:
    return rb_float_new(y);
    break;
  }
  return Qnil;
}

/******/

const gsl_interp_type* get_interp_type(VALUE t)
{
  int type;
  char name[32];
  switch (TYPE(t)) {
  case T_FIXNUM:
    type = FIX2INT(t);
    switch (type) {
    case GSL_INTERP_LINEAR: return gsl_interp_linear; break;
#ifdef GSL_1_1_LATER
    case GSL_INTERP_POLYNOMIAL: return gsl_interp_polynomial; break;
#endif
    case GSL_INTERP_CSPLINE: return gsl_interp_cspline; break;
    case GSL_INTERP_CSPLINE_PERIODIC: return gsl_interp_cspline_periodic; break;
    case GSL_INTERP_AKIMA: return gsl_interp_akima; break;
    case GSL_INTERP_AKIMA_PERIODIC: return gsl_interp_akima_periodic; break;
    default:
      rb_raise(rb_eTypeError, "unknown type %d\n", type);
      break;
    }
    break;
  case T_STRING:
    strcpy(name, STR2CSTR(t));
    if (str_tail_grep(name, "linear") == 0) {
      return gsl_interp_linear;
#ifdef GSL_1_1_LATER
    } else if (str_tail_grep(name, "polynomial") == 0) {
      return gsl_interp_polynomial;
#endif
    } else if (str_tail_grep(name, "cspline") == 0) {
      return gsl_interp_cspline;
    } else if (str_tail_grep(name, "cspline_periodic") == 0) {
      return gsl_interp_cspline_periodic;
    } else if (str_tail_grep(name, "akima") == 0) {
      return gsl_interp_akima;
    } else if (str_tail_grep(name, "akima_periodic") == 0) {
      return gsl_interp_akima_periodic;
    } else {
      rb_raise(rb_eTypeError, "Unknown type");
    }
    break;
  default:
    rb_raise(rb_eTypeError, "Unknown type");
    break;
  }
}

static VALUE rb_gsl_interp_info(VALUE obj)
{
  rb_gsl_interp *p;
  char buf[256];
  Data_Get_Struct(obj, rb_gsl_interp, p);
  sprintf(buf, "Class:      %s\n", rb_class2name(CLASS_OF(obj)));
  sprintf(buf, "%sSuperClass: %s\n", buf, rb_class2name(RCLASS_SUPER(CLASS_OF(obj))));
  sprintf(buf, "%sType:       %s\n", buf, gsl_interp_name(p->p));
  sprintf(buf, "%sxmin:       %f\n", buf, p->p->xmin);
  sprintf(buf, "%sxmax:       %f\n", buf, p->p->xmax);
  sprintf(buf, "%sSize:       %d\n", buf, (int) p->p->size);
  return rb_str_new2(buf);
}

static void rb_gsl_interp_define_const(VALUE klass)
{
  rb_define_const(klass, "LINEAR", INT2FIX(GSL_INTERP_LINEAR));
  rb_define_const(klass, "CSPLINE", INT2FIX(GSL_INTERP_CSPLINE));
  rb_define_const(klass, "CSPLINE_PERIODIC", INT2FIX(GSL_INTERP_CSPLINE_PERIODIC));
  rb_define_const(klass, "AKIMA", INT2FIX(GSL_INTERP_AKIMA));
  rb_define_const(klass, "AKIMA_PERIODIC", INT2FIX(GSL_INTERP_AKIMA_PERIODIC));

  rb_define_const(klass, "Linear", INT2FIX(GSL_INTERP_LINEAR));
#ifdef GSL_1_1_LATER
  rb_define_const(klass, "POLYNOMIAL", INT2FIX(GSL_INTERP_POLYNOMIAL));
  rb_define_const(klass, "Polynomial", INT2FIX(GSL_INTERP_POLYNOMIAL));
#endif
  rb_define_const(klass, "CSpline", INT2FIX(GSL_INTERP_CSPLINE));
  rb_define_const(klass, "CSpline_Periodic", INT2FIX(GSL_INTERP_CSPLINE_PERIODIC));
  rb_define_const(klass, "Akima", INT2FIX(GSL_INTERP_AKIMA));
  rb_define_const(klass, "Akima_Periodic", INT2FIX(GSL_INTERP_AKIMA_PERIODIC));
}

void Init_gsl_interp(VALUE module)
{
  VALUE cgsl_interp;

  cgsl_interp = rb_define_class_under(module, "Interp", cGSL_Object);
  cgsl_interp_accel = rb_define_class_under(cgsl_interp, "Accel", cGSL_Object);

  rb_define_singleton_method(cgsl_interp, "alloc", rb_gsl_interp_new, -1);

  rb_gsl_interp_define_const(cgsl_interp);

  /*****/

  rb_define_singleton_method(cgsl_interp, "bsearch", rb_gsl_interp_bsearch, -1);
  rb_define_method(cgsl_vector, "bsearch", rb_gsl_interp_bsearch, -1);
  rb_define_method(cgsl_interp, "name", rb_gsl_interp_name, 0);
  rb_define_alias(cgsl_interp, "type", "name");
  rb_define_method(cgsl_interp, "min_size", rb_gsl_interp_min_size, 0);
  rb_define_method(cgsl_interp, "init", rb_gsl_interp_init, 2);
  rb_define_method(cgsl_interp, "accel", rb_gsl_interp_accel, 0);
  rb_define_method(cgsl_interp, "eval", rb_gsl_interp_eval, 3);
  rb_define_alias(cgsl_interp, "[]", "eval");
  rb_define_method(cgsl_interp, "eval_e", rb_gsl_interp_eval_e, 3);
  rb_define_method(cgsl_interp, "eval_deriv", rb_gsl_interp_eval_deriv, 3);
  rb_define_alias(cgsl_interp, "deriv", "eval_deriv");
  rb_define_method(cgsl_interp, "eval_deriv_e", rb_gsl_interp_eval_deriv_e, 3);
  rb_define_method(cgsl_interp, "eval_deriv2", rb_gsl_interp_eval_deriv2, 3);
  rb_define_alias(cgsl_interp, "deriv2", "eval_deriv2");
  rb_define_method(cgsl_interp, "eval_deriv2_e", rb_gsl_interp_eval_deriv2_e, 3);
  rb_define_method(cgsl_interp, "eval_integ", rb_gsl_interp_eval_integ, 4);
  rb_define_alias(cgsl_interp, "integ", "eval_integ");
  rb_define_method(cgsl_interp, "eval_integ_e", rb_gsl_interp_eval_integ_e, 4);


  /*****/
  rb_define_method(cgsl_interp_accel, "find", rb_gsl_interp_accel_find, 2);

  rb_define_method(cgsl_interp, "find", rb_gsl_interp_find, 2);
  rb_define_alias(cgsl_interp, "accel_find", "find");

  rb_define_method(cgsl_interp, "info", rb_gsl_interp_info, 0);
}
