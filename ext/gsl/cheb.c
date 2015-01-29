/*
  cheb.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl.h"
#include "include/rb_gsl_common.h"
#include "include/rb_gsl_array.h"
#include "include/rb_gsl_function.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>
#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif

static VALUE cgsl_cheb;

static VALUE rb_gsl_cheb_new(VALUE klass, VALUE nn)
{
  gsl_cheb_series *p = NULL;
  CHECK_FIXNUM(nn);
  p = gsl_cheb_alloc(FIX2INT(nn));
  return Data_Wrap_Struct(klass, 0, gsl_cheb_free, p);
}

static VALUE rb_gsl_cheb_order(VALUE obj)
{
  gsl_cheb_series *p = NULL;
  Data_Get_Struct(obj, gsl_cheb_series, p);
  return INT2FIX(p->order);
}

static VALUE rb_gsl_cheb_a(VALUE obj)
{
  gsl_cheb_series *p = NULL;
  Data_Get_Struct(obj, gsl_cheb_series, p);
  return rb_float_new(p->a);
}

static VALUE rb_gsl_cheb_b(VALUE obj)
{
  gsl_cheb_series *p = NULL;
  Data_Get_Struct(obj, gsl_cheb_series, p);
  return rb_float_new(p->b);
}

static VALUE rb_gsl_cheb_coef(VALUE obj)
{
  gsl_cheb_series *p = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_cheb_series, p);
  v = gsl_vector_view_alloc();
  v->vector.data = p->c;
  v->vector.size = p->order + 1;
  v->vector.stride = 1;
  v->vector.owner = 0;
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_cheb_f(VALUE obj)
{
  gsl_cheb_series *p = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_cheb_series, p);
  v = gsl_vector_view_alloc();
  v->vector.data = p->f;
  v->vector.size = p->order + 1;
  v->vector.stride = 1;
  v->vector.owner = 0;
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_cheb_init(VALUE obj, VALUE ff, VALUE aa, VALUE bb)
{
  gsl_cheb_series *p = NULL;
  gsl_function *fff = NULL;
  double a, b;
  CHECK_FUNCTION(ff);
  Need_Float(aa);  Need_Float(bb);
  Data_Get_Struct(obj, gsl_cheb_series, p);
  Data_Get_Struct(ff, gsl_function, fff);
  a = NUM2DBL(aa);
  b = NUM2DBL(bb);
  gsl_cheb_init(p, fff, a, b);
  return obj;
}

static VALUE rb_gsl_cheb_eval(VALUE obj, VALUE xx)
{
  gsl_cheb_series *p = NULL;
  VALUE x, ary;
  size_t i, j, n;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2;
#endif
  Data_Get_Struct(obj, gsl_cheb_series, p);
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new(gsl_cheb_eval(p, NUM2DBL(xx)));
    break;
  case T_ARRAY:
    //    n = RARRAY(xx)->len;
    n = RARRAY_LEN(xx);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(xx);
      rb_ary_store(ary, i, rb_float_new(gsl_cheb_eval(p, NUM2DBL(x))));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      GetNArray(xx, na);
      ptr1 = (double*) na->ptr;
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary,double*);
      for (i = 0; i < n; i++) ptr2[i] = gsl_cheb_eval(p, ptr1[i]);
      return ary;
    }
#endif
    if (VECTOR_P(xx)) {
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
        gsl_vector_set(vnew, i, gsl_cheb_eval(p, gsl_vector_get(v, i)));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MATRIX_P(xx)) {
      Data_Get_Struct(xx, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
        for (j = 0; j < m->size2; j++) {
          gsl_matrix_set(mnew, i, j, gsl_cheb_eval(p, gsl_matrix_get(m, i, j)));
        }
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  }
  return Qnil;   /* never reach here */
}

static VALUE rb_gsl_cheb_eval_err(VALUE obj, VALUE xx)
{
  gsl_cheb_series *p = NULL;
  double result, err;
  VALUE x, ary, aerr;
  size_t n, i, j;
  gsl_vector *v = NULL, *vnew = NULL, *verr = NULL;
  gsl_matrix *m = NULL, *mnew = NULL, *merr = NULL;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2, *ptr3;
#endif
  Data_Get_Struct(obj, gsl_cheb_series, p);
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    gsl_cheb_eval_err(p, NUM2DBL(xx), &result, &err);
    return rb_ary_new3(2, rb_float_new(result), rb_float_new(err));
    break;
  case T_ARRAY:
    //    n = RARRAY(xx)->len;
    n = RARRAY_LEN(xx);
    ary = rb_ary_new2(n);
    aerr = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(xx);
      gsl_cheb_eval_err(p, NUM2DBL(x), &result, &err);
      rb_ary_store(ary, i, rb_float_new(result));
      rb_ary_store(aerr, i, rb_float_new(err));
    }
    return rb_ary_new3(2, ary, aerr);
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      GetNArray(xx, na);
      ptr1 = (double*) na->ptr;
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      aerr = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary,double*);
      ptr3 = NA_PTR_TYPE(aerr,double*);
      for (i = 0; i < n; i++) {
        gsl_cheb_eval_err(p, ptr1[i], &result, &err);
        ptr2[i] = result;
        ptr3[i] = err;
      }
      return rb_ary_new3(2, ary, aerr);
    }
#endif
    if (VECTOR_P(xx)) {
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      verr = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
        gsl_cheb_eval_err(p, gsl_vector_get(v, i), &result, &err);
        gsl_vector_set(vnew, i, result);
        gsl_vector_set(verr, i, err);
      }
      return rb_ary_new3(2,
       Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew),
       Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, verr));
    } else if (MATRIX_P(xx)) {
      Data_Get_Struct(xx, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      merr = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
        for (j = 0; j < m->size2; j++) {
          gsl_cheb_eval_err(p, gsl_matrix_get(m, i, j), &result, &err);
          gsl_matrix_set(mnew, i, j, result);
          gsl_matrix_set(merr, i, j, err);
        }
      }
      return rb_ary_new3(2,
       Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew),
       Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, merr));
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  }
  return Qnil;   /* never reach here */
}

static VALUE rb_gsl_cheb_eval_n(VALUE obj, VALUE nn, VALUE xx)
{
  gsl_cheb_series *p = NULL;
  VALUE x, ary;
  size_t n, order, i, j;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2;
#endif
  CHECK_FIXNUM(nn);
  order = FIX2INT(nn);
  Data_Get_Struct(obj, gsl_cheb_series, p);
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new(gsl_cheb_eval_n(p, order, NUM2DBL(xx)));
    break;
  case T_ARRAY:
    //    n = RARRAY(xx)->len;
    n = RARRAY_LEN(xx);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(xx);
      rb_ary_store(ary, i, rb_float_new(gsl_cheb_eval_n(p, order, NUM2DBL(x))));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      GetNArray(xx, na);
      ptr1 = (double*) na->ptr;
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary,double*);
      for (i = 0; i < n; i++) ptr2[i] = gsl_cheb_eval_n(p, order, ptr1[i]);
      return ary;
    }
#endif
    if (VECTOR_P(xx)) {
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
        gsl_vector_set(vnew, i, gsl_cheb_eval_n(p, order, gsl_vector_get(v, i)));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MATRIX_P(xx)) {
      Data_Get_Struct(xx, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
        for (j = 0; j < m->size2; j++) {
          gsl_matrix_set(mnew, i, j,
            gsl_cheb_eval_n(p, order, gsl_matrix_get(m, i, j)));
        }
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  }
  return Qnil;   /* never reach here */
}

static VALUE rb_gsl_cheb_eval_n_err(VALUE obj, VALUE nn, VALUE xx)
{
  gsl_cheb_series *p = NULL;
  double result, err;
  VALUE x, ary, aerr;
  size_t n, order, i, j;
  gsl_vector *v, *vnew, *verr;
  gsl_matrix *m, *mnew, *merr;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2, *ptr3;
#endif
  CHECK_FIXNUM(nn);
  order = FIX2INT(nn);
  Data_Get_Struct(obj, gsl_cheb_series, p);
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    gsl_cheb_eval_n_err(p, order, NUM2DBL(xx), &result, &err);
    return rb_ary_new3(2, rb_float_new(result), rb_float_new(err));
    break;
  case T_ARRAY:
    //    n = RARRAY(xx)->len;
    n = RARRAY_LEN(xx);
    ary = rb_ary_new2(n);
    aerr = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(xx);
      gsl_cheb_eval_n_err(p, order, NUM2DBL(x), &result, &err);
      rb_ary_store(ary, i, rb_float_new(result));
      rb_ary_store(aerr, i, rb_float_new(err));
    }
    return rb_ary_new3(2, ary, aerr);
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      GetNArray(xx, na);
      ptr1 = (double*) na->ptr;
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      aerr = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary,double*);
      ptr3 = NA_PTR_TYPE(aerr,double*);
      for (i = 0; i < n; i++) {
        gsl_cheb_eval_n_err(p, order, ptr1[i], &result, &err);
        ptr2[i] = result;
        ptr3[i] = err;
      }
      return rb_ary_new3(2, ary, aerr);
    }
#endif
    if (VECTOR_P(xx)) {
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      verr = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
        gsl_cheb_eval_n_err(p, order, gsl_vector_get(v, i), &result, &err);
        gsl_vector_set(vnew, i, result);
        gsl_vector_set(verr, i, err);
      }
      return rb_ary_new3(2,
       Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew),
       Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, verr));
    } else if (MATRIX_P(xx)) {
      Data_Get_Struct(xx, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      merr = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
        for (j = 0; j < m->size2; j++) {
          gsl_cheb_eval_n_err(p, order, gsl_matrix_get(m, i, j), &result, &err);
          gsl_matrix_set(mnew, i, j, result);
          gsl_matrix_set(merr, i, j, err);
        }
      }
      return rb_ary_new3(2,
       Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew),
       Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, merr));
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  }
  return Qnil;   /* never reach here */
}

static VALUE rb_gsl_cheb_calc_deriv(int argc, VALUE *argv, VALUE obj)
{
  gsl_cheb_series *deriv = NULL, *cs = NULL;
  VALUE retval;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    switch (argc) {
    case 1:
      if (!rb_obj_is_kind_of(argv[0], cgsl_cheb))
        rb_raise(rb_eTypeError, "wrong argument type %s (Cheb expected)",
          rb_class2name(CLASS_OF(argv[0])));
      Data_Get_Struct(argv[0], gsl_cheb_series, cs);
      deriv = gsl_cheb_alloc(cs->order);
      retval = Data_Wrap_Struct(CLASS_OF(argv[0]), 0, gsl_cheb_free, deriv);
      break;
    case 2:
      if (!rb_obj_is_kind_of(argv[0], cgsl_cheb))
        rb_raise(rb_eTypeError, "argv[0] wrong argument type %s (Cheb expected)",
          rb_class2name(CLASS_OF(argv[0])));
      if (!rb_obj_is_kind_of(argv[1], cgsl_cheb))
        rb_raise(rb_eTypeError, "argv[1] wrong argument type %s (Cheb expected)",
          rb_class2name(CLASS_OF(argv[1])));
      Data_Get_Struct(argv[0], gsl_cheb_series, deriv);
      Data_Get_Struct(argv[1], gsl_cheb_series, cs);
      retval = argv[0];
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
      break;
    }
    break;
  default:
    Data_Get_Struct(obj, gsl_cheb_series, cs);
    switch (argc) {
    case 0:
      deriv = gsl_cheb_alloc(cs->order);
      retval = Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_cheb_free, deriv);
      break;
    case 1:
      if (!rb_obj_is_kind_of(argv[0], cgsl_cheb))
        rb_raise(rb_eTypeError, "argv[0] wrong argument type %s (Cheb expected)",
          rb_class2name(CLASS_OF(argv[0])));
      Data_Get_Struct(argv[0], gsl_cheb_series, deriv);
      retval = argv[0];
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
      break;
    }
    break;
  }
  gsl_cheb_calc_deriv(deriv, cs);
  return retval;
}

static VALUE rb_gsl_cheb_calc_integ(int argc, VALUE *argv, VALUE obj)
{
  gsl_cheb_series *deriv = NULL, *cs = NULL;
  VALUE retval;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    switch (argc) {
    case 1:
      if (!rb_obj_is_kind_of(argv[0], cgsl_cheb))
        rb_raise(rb_eTypeError, "wrong argument type %s (Cheb expected)",
          rb_class2name(CLASS_OF(argv[0])));
      Data_Get_Struct(argv[0], gsl_cheb_series, cs);
      deriv = gsl_cheb_alloc(cs->order);
      retval = Data_Wrap_Struct(CLASS_OF(argv[0]), 0, gsl_cheb_free, deriv);
      break;
    case 2:
      if (!rb_obj_is_kind_of(argv[0], cgsl_cheb))
        rb_raise(rb_eTypeError, "argv[0] wrong argument type %s (Cheb expected)",
          rb_class2name(CLASS_OF(argv[0])));
      if (!rb_obj_is_kind_of(argv[1], cgsl_cheb))
        rb_raise(rb_eTypeError, "argv[1] wrong argument type %s (Cheb expected)",
          rb_class2name(CLASS_OF(argv[1])));
      Data_Get_Struct(argv[0], gsl_cheb_series, deriv);
      Data_Get_Struct(argv[1], gsl_cheb_series, cs);
      retval = argv[0];
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
      break;
    }
    break;
  default:
    Data_Get_Struct(obj, gsl_cheb_series, cs);
    switch (argc) {
    case 0:
      deriv = gsl_cheb_alloc(cs->order);
      retval = Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_cheb_free, deriv);
      break;
    case 1:
      if (!rb_obj_is_kind_of(argv[0], cgsl_cheb))
        rb_raise(rb_eTypeError, "argv[0] wrong argument type %s (Cheb expected)",
          rb_class2name(CLASS_OF(argv[0])));
      Data_Get_Struct(argv[0], gsl_cheb_series, deriv);
      retval = argv[0];
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
      break;
    }
    break;
  }
  gsl_cheb_calc_integ(deriv, cs);
  return retval;
}

void Init_gsl_cheb(VALUE module)
{
  cgsl_cheb = rb_define_class_under(module, "Cheb", cGSL_Object);
  rb_define_singleton_method(cgsl_cheb, "new", rb_gsl_cheb_new, 1);
  rb_define_singleton_method(cgsl_cheb, "alloc", rb_gsl_cheb_new, 1);
  rb_define_method(cgsl_cheb, "order", rb_gsl_cheb_order, 0);
  rb_define_method(cgsl_cheb, "a", rb_gsl_cheb_a, 0);
  rb_define_method(cgsl_cheb, "b", rb_gsl_cheb_b, 0);
  rb_define_method(cgsl_cheb, "coef", rb_gsl_cheb_coef, 0);
  rb_define_alias(cgsl_cheb, "c", "coef");
  rb_define_method(cgsl_cheb, "f", rb_gsl_cheb_f, 0);
  rb_define_method(cgsl_cheb, "init", rb_gsl_cheb_init, 3);
  rb_define_method(cgsl_cheb, "eval", rb_gsl_cheb_eval, 1);
  rb_define_method(cgsl_cheb, "eval_err", rb_gsl_cheb_eval_err, 1);
  rb_define_method(cgsl_cheb, "eval_n", rb_gsl_cheb_eval_n, 2);
  rb_define_method(cgsl_cheb, "eval_n_err", rb_gsl_cheb_eval_n_err, 2);
  rb_define_method(cgsl_cheb, "calc_deriv", rb_gsl_cheb_calc_deriv, -1);
  rb_define_alias(cgsl_cheb, "deriv", "calc_deriv");
  rb_define_method(cgsl_cheb, "calc_integ", rb_gsl_cheb_calc_integ, -1);
  rb_define_alias(cgsl_cheb, "integ", "calc_integ");

  rb_define_singleton_method(cgsl_cheb, "calc_deriv", rb_gsl_cheb_calc_deriv, -1);
  rb_define_singleton_method(cgsl_cheb, "calc_integ", rb_gsl_cheb_calc_integ, -1);
}
