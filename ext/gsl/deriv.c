/*
  deriv.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#ifdef GSL_1_4_9_LATER
#include "include/rb_gsl_common.h"
#include "include/rb_gsl_function.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#define RB_GSL_DERIV_H_DEFAULT (1e-8)
#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif

static int get_func2(int argc, VALUE *argv, VALUE obj, VALUE *ff, VALUE *xx, VALUE *hh)
{
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc == 3) {
      CHECK_FUNCTION(argv[0]);
      Need_Float(argv[2]);
      *ff = argv[0];
      *xx = argv[1];
      *hh = argv[2];
    } else if (argc == 2) {
      CHECK_FUNCTION(argv[0]);
      *ff = argv[0];
      *xx = argv[1];
      *hh = rb_float_new(RB_GSL_DERIV_H_DEFAULT);
    } else {
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    }
    break;
  default:
    if (argc == 2) {
      Need_Float(argv[1]);
      *ff = obj;
      *xx = argv[0];
      *hh = argv[1];
    } else if (argc == 1) {
      *ff = obj;
      *xx = argv[0];
      *hh = rb_float_new(RB_GSL_DERIV_H_DEFAULT);
    } else {
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    }
    break;
  }
  return 0;
}

#ifdef RB_GSL_DERIV_H_DEFAULT
#undef RB_GSL_DERIV_H_DEFAULT
#endif

static VALUE rb_gsl_deriv_eval(VALUE obj, VALUE xx, VALUE hh,
             int (*deriv)(const gsl_function *,
              double, double,
              double *, double *))
{
  gsl_function *f = NULL;
  double result, abserr, h;
  VALUE x, ary, aerr;
  gsl_vector *v = NULL, *vnew = NULL, *verr = NULL;
  gsl_matrix *m = NULL, *mnew = NULL, *merr = NULL;
  size_t n, i, j;
  int status;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2, *ptr3;
  VALUE ary2, ary3;
#endif
  Need_Float(hh);
  Data_Get_Struct(obj, gsl_function, f);
  h = NUM2DBL(hh);
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    status = (*deriv)(f, NUM2DBL(xx), h, &result, &abserr);
    return rb_ary_new3(3, rb_float_new(result), rb_float_new(abserr), INT2FIX(status));
    break;
  case T_ARRAY:
    //    n = RARRAY(xx)->len;
    n = RARRAY_LEN(xx);
    ary = rb_ary_new2(n);
    aerr = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(x);
      (*deriv)(f, NUM2DBL(x), h, &result, &abserr);
      rb_ary_store(ary, i, rb_float_new(result));
      rb_ary_store(aerr, i, rb_float_new(abserr));
    }
    return rb_ary_new3(2, ary, aerr);
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      GetNArray(xx, na);
      n = na->total;
      ptr1 = (double*) na->ptr;
      ary2 = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ary3 = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary2, double*);
      ptr3 = NA_PTR_TYPE(ary3, double*);
      for (i = 0; i < n; i++) {
        (*deriv)(f, ptr1[i], h, &result, &abserr);
        ptr2[i] = result;
        ptr3[i] = abserr;
      }
      return rb_ary_new3(2, ary2, ary3);
    }
#endif
    if (VECTOR_P(xx)) {
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      verr = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
        (*deriv)(f, gsl_vector_get(v, i), h, &result, &abserr);
        gsl_vector_set(vnew, i, result);
        gsl_vector_set(verr, i, abserr);
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
          (*deriv)(f, gsl_matrix_get(m, i, j), h, &result, &abserr);
          gsl_matrix_set(mnew, i, j, result);
          gsl_matrix_set(merr, i, j, abserr);
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
  return Qnil; /* never reach here */
}

static VALUE rb_gsl_deriv_central(int argc, VALUE *argv, VALUE obj)
{
  VALUE ff, xx, hh;
  get_func2(argc, argv, obj, &ff, &xx, &hh);
  return rb_gsl_deriv_eval(ff, xx, hh, gsl_deriv_central);
}

static VALUE rb_gsl_deriv_forward(int argc, VALUE *argv, VALUE obj)
{
  VALUE ff, xx, hh;
  get_func2(argc, argv, obj, &ff, &xx, &hh);
  return rb_gsl_deriv_eval(ff, xx, hh, gsl_deriv_forward);
}

static VALUE rb_gsl_deriv_backward(int argc, VALUE *argv, VALUE obj)
{
  VALUE ff, xx, hh;
  get_func2(argc, argv, obj, &ff, &xx, &hh);
  return rb_gsl_deriv_eval(ff, xx, hh, gsl_deriv_backward);
}

void Init_gsl_deriv(VALUE module)
{
  VALUE mgsl_deriv;

  mgsl_deriv = rb_define_module_under(module, "Deriv");

  rb_define_method(cgsl_function, "deriv_central", rb_gsl_deriv_central, -1);
  rb_define_method(cgsl_function, "deriv_forward", rb_gsl_deriv_forward, -1);
  rb_define_method(cgsl_function, "deriv_backward", rb_gsl_deriv_backward, -1);

  rb_define_singleton_method(mgsl_deriv, "central", rb_gsl_deriv_central, -1);
  rb_define_singleton_method(mgsl_deriv, "forward", rb_gsl_deriv_forward, -1);
  rb_define_singleton_method(mgsl_deriv, "backward", rb_gsl_deriv_backward, -1);
}

#endif
