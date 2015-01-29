/*
  diff.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_common.h"
#include "include/rb_gsl_function.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_diff.h>

static int get_func(int argc, VALUE *argv, VALUE obj, VALUE *ff, VALUE *xx);

static int get_func(int argc, VALUE *argv, VALUE obj, VALUE *ff, VALUE *xx)
{
   switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2) rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)",
          argc);
    CHECK_FUNCTION(argv[0]);
    *ff = argv[0];
    *xx = argv[1];
    break;
  default:
    if (argc != 1) rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)",
          argc);
    *ff = obj;
    *xx = argv[0];
    break;
  }
  return 0;
}

static VALUE rb_gsl_diff_eval(VALUE obj, VALUE xx,
            int (*diff)(const gsl_function *, double, double *, double *))
{
  gsl_function *f = NULL;
  double result, abserr;
  VALUE x, ary, aerr;
  gsl_vector *v = NULL, *vnew = NULL, *verr = NULL;
  gsl_matrix *m = NULL, *mnew = NULL, *merr = NULL;
  size_t n, i, j;
  int status;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2, *ptr3;
  struct NARRAY *na;
  VALUE ary2, ary3;
#endif
  Data_Get_Struct(obj, gsl_function, f);
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    status = (*diff)(f, NUM2DBL(xx), &result, &abserr);
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
      (*diff)(f, NUM2DBL(x), &result, &abserr);
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
        (*diff)(f, ptr1[i], &result, &abserr);  
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
        (*diff)(f, gsl_vector_get(v, i), &result, &abserr);  
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
          (*diff)(f, gsl_matrix_get(m, i, j), &result, &abserr);  
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
  return Qnil;  /* never reach here */
}

static VALUE rb_gsl_diff_central(int argc, VALUE *argv, VALUE obj)
{
  VALUE ff, xx;
  get_func(argc, argv, obj, &ff, &xx);
  return rb_gsl_diff_eval(ff, xx, gsl_diff_central);
}

static VALUE rb_gsl_diff_forward(int argc, VALUE *argv, VALUE obj)
{
  VALUE ff, xx;
  get_func(argc, argv, obj, &ff, &xx);
  return rb_gsl_diff_eval(ff, xx, gsl_diff_forward);
}

static VALUE rb_gsl_diff_backward(int argc, VALUE *argv, VALUE obj)
{
  VALUE ff, xx;
  get_func(argc, argv, obj, &ff, &xx);
  return rb_gsl_diff_eval(ff, xx, gsl_diff_backward);
}

void Init_gsl_diff(VALUE module)
{
  VALUE mgsl_diff;

  mgsl_diff = rb_define_module_under(module, "Diff");

  rb_define_method(cgsl_function, "diff_central", rb_gsl_diff_central, -1);
  rb_define_method(cgsl_function, "diff_forward", rb_gsl_diff_forward, -1);
  rb_define_method(cgsl_function, "diff_backward", rb_gsl_diff_backward, -1);

  rb_define_singleton_method(mgsl_diff, "central", rb_gsl_diff_central, -1);
  rb_define_singleton_method(mgsl_diff, "forward", rb_gsl_diff_forward, -1);
  rb_define_singleton_method(mgsl_diff, "backward", rb_gsl_diff_backward, -1);
}

