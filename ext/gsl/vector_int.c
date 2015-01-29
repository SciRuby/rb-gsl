/*
  vector_int.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada
        Modified by Seiya Nishizawa        14/Apr/2004

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_array.h"
#include "include/rb_gsl_complex.h"
#ifdef HAVE_NARRAY_H
#include "include/rb_gsl_with_narray.h"
#endif

VALUE rb_gsl_vector_int_inner_product(int argc, VALUE *argv, VALUE obj);
VALUE rb_gsl_vector_int_do_something(VALUE obj, void (*func)(gsl_vector_int*));

static VALUE rb_gsl_vector_int_to_i(VALUE obj)
{
  return obj;
}

VALUE rb_gsl_vector_int_to_f(VALUE obj)
{
  gsl_vector_int *v;
  gsl_vector *vnew;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_int, v);
  vnew = gsl_vector_alloc(v->size);
  for (i = 0; i < v->size; i++)
    gsl_vector_set(vnew, i, (double) gsl_vector_int_get(v, i));
  if (VECTOR_INT_COL_P(obj))
    return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, vnew);
  else
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
}

static VALUE rb_gsl_vector_int_to_complex(VALUE obj)
{
  gsl_vector_int *v;
  gsl_vector_complex *vnew;
  gsl_complex z;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_int, v);
  vnew = gsl_vector_complex_alloc(v->size);
  for (i = 0; i < v->size; i++) {
    GSL_SET_REAL(&z, (double) gsl_vector_int_get(v, i));
    GSL_SET_IMAG(&z, 0.0);
    gsl_vector_complex_set(vnew, i, z);
  }
  if (VECTOR_INT_COL_P(obj))
    return Data_Wrap_Struct(cgsl_vector_complex_col, 0, gsl_vector_complex_free, vnew);
  else
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
}

static VALUE rb_gsl_vector_int_coerce(VALUE obj, VALUE other)
{
  gsl_vector_int *v = NULL, *vnew = NULL;
  VALUE vv;
  Data_Get_Struct(obj, gsl_vector_int, v);
  switch (TYPE(other)) {
  case T_FIXNUM:
    vnew = gsl_vector_int_alloc(v->size);
    if (vnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_int_alloc failed");
    gsl_vector_int_set_all(vnew, FIX2INT(other));
    vv = Data_Wrap_Struct(VECTOR_INT_ROW_COL(obj), 0, gsl_vector_int_free, vnew);
    return rb_ary_new3(2, vv, obj);
    break;
  default:
    return rb_ary_new3(2, other, rb_gsl_vector_int_to_f(obj));
  }
  //  return rb_ary_new3(2, other, rb_gsl_vector_int_to_f(obj));
}

VALUE rb_gsl_vector_add(VALUE obj, VALUE b);
VALUE rb_gsl_vector_sub(VALUE obj, VALUE b);
VALUE rb_gsl_vector_mul(VALUE obj, VALUE b);
VALUE rb_gsl_vector_div(VALUE obj, VALUE b);
static VALUE rb_gsl_vector_int_add(VALUE obj, VALUE b)
{
  gsl_vector_int *v, *vnew, *vb;
  switch (TYPE(b)) {
  case T_FIXNUM:
    return rb_gsl_vector_int_add_constant(obj, b);
    break;
  case T_FLOAT:
    return rb_gsl_vector_add_constant(rb_gsl_vector_int_to_f(obj), b);
    break;
  default:
   if (rb_obj_is_kind_of(b, cgsl_vector_int)) {
     Data_Get_Struct(obj, gsl_vector_int, v);
     Data_Get_Struct(b, gsl_vector_int, vb);
     vnew = gsl_vector_int_alloc(v->size);
     gsl_vector_int_memcpy(vnew, v);
     gsl_vector_int_add(vnew, vb);
     return Data_Wrap_Struct(VECTOR_INT_ROW_COL(obj), 0, gsl_vector_int_free, vnew);
   } else {
     return rb_gsl_vector_add(rb_gsl_vector_int_to_f(obj), b);
   }
   break;
  }
}


static VALUE rb_gsl_vector_int_sub(VALUE obj, VALUE b)
{
  gsl_vector_int *v, *vnew, *vb;
  switch (TYPE(b)) {
  case T_FIXNUM:
    return rb_gsl_vector_int_add_constant(obj, INT2FIX(-FIX2INT(b)));
    break;
  case T_FLOAT:
    return rb_gsl_vector_add_constant(rb_gsl_vector_int_to_f(obj), rb_float_new(-NUM2DBL(b)));
    break;
  default:
    if (rb_obj_is_kind_of(b, cgsl_vector_int)) {
     Data_Get_Struct(obj, gsl_vector_int, v);
     Data_Get_Struct(b, gsl_vector_int, vb);
     vnew = gsl_vector_int_alloc(v->size);
     gsl_vector_int_memcpy(vnew, v);
     gsl_vector_int_sub(vnew, vb);
     return Data_Wrap_Struct(VECTOR_INT_ROW_COL(obj), 0, gsl_vector_int_free, vnew);
    } else {
      return rb_gsl_vector_sub(rb_gsl_vector_int_to_f(obj), b);
    }
    break;
  }
}

gsl_vector_int* mygsl_vector_int_mul_matrix(gsl_vector_int *v, gsl_matrix_int *m);

static VALUE rb_gsl_vector_int_mul(VALUE obj, VALUE b)
{
  VALUE argv[2];
  gsl_vector_int *v, *vnew, *v2;
  gsl_matrix_int *m;
  int val;
  size_t i, j;
  switch (TYPE(b)) {
  case T_FIXNUM:
  case T_FLOAT:
    return rb_gsl_vector_int_scale(obj, b);
    break;
  default:
    if (VECTOR_INT_ROW_P(obj) && VECTOR_INT_COL_P(b)) {
      argv[0] = obj;
      argv[1] = b;
      return rb_gsl_vector_int_inner_product(2, argv, CLASS_OF(obj));
    } else  if (VECTOR_INT_ROW_P(obj) && MATRIX_INT_P(b)) {
      Data_Get_Struct(obj, gsl_vector_int, v);
      Data_Get_Struct(b, gsl_matrix_int, m);
      vnew = mygsl_vector_int_mul_matrix(v, m);
      return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, vnew);
    } else if (VECTOR_INT_COL_P(obj) && VECTOR_INT_ROW_P(b)) {
      Data_Get_Struct(obj, gsl_vector_int, v);
      Data_Get_Struct(b, gsl_vector_int, v2);
      if (v->size != v2->size) rb_raise(rb_eIndexError, "Vector sizes does not match.");
      m = gsl_matrix_int_alloc(v->size, v2->size);
      for (i = 0; i < v->size; i++) {
  for (j = 0; j < v2->size; j++) {
    val = gsl_vector_int_get(v, i)*gsl_vector_int_get(v2, j);
    gsl_matrix_int_set(m, i, j, val);
  }
      }
      return Data_Wrap_Struct(cgsl_matrix_int, 0, gsl_matrix_int_free, m);
    } else {
      return rb_gsl_vector_mul(rb_gsl_vector_int_to_f(obj), b);
    }
    break;
  }
}

static VALUE rb_gsl_vector_int_div(VALUE obj, VALUE b)
{
  return rb_gsl_vector_div(rb_gsl_vector_int_to_f(obj), b);
}

void Init_gsl_vector_int(VALUE module)
{
  rb_define_method(cgsl_vector_int, "to_f", rb_gsl_vector_int_to_f, 0);
  rb_define_method(cgsl_vector_int, "to_i", rb_gsl_vector_int_to_i, 0);
  rb_define_method(cgsl_vector_int, "to_complex", rb_gsl_vector_int_to_complex, 0);

  /*****/
  rb_define_method(cgsl_vector_int, "coerce", rb_gsl_vector_int_coerce, 1);

  rb_define_method(cgsl_vector_int, "add", rb_gsl_vector_int_add, 1);
  rb_define_method(cgsl_vector_int, "sub", rb_gsl_vector_int_sub, 1);
  rb_define_method(cgsl_vector_int, "mul", rb_gsl_vector_int_mul, 1);
  rb_define_method(cgsl_vector_int, "div", rb_gsl_vector_int_div, 1);

  rb_define_alias(cgsl_vector_int, "+", "add");
  rb_define_alias(cgsl_vector_int, "-", "sub");
  rb_define_alias(cgsl_vector_int, "*", "mul");
  rb_define_alias(cgsl_vector_int, "/", "div");

  Init_gsl_vector_int_init(module);
}

