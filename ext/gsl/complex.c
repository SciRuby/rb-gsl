/*
  complex.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "include/rb_gsl_complex.h"
#include "include/rb_gsl_array.h"

VALUE cgsl_complex;

// Set real/imag components of gsl_complex from Ruby object.
// If z is NULL, it is as if it pointed to a temporary (0+0i) value,
// otherwise the gsl_complex pointed to z is modified.
// Returns the resulting gsl_complex value in both cases.
//
// NOTE: This function does not always set both components of *z, so if *z is
// non-null, it should be initialized before calling this function.
gsl_complex rb_gsl_obj_to_gsl_complex(VALUE obj, gsl_complex *z)
{
  VALUE vre, vim;
  gsl_complex tmp, *zz;

  if(!z) {
    z = &tmp;
    GSL_SET_COMPLEX(z, 0.0, 0.0);
  }

  if(obj == Qnil) {
    return *z;
  }

  switch(TYPE(obj)) {
  case T_ARRAY:
    vre = rb_ary_entry(obj,0);
    vim = rb_ary_entry(obj,1);
    if (!NIL_P(vre)) GSL_SET_REAL(z, NUM2DBL(vre));
    if (!NIL_P(vim)) GSL_SET_IMAG(z, NUM2DBL(vim));
    break;
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    *z = gsl_complex_rect(NUM2DBL(obj), 0.0);
    break;
  default:
    if (rb_obj_is_kind_of(obj, cgsl_complex)) {
      Data_Get_Struct(obj, gsl_complex, zz);
      *z = *zz;
    } else {
      rb_raise(rb_eTypeError,
          "wrong type %s, (nil, Array, Float, Integer, or GSL::Complex expected)",
          rb_class2name(CLASS_OF(obj)));
    }
    break;
  }
  return *z;
}

static VALUE rb_gsl_complex_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_complex *c = NULL;
  VALUE obj;
  obj = Data_Make_Struct(klass, gsl_complex, 0, free, c);
  switch (argc) {
  case 1:
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      *c = ary2complex(argv[0]);
      break;
    case T_FLOAT:
    case T_FIXNUM:
    case T_BIGNUM:
      Need_Float(argv[0]);
      *c = gsl_complex_rect(NUM2DBL(argv[0]), 0.0);
      break;
    default:
      rb_raise(rb_eTypeError, "wrong argument type %s",
         rb_class2name(CLASS_OF(argv[0])));
    }
    break;
  case 2:
    Need_Float(argv[0]); Need_Float(argv[1]);
    *c = gsl_complex_rect(NUM2DBL(argv[0]),  NUM2DBL(argv[1]));
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  }
  return obj;
}

gsl_complex* make_complex(double re, double im)
{
  gsl_complex *c = NULL;
  c = ALLOC(gsl_complex);
  *c = gsl_complex_rect(re, im);
  return c;
}

static VALUE rb_gsl_complex_polar(VALUE klass, VALUE r, VALUE theta)
{
  VALUE obj;
  gsl_complex *c = NULL;
  Need_Float(r); Need_Float(theta);
  obj = Data_Make_Struct(klass, gsl_complex, 0, free, c);
  *c = gsl_complex_polar(NUM2DBL(r), NUM2DBL(theta));
  return obj;
}

static VALUE rb_gsl_complex_get(VALUE obj, VALUE ii)
{
  gsl_complex *c = NULL;
  int i;
  CHECK_FIXNUM(ii);
  Data_Get_Struct(obj, gsl_complex, c);
  i = FIX2INT(ii);
  switch (i) {
  case 0:
    return  rb_float_new(GSL_REAL(*c));
    break;
  case 1:
    return  rb_float_new(GSL_IMAG(*c));
    break;
  default:
    rb_raise(rb_eArgError, "wrong argument (%d for 0 or 1)", i);
    break;
  }
  return Qnil;
}

static VALUE rb_gsl_complex_real(VALUE obj)
{
  gsl_complex *c = NULL;
  Data_Get_Struct(obj, gsl_complex, c);
  return  rb_float_new(GSL_REAL(*c));
}

static VALUE rb_gsl_complex_imag(VALUE obj)
{
  gsl_complex *c = NULL;
  Data_Get_Struct(obj, gsl_complex, c);
  return  rb_float_new(GSL_IMAG(*c));
}

static VALUE rb_gsl_complex_print(VALUE obj)
{
  gsl_complex *c = NULL;
  Data_Get_Struct(obj, gsl_complex, c);
  fprintf(stdout, "[%4.3e %4.3e] \n", GSL_REAL(*c), GSL_IMAG(*c));
  return obj;
}

static VALUE rb_gsl_complex_printf(VALUE obj, VALUE s)
{
  gsl_complex *c = NULL;
  char tmp[32], format[64];
  Check_Type(s, T_STRING);
  Data_Get_Struct(obj, gsl_complex, c);
  strcpy(tmp, STR2CSTR(s));
  sprintf(format, "%s %s\n", tmp, tmp);
  fprintf(stdout, format, GSL_REAL(*c), GSL_IMAG(*c));
  return obj;
}

static VALUE rb_gsl_complex_return_double(double (*func)(gsl_complex), VALUE obj);
static VALUE rb_gsl_complex_return_double(double (*func)(gsl_complex), VALUE obj)
{
  gsl_complex *c = NULL;
  Data_Get_Struct(obj, gsl_complex, c);
  return  rb_float_new((*func)(*c));
}

static VALUE rb_gsl_complex_arg(VALUE obj)
{
  return rb_gsl_complex_return_double(gsl_complex_arg, obj);
}

static VALUE rb_gsl_complex_abs(VALUE obj)
{
  return rb_gsl_complex_return_double(gsl_complex_abs, obj);
}

static VALUE rb_gsl_complex_abs2(VALUE obj)
{
  return rb_gsl_complex_return_double(gsl_complex_abs2, obj);
}

static VALUE rb_gsl_complex_logabs(VALUE obj)
{
  return rb_gsl_complex_return_double(gsl_complex_logabs, obj);
}

static VALUE rb_gsl_complex_arithmetics2(gsl_complex (*func)(gsl_complex, double),
           VALUE obj, VALUE xx);

static VALUE rb_gsl_complex_arithmetics2(gsl_complex (*func)(gsl_complex, double),
          VALUE obj, VALUE xx)
{
  gsl_complex *a = NULL, *c = NULL, tmp;
  VALUE obj2;
  double x;
  Need_Float(xx);
  Data_Get_Struct(obj, gsl_complex, a);
  x = NUM2DBL(xx);
  tmp = (*func)(*a, x);
  obj2 = Data_Make_Struct(cgsl_complex, gsl_complex, 0, free, c);
  *c = tmp;
  return obj2;
}

static VALUE rb_gsl_complex_add_real(VALUE obj, VALUE xx)
{
  return rb_gsl_complex_arithmetics2(gsl_complex_add_real, obj, xx);
}

static VALUE rb_gsl_complex_sub_real(VALUE obj, VALUE xx)
{
  return rb_gsl_complex_arithmetics2(gsl_complex_sub_real, obj, xx);
}

static VALUE rb_gsl_complex_mul_real(VALUE obj, VALUE xx)
{
  return rb_gsl_complex_arithmetics2(gsl_complex_mul_real, obj, xx);
}

static VALUE rb_gsl_complex_div_real(VALUE obj, VALUE xx)
{
  return rb_gsl_complex_arithmetics2(gsl_complex_div_real, obj, xx);
}

static VALUE rb_gsl_complex_add_imag(VALUE obj, VALUE xx)
{
  return rb_gsl_complex_arithmetics2(gsl_complex_add_imag, obj, xx);
}

static VALUE rb_gsl_complex_sub_imag(VALUE obj, VALUE xx)
{
  return rb_gsl_complex_arithmetics2(gsl_complex_sub_imag, obj, xx);
}

static VALUE rb_gsl_complex_mul_imag(VALUE obj, VALUE xx)
{
  return rb_gsl_complex_arithmetics2(gsl_complex_mul_imag, obj, xx);
}

static VALUE rb_gsl_complex_div_imag(VALUE obj, VALUE xx)
{
  return rb_gsl_complex_arithmetics2(gsl_complex_div_imag, obj, xx);
}

static VALUE rb_gsl_complex_operate(gsl_complex (*func)(gsl_complex), VALUE obj);
static VALUE rb_gsl_complex_operate2(gsl_complex (*func)(gsl_complex), int argc, VALUE *argv, VALUE obj);

static VALUE rb_gsl_complex_operate(gsl_complex (*func)(gsl_complex), VALUE obj)
{
  gsl_complex *c = NULL, *cnew = NULL;
  Data_Get_Struct(obj, gsl_complex, c);
  cnew = ALLOC(gsl_complex);
  *cnew = (*func)(*c);
  return Data_Wrap_Struct(cgsl_complex, 0, free, cnew);
}

static VALUE rb_gsl_complex_operate2(gsl_complex (*func)(gsl_complex), int argc, VALUE *argv, VALUE obj)
{
  gsl_complex *c = NULL, *cnew = NULL, tmp;
  gsl_vector_complex *v, *vnew;
  gsl_matrix_complex *m, *mnew;
  size_t i, j;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    switch (argc) {
    case 1:
      switch (TYPE(argv[0])) {
      case T_ARRAY:
        tmp = ary2complex(argv[0]);
        c = &tmp;
        break;
      default:
        if (VECTOR_COMPLEX_P(argv[0])) {
          Data_Get_Struct(argv[0], gsl_vector_complex, v);
          vnew = gsl_vector_complex_alloc(v->size);
          for (i = 0; i < v->size; i++) {
            c = GSL_COMPLEX_AT(v, i);
            gsl_vector_complex_set(vnew, i, (*func)(*c));
          }
          return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_free, vnew);
        } else if (MATRIX_COMPLEX_P(obj)) {
          Data_Get_Struct(obj, gsl_matrix_complex, m);
          mnew = gsl_matrix_complex_alloc(m->size1, m->size2);
          for (i = 0; i < m->size1; i++) {
            for (j = 0; j < m->size2; j++) {
              gsl_matrix_complex_set(mnew, i, j,
                   (*func)(gsl_matrix_complex_get(m, i, j)));
            }
          }
          return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, mnew);
        } else {
          CHECK_COMPLEX(argv[0]);
          Data_Get_Struct(argv[0], gsl_complex, c);
        }
        break;
      }
      break;
    case 2:
      c = &tmp;
      GSL_SET_REAL(c, NUM2DBL(argv[0]));
      GSL_SET_IMAG(c, NUM2DBL(argv[1]));
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
      break;
    }
    break;
  default:
    Data_Get_Struct(obj, gsl_complex, c);
    break;
  }
  cnew = ALLOC(gsl_complex);
  *cnew = (*func)(*c);
  return Data_Wrap_Struct(cgsl_complex, 0, free, cnew);
}

static VALUE rb_gsl_complex_conjugate(VALUE obj)
{
  return rb_gsl_complex_operate(gsl_complex_conjugate, obj);
}

static VALUE rb_gsl_complex_inverse(VALUE obj)
{
  return rb_gsl_complex_operate(gsl_complex_inverse, obj);
}

static VALUE rb_gsl_complex_negative(VALUE obj)
{
  return rb_gsl_complex_operate(gsl_complex_negative, obj);
}

/* singleton */
static VALUE rb_gsl_complex_sqrt_real(VALUE obj, VALUE x)
{
  gsl_complex *cnew = NULL, tmp;
  Need_Float(x);
  tmp = gsl_complex_sqrt_real(NUM2DBL(x));
  cnew = ALLOC(gsl_complex);
  *cnew = tmp;
  return Data_Wrap_Struct(cgsl_complex, 0, free, cnew);
}

static VALUE rb_gsl_complex_sqrt(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex *z, tmp;
  switch (TYPE(obj)) {
  case T_MODULE: case T_CLASS: case T_OBJECT:
    switch (argc) {
    case 1:
      switch (TYPE(argv[0])) {
      case T_FIXNUM: case T_FLOAT:
        return rb_gsl_complex_sqrt_real(obj, argv[0]);
        break;
      case T_ARRAY:
        tmp = ary2complex(argv[0]);
        z = ALLOC(gsl_complex);
        *z = gsl_complex_sqrt(tmp);
        return Data_Wrap_Struct(cgsl_complex, 0, free, z);
        break;
      default:
        CHECK_COMPLEX(argv[0]);
        Data_Get_Struct(argv[0], gsl_complex, z);
        tmp = *z;
        z = ALLOC(gsl_complex);
        *z = gsl_complex_sqrt(tmp);
        return Data_Wrap_Struct(cgsl_complex, 0, free, z);
        break;
      }
      break;
    case 2:
      z = ALLOC(gsl_complex);
      GSL_SET_REAL(&tmp, NUM2DBL(argv[0]));
      GSL_SET_IMAG(&tmp, NUM2DBL(argv[1]));
      *z = gsl_complex_sqrt(tmp);
      return Data_Wrap_Struct(cgsl_complex, 0, free, z);
      break;
    default:
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
      break;
    }
    break;
  default:
    return rb_gsl_complex_operate2(gsl_complex_sqrt, argc, argv, obj);
    break;
  }
}

VALUE rb_gsl_complex_pow(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex *c = NULL, *a = NULL, *cnew = NULL, tmpc, tmpa;
  gsl_vector_complex *v, *vnew;
  gsl_matrix_complex *m, *mnew;
  size_t i, j;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2)
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)", argc);
    switch (TYPE(argv[1])) {
    case T_ARRAY:
      tmpa = ary2complex(argv[1]);
      a = &tmpa;
      break;
    case T_FLOAT: case T_FIXNUM:
      return rb_gsl_complex_pow_real(argc, argv, obj);
      break;
    default:
      CHECK_COMPLEX(argv[1]);
      Data_Get_Struct(argv[1], gsl_complex, a);
      break;
    }

    switch (TYPE(argv[0])) {
    case T_ARRAY:
      tmpc = ary2complex(argv[0]);
      c = &tmpc;
      break;
    default:
      if (VECTOR_COMPLEX_P(argv[0])) {
        Data_Get_Struct(argv[0], gsl_vector_complex, v);
        vnew = gsl_vector_complex_alloc(v->size);
        for (i = 0; i < v->size; i++) {
          c = GSL_COMPLEX_AT(v, i);
          gsl_vector_complex_set(vnew, i, gsl_complex_pow(*c, *a));
        }
        return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
      }
      if (MATRIX_COMPLEX_P(argv[0])) {
        Data_Get_Struct(argv[0], gsl_matrix_complex, m);
        mnew = gsl_matrix_complex_alloc(m->size1, m->size2);
        for (i = 0; i < m->size1; i++) {
          for (j = 0; j < m->size2; j++) {
            c = gsl_matrix_complex_ptr(m, i, j);
            gsl_matrix_complex_set(mnew, i, j, gsl_complex_pow(*c, *a));
          }
        }
        return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, mnew);
      }
      CHECK_COMPLEX(argv[0]);
      Data_Get_Struct(argv[0], gsl_complex, c);
      break;
    }

    break;
  default:
    if (argc != 1)
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
    CHECK_COMPLEX(obj);
    Data_Get_Struct(obj, gsl_complex, c);
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      tmpa = ary2complex(argv[0]);
      a = &tmpa;
      break;
    case T_FLOAT: case T_FIXNUM:
      return rb_gsl_complex_pow_real(argc, argv, obj);
      break;
    default:
      CHECK_COMPLEX(argv[0]);
      Data_Get_Struct(argv[0], gsl_complex, a);
      break;
    }
    break;
  }
  cnew = ALLOC(gsl_complex);
  *cnew = gsl_complex_pow(*c, *a);
  return Data_Wrap_Struct(cgsl_complex, 0, free, cnew);
}

VALUE rb_gsl_complex_pow_real(int argc, VALUE *argv, VALUE obj)
{
  double a = 1;
  gsl_complex *c = NULL, *cnew = NULL, tmpc;
  gsl_vector_complex *v, *vnew;
  gsl_matrix_complex *m, *mnew;
  size_t i, j;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2)
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)", argc);
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      tmpc = ary2complex(argv[0]);
      c = &tmpc;
      a = NUM2DBL(argv[1]);
      break;
    default:
      if (VECTOR_COMPLEX_P(argv[0])) {
  Data_Get_Struct(argv[0], gsl_vector_complex, v);
  vnew = gsl_vector_complex_alloc(v->size);
  a = NUM2DBL(argv[1]);
  for (i = 0; i < v->size; i++) {
    c = GSL_COMPLEX_AT(v, i);
    tmpc = gsl_complex_pow_real(*c, a);
    gsl_vector_complex_set(vnew, i, tmpc);
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_free, vnew);
      }
      if (MATRIX_COMPLEX_P(argv[0])) {
  Data_Get_Struct(argv[0], gsl_matrix_complex, m);
  mnew = gsl_matrix_complex_alloc(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      tmpc = gsl_complex_pow_real(gsl_matrix_complex_get(m, i, j), a);
      gsl_matrix_complex_set(mnew, i, j, tmpc);
    }
  }
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, mnew);
      }
      CHECK_COMPLEX(argv[0]);
      Data_Get_Struct(argv[0], gsl_complex, c);
      break;
    }
    Need_Float(argv[1]);
    a = NUM2DBL(argv[1]);
    break;
  default:
    if (argc != 1)
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
    CHECK_COMPLEX(obj);
    Need_Float(argv[0]);
    Data_Get_Struct(obj, gsl_complex, c);
    a = NUM2DBL(argv[0]);
    break;
  }
  cnew = ALLOC(gsl_complex);
  *cnew = gsl_complex_pow_real(*c, a);
  return Data_Wrap_Struct(cgsl_complex, 0, free, cnew);
}

static VALUE rb_gsl_complex_exp(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_exp, argc, argv, obj);
}

static VALUE rb_gsl_complex_log(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_log, argc, argv, obj);
}

static VALUE rb_gsl_complex_log10(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_log10, argc, argv, obj);
}

static VALUE rb_gsl_complex_log_b(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex *c = NULL, *a = NULL, *cnew = NULL, tmpc, tmpa;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 2)
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 2)", argc);
    switch (TYPE(argv[1])) {
    case T_ARRAY:
      tmpa = ary2complex(argv[1]);
      a = &tmpa;
      break;
    default:
      CHECK_COMPLEX(argv[1]);
      Data_Get_Struct(argv[1], gsl_complex, a);
      break;
    }

    switch (TYPE(argv[0])) {
    case T_ARRAY:
      tmpc = ary2complex(argv[0]);
      c = &tmpc;
      break;
    default:
      CHECK_COMPLEX(argv[0]);
      Data_Get_Struct(argv[0], gsl_complex, c);
      break;
    }

    break;
  default:
    if (argc != 1)
      rb_raise(rb_eArgError, "wrong number of arguments (%d for 1)", argc);
    CHECK_COMPLEX(obj);
    Data_Get_Struct(obj, gsl_complex, c);
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      tmpa = ary2complex(argv[0]);
      a = &tmpa;
      break;
    default:
      CHECK_COMPLEX(argv[0]);
      Data_Get_Struct(argv[0], gsl_complex, a);
      break;
    }
    break;
  }
  cnew = ALLOC(gsl_complex);
  *cnew = gsl_complex_log_b(*c, *a);
  return Data_Wrap_Struct(cgsl_complex, 0, free, cnew);
}

static VALUE rb_gsl_complex_sin(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_sin, argc, argv, obj);
}

static VALUE rb_gsl_complex_cos(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_cos, argc, argv, obj);
}

static VALUE rb_gsl_complex_tan(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_tan, argc, argv, obj);
}

static VALUE rb_gsl_complex_sec(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_sec, argc, argv, obj);
}

static VALUE rb_gsl_complex_csc(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_csc, argc, argv, obj);
}

static VALUE rb_gsl_complex_cot(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_cot, argc, argv, obj);
}

/*****/
static VALUE rb_gsl_complex_arcsin(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_arcsin, argc, argv, obj);
}

static VALUE rb_gsl_complex_arcsin_real(VALUE obj, VALUE xx)
{
  gsl_complex tmp, *cnew = NULL;
  double x;
  Need_Float(xx);
  x = NUM2DBL(xx);
  tmp = gsl_complex_arcsin_real(x);
  cnew = ALLOC(gsl_complex);
  *cnew = tmp;
  return Data_Wrap_Struct(cgsl_complex, 0, free, cnew);
}

static VALUE rb_gsl_complex_arccos(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_arccos, argc, argv, obj);
}

static VALUE rb_gsl_complex_arccos_real(VALUE obj, VALUE xx)
{
  gsl_complex tmp, *cnew = NULL;
  double x;
  Need_Float(xx);
  x = NUM2DBL(xx);
  tmp = gsl_complex_arccos_real(x);
  cnew = ALLOC(gsl_complex);
  *cnew = tmp;
  return Data_Wrap_Struct(cgsl_complex, 0, free, cnew);
}

static VALUE rb_gsl_complex_arctan(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_arctan, argc, argv, obj);
}

static VALUE rb_gsl_complex_arcsec(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_arcsec, argc, argv, obj);
}

static VALUE rb_gsl_complex_arcsec_real(VALUE obj, VALUE xx)
{
  gsl_complex tmp, *cnew = NULL;
  double x;
  Need_Float(xx);
  x = NUM2DBL(xx);
  tmp = gsl_complex_arcsec_real(x);
  cnew = ALLOC(gsl_complex);
  *cnew = tmp;
  return Data_Wrap_Struct(cgsl_complex, 0, free, cnew);
}

static VALUE rb_gsl_complex_arccsc(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_arccsc, argc, argv, obj);
}

static VALUE rb_gsl_complex_arccsc_real(VALUE obj, VALUE xx)
{
  gsl_complex tmp, *cnew = NULL;
  double x;
  Need_Float(xx);
  x = NUM2DBL(xx);
  tmp = gsl_complex_arccsc_real(x);
  cnew = ALLOC(gsl_complex);
  *cnew = tmp;
  return Data_Wrap_Struct(cgsl_complex, 0, free, cnew);
}

static VALUE rb_gsl_complex_arccot(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_arccot, argc, argv, obj);
}

static VALUE rb_gsl_complex_sinh(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_sinh, argc, argv, obj);
}

static VALUE rb_gsl_complex_cosh(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_cosh, argc, argv, obj);
}

static VALUE rb_gsl_complex_tanh(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_tanh, argc, argv, obj);
}

static VALUE rb_gsl_complex_sech(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_sech, argc, argv, obj);
}

static VALUE rb_gsl_complex_csch(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_csch, argc, argv, obj);
}

static VALUE rb_gsl_complex_coth(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_coth, argc, argv, obj);
}

static VALUE rb_gsl_complex_arcsinh(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_arcsinh, argc, argv, obj);
}

static VALUE rb_gsl_complex_arccosh_real(VALUE obj, VALUE xx)
{
  gsl_complex tmp, *cnew = NULL;
  double x;
  Need_Float(xx);
  x = NUM2DBL(xx);
  tmp = gsl_complex_arccosh_real(x);
  cnew = ALLOC(gsl_complex);
  *cnew = tmp;
  return Data_Wrap_Struct(cgsl_complex, 0, free, cnew);
}

static VALUE rb_gsl_complex_arccosh(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_arccosh, argc, argv, obj);
}

static VALUE rb_gsl_complex_arctanh(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_arctanh, argc, argv, obj);
}

static VALUE rb_gsl_complex_arctanh_real(VALUE obj, VALUE xx)
{
  gsl_complex tmp, *cnew = NULL;
  double x;
  Need_Float(xx);
  x = NUM2DBL(xx);
  tmp = gsl_complex_arctanh_real(x);
  cnew = ALLOC(gsl_complex);
  *cnew = tmp;
  return Data_Wrap_Struct(cgsl_complex, 0, free, cnew);
}

static VALUE rb_gsl_complex_arcsech(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_arcsech, argc, argv, obj);
}

static VALUE rb_gsl_complex_arccsch(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_arccsch, argc, argv, obj);
}

static VALUE rb_gsl_complex_arccoth(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_complex_operate2(gsl_complex_arccoth, argc, argv, obj);
}

int rbgsl_complex_equal(const gsl_complex *z1, const gsl_complex *z2,
      double eps)
{
  if (gsl_fcmp(z1->dat[0], z2->dat[0], eps) != 0) return 0;
  if (gsl_fcmp(z1->dat[1], z2->dat[1], eps) != 0) return 0;
  return 1;
}

int rbgsl_complex_zero(const gsl_complex *z1)
{
  if (z1->dat[0] != 0.0) return 0;
  if (z1->dat[1] != 0.0) return 0;
  return 1;
}

static VALUE rb_gsl_complex_equal(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex *z1 = NULL, *z2 = NULL;
  double eps = 1e-8;

  switch (argc) {
  case 1:
    CHECK_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_complex, z2);
    break;
  case 2:
    CHECK_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_complex, z2);
    eps = NUM2DBL(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of argumsnts (%d for 1 or 2)", argc);
  }
  Data_Get_Struct(obj, gsl_complex, z1);
  if (rbgsl_complex_equal(z1, z2, eps)) return Qtrue;
  else return Qfalse;
}

static VALUE rb_gsl_complex_zero(VALUE obj)
{
  gsl_complex *z = NULL;
  Data_Get_Struct(obj, gsl_complex, z);
  if (rbgsl_complex_zero(z)) return Qtrue;
  else return Qfalse;
}

static VALUE rb_gsl_complex_to_s(VALUE obj)
{
  char buf[256];
  gsl_complex *z;
  Data_Get_Struct(obj, gsl_complex, z);
  sprintf(buf, "[ %4.3e %4.3e ]", GSL_REAL(*z), GSL_IMAG(*z));
  return rb_str_new2(buf);
}

static VALUE rb_gsl_complex_inspect(VALUE obj)
{
  char buf[256];
  VALUE str;
  sprintf(buf, "%s\n", rb_class2name(CLASS_OF(obj)));
  str = rb_str_new2(buf);
  return rb_str_concat(str, rb_gsl_complex_to_s(obj));
}

void Init_gsl_complex(VALUE module)
{
  cgsl_complex = rb_define_class_under(module, "Complex",  rb_cNumeric);
  rb_define_singleton_method(cgsl_complex, "alloc", rb_gsl_complex_new, -1);
  rb_define_singleton_method(cgsl_complex, "rect", rb_gsl_complex_new, -1);
  rb_define_singleton_method(cgsl_complex, "[]", rb_gsl_complex_new, -1);
  rb_define_singleton_method(cgsl_complex, "polar", rb_gsl_complex_polar, 2);

  rb_define_method(cgsl_complex, "equal?", rb_gsl_complex_equal, -1);
  rb_define_alias(cgsl_complex, "==", "equal?");
  rb_define_method(cgsl_complex, "zero?", rb_gsl_complex_zero, 0);

  rb_define_method(cgsl_complex, "get", rb_gsl_complex_get, 1);
  rb_define_alias(cgsl_complex, "[]", "get");
  rb_define_method(cgsl_complex, "real", rb_gsl_complex_real, 0);
  rb_define_alias(cgsl_complex, "re", "real");
  rb_define_alias(cgsl_complex, "REAL", "real");
  rb_define_method(cgsl_complex, "imag", rb_gsl_complex_imag, 0);
  rb_define_alias(cgsl_complex, "im", "imag");
  rb_define_alias(cgsl_complex, "IMAG", "imag");
  rb_define_method(cgsl_complex, "print", rb_gsl_complex_print, 0);
  rb_define_method(cgsl_complex, "printf", rb_gsl_complex_printf, 1);
  rb_define_alias(cgsl_complex, "show", "print");

  rb_define_method(cgsl_complex, "arg", rb_gsl_complex_arg, 0);
  rb_define_alias(cgsl_complex, "angle", "arg");
  rb_define_alias(cgsl_complex, "phase", "arg");
  rb_define_method(cgsl_complex, "abs", rb_gsl_complex_abs, 0);
  rb_define_alias(cgsl_complex, "amp", "abs");
  rb_define_alias(cgsl_complex, "mag", "abs");
  rb_define_method(cgsl_complex, "abs2", rb_gsl_complex_abs2, 0);
  rb_define_method(cgsl_complex, "logabs", rb_gsl_complex_logabs, 0);

  rb_define_method(cgsl_complex, "add_real", rb_gsl_complex_add_real, 1);
  rb_define_method(cgsl_complex, "sub_real", rb_gsl_complex_sub_real, 1);
  rb_define_method(cgsl_complex, "mul_real", rb_gsl_complex_mul_real, 1);
  rb_define_method(cgsl_complex, "div_real", rb_gsl_complex_div_real, 1);
  rb_define_method(cgsl_complex, "add_imag", rb_gsl_complex_add_imag, 1);
  rb_define_method(cgsl_complex, "sub_imag", rb_gsl_complex_sub_imag, 1);
  rb_define_method(cgsl_complex, "mul_imag", rb_gsl_complex_mul_imag, 1);
  rb_define_method(cgsl_complex, "div_imag", rb_gsl_complex_div_imag, 1);

  rb_define_method(cgsl_complex, "conjugate", rb_gsl_complex_conjugate, 0);
  rb_define_alias(cgsl_complex, "conj", "conjugate");
  rb_define_method(cgsl_complex, "inverse", rb_gsl_complex_inverse, 0);
  rb_define_method(cgsl_complex, "negative", rb_gsl_complex_negative, 0);

  rb_define_singleton_method(cgsl_complex, "sqrt", rb_gsl_complex_sqrt, -1);
  rb_define_singleton_method(cgsl_complex, "sqrt_real", rb_gsl_complex_sqrt_real, 1);
  rb_define_singleton_method(cgsl_complex, "pow", rb_gsl_complex_pow, -1);
  rb_define_singleton_method(cgsl_complex, "pow_real", rb_gsl_complex_pow_real, -1);

  rb_define_singleton_method(cgsl_complex, "exp", rb_gsl_complex_exp, -1);
  rb_define_singleton_method(cgsl_complex, "log", rb_gsl_complex_log, -1);
  rb_define_singleton_method(cgsl_complex, "log10", rb_gsl_complex_log10, -1);
  rb_define_singleton_method(cgsl_complex, "log_b", rb_gsl_complex_log_b, -1);

  rb_define_singleton_method(cgsl_complex, "sin", rb_gsl_complex_sin, -1);
  rb_define_singleton_method(cgsl_complex, "cos", rb_gsl_complex_cos, -1);
  rb_define_singleton_method(cgsl_complex, "tan", rb_gsl_complex_tan, -1);
  rb_define_singleton_method(cgsl_complex, "sec", rb_gsl_complex_sec, -1);
  rb_define_singleton_method(cgsl_complex, "csc", rb_gsl_complex_csc, -1);
  rb_define_singleton_method(cgsl_complex, "cot", rb_gsl_complex_cot, -1);

  rb_define_singleton_method(cgsl_complex, "arcsin", rb_gsl_complex_arcsin, -1);
  rb_define_singleton_method(cgsl_complex, "arcsin_real", rb_gsl_complex_arcsin_real, 1);
  rb_define_singleton_method(cgsl_complex, "arccos", rb_gsl_complex_arccos, -1);
  rb_define_singleton_method(cgsl_complex, "arccos_real", rb_gsl_complex_arccos_real, 1);
  rb_define_singleton_method(cgsl_complex, "arctan", rb_gsl_complex_arctan, -1);
  rb_define_singleton_method(cgsl_complex, "arcsec", rb_gsl_complex_arcsec, -1);
  rb_define_singleton_method(cgsl_complex, "arcsec_real", rb_gsl_complex_arcsec_real, 1);
  rb_define_singleton_method(cgsl_complex, "arccsc", rb_gsl_complex_arccsc, -1);
  rb_define_singleton_method(cgsl_complex, "arccsc_real", rb_gsl_complex_arccsc_real, 1);
  rb_define_singleton_method(cgsl_complex, "arccot", rb_gsl_complex_arccot, -1);

  rb_define_singleton_method(cgsl_complex, "sinh", rb_gsl_complex_sinh, -1);
  rb_define_singleton_method(cgsl_complex, "cosh", rb_gsl_complex_cosh, -1);
  rb_define_singleton_method(cgsl_complex, "tanh", rb_gsl_complex_tanh, -1);
  rb_define_singleton_method(cgsl_complex, "sech", rb_gsl_complex_sech, -1);
  rb_define_singleton_method(cgsl_complex, "csch", rb_gsl_complex_csch, -1);
  rb_define_singleton_method(cgsl_complex, "coth", rb_gsl_complex_coth, -1);

  rb_define_singleton_method(cgsl_complex, "arcsinh", rb_gsl_complex_arcsinh, -1);
  rb_define_singleton_method(cgsl_complex, "arccosh", rb_gsl_complex_arccosh, -1);
  rb_define_singleton_method(cgsl_complex, "arccosh_real", rb_gsl_complex_arccosh_real, 1);
  rb_define_singleton_method(cgsl_complex, "arctanh", rb_gsl_complex_arctanh, -1);
  rb_define_singleton_method(cgsl_complex, "arctanh_real", rb_gsl_complex_arctanh_real, 1);
  rb_define_singleton_method(cgsl_complex, "arcsech", rb_gsl_complex_arcsech, -1);
  rb_define_singleton_method(cgsl_complex, "arccsch", rb_gsl_complex_arccsch, -1);
  rb_define_singleton_method(cgsl_complex, "arccoth", rb_gsl_complex_arccoth, -1);

  /***/

  rb_define_method(cgsl_complex, "pow", rb_gsl_complex_pow, -1);
  rb_define_method(cgsl_complex, "pow_real", rb_gsl_complex_pow_real, -1);

  rb_define_method(cgsl_complex, "exp", rb_gsl_complex_exp, -1);
  rb_define_method(cgsl_complex, "log", rb_gsl_complex_log, -1);
  rb_define_method(cgsl_complex, "log10", rb_gsl_complex_log10, -1);
  rb_define_method(cgsl_complex, "log_b", rb_gsl_complex_log_b, -1);

  rb_define_method(cgsl_complex, "sin", rb_gsl_complex_sin, -1);
  rb_define_method(cgsl_complex, "cos", rb_gsl_complex_cos, -1);
  rb_define_method(cgsl_complex, "tan", rb_gsl_complex_tan, -1);
  rb_define_method(cgsl_complex, "sec", rb_gsl_complex_sec, -1);
  rb_define_method(cgsl_complex, "csc", rb_gsl_complex_csc, -1);
  rb_define_method(cgsl_complex, "cot", rb_gsl_complex_cot, -1);

  rb_define_method(cgsl_complex, "arcsin", rb_gsl_complex_arcsin, -1);
  rb_define_method(cgsl_complex, "arccos", rb_gsl_complex_arccos, -1);
  rb_define_method(cgsl_complex, "arctan", rb_gsl_complex_arctan, -1);
  rb_define_method(cgsl_complex, "arcsec", rb_gsl_complex_arcsec, -1);
  rb_define_method(cgsl_complex, "arccsc", rb_gsl_complex_arccsc, -1);
  rb_define_method(cgsl_complex, "arccot", rb_gsl_complex_arccot, -1);

  rb_define_method(cgsl_complex, "sinh", rb_gsl_complex_sinh, -1);
  rb_define_method(cgsl_complex, "cosh", rb_gsl_complex_cosh, -1);
  rb_define_method(cgsl_complex, "tanh", rb_gsl_complex_tanh, -1);
  rb_define_method(cgsl_complex, "sech", rb_gsl_complex_sech, -1);
  rb_define_method(cgsl_complex, "csch", rb_gsl_complex_csch, -1);
  rb_define_method(cgsl_complex, "coth", rb_gsl_complex_coth, -1);

  rb_define_method(cgsl_complex, "arcsinh", rb_gsl_complex_arcsinh, -1);
  rb_define_method(cgsl_complex, "arccosh", rb_gsl_complex_arccosh, -1);
  rb_define_method(cgsl_complex, "arctanh", rb_gsl_complex_arctanh, -1);
  rb_define_method(cgsl_complex, "arcsech", rb_gsl_complex_arcsech, -1);
  rb_define_method(cgsl_complex, "arccsch", rb_gsl_complex_arccsch, -1);
  rb_define_method(cgsl_complex, "arccoth", rb_gsl_complex_arccoth, -1);

  rb_define_method(cgsl_complex, "to_s", rb_gsl_complex_to_s, 0);
  rb_define_method(cgsl_complex, "inspect", rb_gsl_complex_inspect, 0);
}
