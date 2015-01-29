/*
  root.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
  (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_array.h"
#include "include/rb_gsl_common.h"
#include "include/rb_gsl_function.h"
#include "include/rb_gsl_root.h"

EXTERN VALUE cgsl_function_fdf;

enum {
  GSL_ROOT_FSOLVER_BISECTION,
  GSL_ROOT_FSOLVER_FALSEPOS,
  GSL_ROOT_FSOLVER_BRENT,
  GSL_ROOT_FDFSOLVER_NEWTON,
  GSL_ROOT_FDFSOLVER_SECANT,
  GSL_ROOT_FDFSOLVER_STEFFENSON,
};

/*****/

static VALUE rb_gsl_fsolver_new(VALUE klass, VALUE t)
{
  gsl_root_fsolver *s = NULL;
  const gsl_root_fsolver_type *T;
  char name[32];
  switch (TYPE(t)) {
  case T_STRING:
    strcpy(name, STR2CSTR(t));
    if (!str_tail_grep(name, "bisection")) {
      T = gsl_root_fsolver_bisection;
    } else if (!str_tail_grep(name, "falsepos")) {
      T = gsl_root_fsolver_falsepos;
    } else if (!str_tail_grep(name, "brent")) {
      T = gsl_root_fsolver_brent;
    } else {
      rb_raise(rb_eTypeError,
         "type must be \"bisection\" or \"falsepos\", or \"brent\".");
    }
    break;
  case T_FIXNUM:
    switch (FIX2INT(t)) {
    case GSL_ROOT_FSOLVER_BISECTION:
      T = gsl_root_fsolver_bisection;
      break;
    case GSL_ROOT_FSOLVER_FALSEPOS:
      T = gsl_root_fsolver_falsepos;
      break;
    case GSL_ROOT_FSOLVER_BRENT:
      T = gsl_root_fsolver_brent;
      break;
    default:
      rb_raise(rb_eTypeError, "type must be BISECTION or FALSEPOS, or BRENT.");
      break;
    }
    break;
  default:
    rb_raise(rb_eTypeError, "wrong argument type %s (String or Fixnum expected)",
       rb_class2name(CLASS_OF(t)));
    break;
  }
  s = gsl_root_fsolver_alloc(T);
  return Data_Wrap_Struct(klass, 0, gsl_root_fsolver_free, s);
}

static VALUE rb_gsl_fsolver_set(VALUE obj, VALUE func, VALUE xl, VALUE xh)
{
  gsl_root_fsolver *s = NULL;
  gsl_function *fff = NULL;
  double xlow, xup;
  Need_Float(xl); Need_Float(xh);
  CHECK_FUNCTION(func);
  Data_Get_Struct(obj, gsl_root_fsolver, s);
  Data_Get_Struct(func, gsl_function, fff);
  xlow = NUM2DBL(xl);
  xup = NUM2DBL(xh);
  gsl_root_fsolver_set(s, fff, xlow, xup);
  return obj;
}

static VALUE rb_gsl_fsolver_iterate(VALUE obj)
{
  gsl_root_fsolver *s = NULL;
  Data_Get_Struct(obj, gsl_root_fsolver, s);
  return INT2FIX(gsl_root_fsolver_iterate(s));
}

static VALUE rb_gsl_fsolver_root(VALUE obj)
{
  gsl_root_fsolver *s = NULL;
  Data_Get_Struct(obj, gsl_root_fsolver, s);
  return rb_float_new(gsl_root_fsolver_root(s));
}

static VALUE rb_gsl_fsolver_x_lower(VALUE obj)
{
  gsl_root_fsolver *s = NULL;
  Data_Get_Struct(obj, gsl_root_fsolver, s);
  return rb_float_new(gsl_root_fsolver_x_lower(s));
}

static VALUE rb_gsl_fsolver_x_upper(VALUE obj)
{
  gsl_root_fsolver *s = NULL;
  Data_Get_Struct(obj, gsl_root_fsolver, s);
  return rb_float_new(gsl_root_fsolver_x_upper(s));
}

static VALUE rb_gsl_fsolver_name(VALUE obj)
{
  gsl_root_fsolver *s = NULL;
  Data_Get_Struct(obj, gsl_root_fsolver, s);
  return rb_str_new2(gsl_root_fsolver_name(s));
}

static VALUE rb_gsl_fsolver_test_interval(VALUE obj, VALUE eabs, VALUE erel)
{
  gsl_root_fsolver *s = NULL;
  Need_Float(eabs); Need_Float(erel);
  Data_Get_Struct(obj, gsl_root_fsolver, s);
  return INT2FIX(gsl_root_test_interval(s->x_lower, s->x_upper,
          NUM2DBL(eabs), NUM2DBL(erel)));
}

static VALUE rb_gsl_root_test_interval(VALUE obj, VALUE xl, VALUE xu, VALUE eabs,
               VALUE erel)
{
  Need_Float(xl); Need_Float(xu);
  Need_Float(eabs); Need_Float(erel);
  return INT2FIX(gsl_root_test_interval(NUM2DBL(xl), NUM2DBL(xu),
          NUM2DBL(eabs), NUM2DBL(erel)));
}

static VALUE rb_gsl_root_test_delta(VALUE obj, VALUE xl, VALUE xu, VALUE eabs,
               VALUE erel)
{
  Need_Float(xl); Need_Float(xu);
  Need_Float(eabs); Need_Float(erel);
  return INT2FIX(gsl_root_test_delta(NUM2DBL(xl), NUM2DBL(xu),
             NUM2DBL(eabs), NUM2DBL(erel)));
}

static VALUE rb_gsl_root_test_residual(VALUE obj, VALUE xl,VALUE eabs)
{
  Need_Float(xl); Need_Float(eabs);
  return INT2FIX(gsl_root_test_residual(NUM2DBL(xl),  NUM2DBL(eabs)));
}

static VALUE rb_gsl_fsolver_solve(int argc, VALUE *argv, VALUE *obj)
{
  gsl_root_fsolver *s = NULL;
  gsl_function *F = NULL;
  double x, xl, xh, epsabs = 0.0, epsrel = 1e-6;
  int status, iter = 0, max_iter = 100;
  switch (argc) {
  case 3:
    Check_Type(argv[2], T_ARRAY);
    epsabs = NUM2DBL(rb_ary_entry(argv[2], 0));
    epsrel = NUM2DBL(rb_ary_entry(argv[2], 1));
    /* no break */
  case 2:
    Check_Type(argv[1], T_ARRAY);
    xl = NUM2DBL(rb_ary_entry(argv[1], 0));
    xh = NUM2DBL(rb_ary_entry(argv[1], 1));
    break;
  default:
    rb_raise(rb_eArgError,
       "Usage: solve(f = Function, range = Array, eps = Array)");
    break;
  }
  CHECK_FUNCTION(argv[0]);
  Data_Get_Struct(argv[0], gsl_function, F);
  Data_Get_Struct(obj, gsl_root_fsolver, s);
  gsl_root_fsolver_set(s, F, xl, xh);
  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    x = gsl_root_fsolver_root (s);
    xl = gsl_root_fsolver_x_lower (s);
    xh = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (xl, xh, epsabs, epsrel);
    if (status == GSL_SUCCESS) break;
  } while (status == GSL_CONTINUE && iter < max_iter);
  return rb_ary_new3(3, rb_float_new(x), INT2FIX(iter), INT2FIX(status));
}

static VALUE rb_gsl_fdfsolver_new(VALUE klass, VALUE t)
{
  gsl_root_fdfsolver *s = NULL;
  const gsl_root_fdfsolver_type *T;
  char name[32];
  switch (TYPE(t)) {
  case T_STRING:
    strcpy(name, STR2CSTR(t));
    if (!str_tail_grep(name, "newton")) {
      T = gsl_root_fdfsolver_newton;
    } else if (!str_tail_grep(name, "secant")) {
      T = gsl_root_fdfsolver_secant;
    } else if (!str_tail_grep(name, "steffenson")) {
      T = gsl_root_fdfsolver_steffenson;
    } else {
      rb_raise(rb_eTypeError, "type must be NEWTON or SECANT, or STEFFENSON.");
    }
    break;
  case T_FIXNUM:
    switch (FIX2INT(t)) {
    case GSL_ROOT_FDFSOLVER_NEWTON:
      T = gsl_root_fdfsolver_newton;
      break;
    case GSL_ROOT_FDFSOLVER_SECANT:
      T = gsl_root_fdfsolver_secant;
      break;
    case GSL_ROOT_FDFSOLVER_STEFFENSON:
      T = gsl_root_fdfsolver_steffenson;
      break;
    default:
      rb_raise(rb_eTypeError, "type must be NEWTON or SECANT, or STEFFENSON.");
      break;
    }
    break;
  default:
    rb_raise(rb_eTypeError, "wrong argument type %s (String or Fixnum expected)",
       rb_class2name(CLASS_OF(t)));
    break;
  }
  s = gsl_root_fdfsolver_alloc(T);
  return Data_Wrap_Struct(klass, 0, gsl_root_fdfsolver_free, s);
}

static VALUE rb_gsl_fdfsolver_set(VALUE obj, VALUE func, VALUE r)
{
  gsl_root_fdfsolver *s = NULL;
  gsl_function_fdf *fff = NULL;
  double root;
  CHECK_FUNCTION_FDF(func);
  Data_Get_Struct(obj, gsl_root_fdfsolver, s);
  Data_Get_Struct(func, gsl_function_fdf, fff);
  root = NUM2DBL(r);
  gsl_root_fdfsolver_set(s, fff, root);
  return obj;
}

static VALUE rb_gsl_fdfsolver_iterate(VALUE obj)
{
  gsl_root_fdfsolver *s = NULL;
  Data_Get_Struct(obj, gsl_root_fdfsolver, s);
  return INT2FIX(gsl_root_fdfsolver_iterate(s));
}

static VALUE rb_gsl_fdfsolver_root(VALUE obj)
{
  gsl_root_fdfsolver *s = NULL;
  Data_Get_Struct(obj, gsl_root_fdfsolver, s);
  return rb_float_new(gsl_root_fdfsolver_root(s));
}

static VALUE rb_gsl_fdfsolver_name(VALUE obj)
{
  gsl_root_fdfsolver *s = NULL;
  Data_Get_Struct(obj, gsl_root_fdfsolver, s);
  return rb_str_new2(gsl_root_fdfsolver_name(s));
}

static VALUE rb_gsl_fdfsolver_solve(int argc, VALUE *argv, VALUE *obj)
{
  gsl_root_fdfsolver *s = NULL;
  double x = 0.0, x0, epsabs = 0.0, epsrel = 1e-6;
  gsl_function_fdf *F = NULL;
  int status, iter = 0, max_iter = 100;
  switch (argc) {
  case 3:
    Check_Type(argv[2], T_ARRAY);
    epsabs = NUM2DBL(rb_ary_entry(argv[2], 0));
    epsrel = NUM2DBL(rb_ary_entry(argv[2], 1));
    /* no break */
  case 2:
    Need_Float(argv[1]);
    x0 = NUM2DBL(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "Usage: solve(f = Function, range = Array, eps = Array)");
    break;
  }
  CHECK_FUNCTION_FDF(argv[0]);
  Data_Get_Struct(argv[0], gsl_function_fdf, F);
  Data_Get_Struct(obj, gsl_root_fdfsolver, s);
  gsl_root_fdfsolver_set(s, F, x0);
  do {
    iter++;
    status = gsl_root_fdfsolver_iterate (s);
    x0 = x;
    x = gsl_root_fdfsolver_root (s);
    status = gsl_root_test_delta(x, x0, epsabs, epsrel);
    if (status == GSL_SUCCESS) break;
  } while (status == GSL_CONTINUE && iter < max_iter);
  return rb_ary_new3(3, rb_float_new(x), INT2FIX(iter), INT2FIX(status));
}

static VALUE rb_gsl_function_rootfinder(int argc, VALUE *argv, VALUE obj)
{
  int status, iter = 0, max_iter = 1000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s = NULL;
  gsl_function *F = NULL;
  double r, a, b, epsabs = 0.0, epsrel = 1e-6;
  Data_Get_Struct(obj, gsl_function, F);
  switch (argc) {
  case 2:
    a = NUM2DBL(argv[0]);
    b = NUM2DBL(argv[1]);
    break;
  case 1:
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      a = NUM2DBL(rb_ary_entry(argv[0], 0));
      b = NUM2DBL(rb_ary_entry(argv[0], 1));
      break;
    default:
      rb_raise(rb_eTypeError, "interval must be given by an array [a, b]");
    }
    break;
  default:
    rb_raise(rb_eArgError, "interval must be given");
    break;
  }
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, F, a, b);
  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    a = gsl_root_fsolver_x_lower (s);
    b = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (a, b, epsabs, epsrel);
    if (status == GSL_SUCCESS) break;
  } while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free(s);
  if (status == GSL_SUCCESS) {
    return rb_ary_new3(3, rb_float_new(r), INT2FIX(iter), INT2FIX(status));
  } else {
    printf("not converged\n");
    return Qfalse;
  }
}

void Init_gsl_root(VALUE module)
{
  VALUE mgsl_root;
  VALUE cgsl_fsolver;
  VALUE cgsl_fdfsolver;

  mgsl_root = rb_define_module_under(module, "Root");

  cgsl_fsolver = rb_define_class_under(mgsl_root, "FSolver", cGSL_Object);
  rb_define_singleton_method(cgsl_fsolver, "alloc", rb_gsl_fsolver_new, 1);

  rb_define_method(cgsl_fsolver, "set", rb_gsl_fsolver_set, 3);
  rb_define_method(cgsl_fsolver, "iterate", rb_gsl_fsolver_iterate, 0);
  rb_define_method(cgsl_fsolver, "root", rb_gsl_fsolver_root, 0);
  rb_define_method(cgsl_fsolver, "name", rb_gsl_fsolver_name, 0);
  rb_define_method(cgsl_fsolver, "x_lower", rb_gsl_fsolver_x_lower, 0);
  rb_define_method(cgsl_fsolver, "x_upper", rb_gsl_fsolver_x_upper, 0);
  rb_define_method(cgsl_fsolver, "test_interval", rb_gsl_fsolver_test_interval, 2);
  rb_define_method(cgsl_fsolver, "solve", rb_gsl_fsolver_solve, -1);

  rb_define_singleton_method(mgsl_root, "test_interval",
          rb_gsl_root_test_interval, 4);
  rb_define_singleton_method(mgsl_root, "test_delta",
          rb_gsl_root_test_delta, 4);
  rb_define_singleton_method(mgsl_root, "test_residual",
          rb_gsl_root_test_residual, 2);

  cgsl_fdfsolver = rb_define_class_under(mgsl_root, "FdfSolver", cGSL_Object);
  rb_define_singleton_method(cgsl_fdfsolver, "alloc", rb_gsl_fdfsolver_new, 1);

  rb_define_method(cgsl_fdfsolver, "set", rb_gsl_fdfsolver_set, 2);
  rb_define_method(cgsl_fdfsolver, "iterate", rb_gsl_fdfsolver_iterate, 0);
  rb_define_method(cgsl_fdfsolver, "root", rb_gsl_fdfsolver_root, 0);
  rb_define_method(cgsl_fdfsolver, "name", rb_gsl_fdfsolver_name, 0);
  rb_define_method(cgsl_fdfsolver, "solve", rb_gsl_fdfsolver_solve, -1);

  rb_define_method(cgsl_function, "fsolve", rb_gsl_function_rootfinder, -1);
  rb_define_alias(cgsl_function, "solve", "fsolve");

  rb_define_const(cgsl_fsolver, "BISECTION", INT2FIX(GSL_ROOT_FSOLVER_BISECTION));
  rb_define_const(cgsl_fsolver, "FALSEPOS", INT2FIX(GSL_ROOT_FSOLVER_FALSEPOS));
  rb_define_const(cgsl_fsolver, "BRENT", INT2FIX(GSL_ROOT_FSOLVER_BRENT));
  rb_define_const(cgsl_fsolver, "Bisection", INT2FIX(GSL_ROOT_FSOLVER_BISECTION));
  rb_define_const(cgsl_fsolver, "Falsepos", INT2FIX(GSL_ROOT_FSOLVER_FALSEPOS));
  rb_define_const(cgsl_fsolver, "Brent", INT2FIX(GSL_ROOT_FSOLVER_BRENT));

  rb_define_const(cgsl_fdfsolver, "NEWTON", INT2FIX(GSL_ROOT_FDFSOLVER_NEWTON));
  rb_define_const(cgsl_fdfsolver, "SECANT", INT2FIX(GSL_ROOT_FDFSOLVER_SECANT));
  rb_define_const(cgsl_fdfsolver, "STEFFENSON", INT2FIX(GSL_ROOT_FDFSOLVER_STEFFENSON));
  rb_define_const(cgsl_fdfsolver, "Newton", INT2FIX(GSL_ROOT_FDFSOLVER_NEWTON));
  rb_define_const(cgsl_fdfsolver, "Secant", INT2FIX(GSL_ROOT_FDFSOLVER_SECANT));
  rb_define_const(cgsl_fdfsolver, "Steffenson", INT2FIX(GSL_ROOT_FDFSOLVER_STEFFENSON));
}
