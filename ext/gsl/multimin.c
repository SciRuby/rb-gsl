/*
  multimin.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl.h"
#include "include/rb_gsl_array.h"
#include "include/rb_gsl_function.h"
#include <gsl/gsl_multimin.h>

#ifndef CHECK_MULTIMIN_FUNCTION
#define CHECK_MULTIMIN_FUNCTION(x) if(CLASS_OF(x)!=cgsl_multimin_function)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (GSL::MultiMin::Function expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

#ifndef CHECK_MULTIMIN_FUNCTION_FDF
#define CHECK_MULTIMIN_FUNCTION_FDF(x) if(CLASS_OF(x)!=cgsl_multimin_function_fdf)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (GSL::MultiMin::Function_fdf expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

VALUE cgsl_multimin_function_fdf;  /* Used also in multimin_fsdf.c */
static VALUE cgsl_multimin_function;

enum {
  GSL_FDFMINIMIZER_CONJUGATE_FR,
  GSL_FDFMINIMIZER_CONJUGATE_PR,
  GSL_FDFMINIMIZER_VECTOR_BFGS,
  GSL_FDFMINIMIZER_STEEPEST_DESCENT,
  GSL_FMINIMIZER_NMSIMPLEX,
#ifdef GSL_1_9_LATER
  GSL_FDFMINIMIZER_VECTOR_BFGS2,
#endif
#ifdef GSL_1_13_LATER
  GSL_FMINIMIZER_NMSIMPLEX2RAND,
#endif
};

static const gsl_multimin_fdfminimizer_type* get_fdfminimizer_type(VALUE t);
#ifdef GSL_1_3_LATER
static const gsl_multimin_fminimizer_type* get_fminimizer_type(VALUE t);
#endif
static void define_const();

static void gsl_multimin_function_free(gsl_multimin_function *f);
static double rb_gsl_multimin_function_f(const gsl_vector *x, void *p);
static void set_function(int i, VALUE *argv, gsl_multimin_function *F);

static void gsl_multimin_function_fdf_free(gsl_multimin_function_fdf *f);

double rb_gsl_multimin_function_fdf_f(const gsl_vector *x, void *p); 
void rb_gsl_multimin_function_fdf_df(const gsl_vector *x, void *p, 
					   gsl_vector *g);
void rb_gsl_multimin_function_fdf_fdf(const gsl_vector *x, void *p, 
					    double *f, gsl_vector *g);
static void set_function_fdf(int i, VALUE *argv, gsl_multimin_function_fdf *F);

/*** multimin_funcion ***/
static void gsl_multimin_function_mark(gsl_multimin_function *F);
static void gsl_multimin_function_fdf_mark(gsl_multimin_function_fdf *F);

static void gsl_multimin_function_mark(gsl_multimin_function *F)
{
  rb_gc_mark((VALUE) F->params);
}

static void gsl_multimin_function_fdf_mark(gsl_multimin_function_fdf *F)
{
  rb_gc_mark((VALUE) F->params);
}

static VALUE rb_gsl_multimin_function_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_multimin_function *F = NULL;
  VALUE ary;
  size_t i;
  F = ALLOC(gsl_multimin_function);
  F->f = &rb_gsl_multimin_function_f;
  ary = rb_ary_new2(2);
  /*  (VALUE) F->params = ary;*/
  F->params = (void *) ary;
  if (rb_block_given_p()) rb_ary_store(ary, 0, rb_block_proc());
  else rb_ary_store(ary, 0, Qnil);
  rb_ary_store(ary, 1, Qnil);
  switch (argc) {
  case 0:
    break;
  case 1:
    set_function(0, argv, F);
    break;
  case 2:
  case 3:
    for (i = 0; (int) i < argc; i++) set_function(i, argv, F);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
  }
  return Data_Wrap_Struct(klass, gsl_multimin_function_mark, gsl_multimin_function_free, F);
}

static void gsl_multimin_function_free(gsl_multimin_function *f)
{
  free((gsl_multimin_function *) f);
}

static VALUE rb_gsl_multimin_function_n(VALUE obj)
{
  gsl_multimin_function *F = NULL;
  Data_Get_Struct(obj, gsl_multimin_function, F);
  return INT2FIX(F->n);
}

static double rb_gsl_multimin_function_f(const gsl_vector *x, void *p)
{
  VALUE vx, vp, proc, result;
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  proc = rb_ary_entry((VALUE) p, 0);
  vp = rb_ary_entry((VALUE) p, 1);
  if (NIL_P(vp)) result = rb_funcall(proc, RBGSL_ID_call, 1, vx);
  else result = rb_funcall(proc, RBGSL_ID_call, 2, vx, vp);
  return NUM2DBL(result);
}

static VALUE rb_gsl_multimin_function_eval(VALUE obj, VALUE vx)
{
  gsl_multimin_function *F = NULL;
  VALUE vp, proc, ary, result;
  Data_Get_Struct(obj, gsl_multimin_function, F);
  ary = (VALUE) F->params;
  proc = rb_ary_entry(ary, 0);
  vp = rb_ary_entry(ary, 1);
  if (NIL_P(vp)) result = rb_funcall(proc, RBGSL_ID_call, 1, vx);
  else result = rb_funcall(proc, RBGSL_ID_call, 2, vx, vp);
  return result;
}

static void set_function(int i, VALUE *argv, gsl_multimin_function *F)
{
  VALUE ary;
  ary = (VALUE) F->params;
  if (TYPE(argv[i]) == T_FIXNUM) F->n = FIX2INT(argv[i]);
  else if (rb_obj_is_kind_of(argv[i], rb_cProc)) 
    rb_ary_store(ary, 0, argv[i]);
  else if (TYPE(argv[i]) == T_ARRAY || rb_obj_is_kind_of(argv[i], cgsl_vector)
		|| TYPE(argv[i]) == T_FIXNUM || TYPE(argv[i]) == T_FLOAT) {
    rb_ary_store(ary, 1, argv[i]);
  } else {
    rb_raise(rb_eTypeError, "wrong type of argument %d (Fixnum or Proc)", i);
  }
}

static VALUE rb_gsl_multimin_function_set_f(int argc, VALUE *argv, VALUE obj)
{
  gsl_multimin_function *F = NULL;
  VALUE ary;
  size_t i;
  Data_Get_Struct(obj, gsl_multimin_function, F);
  ary = (VALUE) F->params;
  if (rb_block_given_p()) rb_ary_store(ary, 0, rb_block_proc());
  switch (argc) {
  case 1:
    set_function(0, argv, F);
    break;
  case 2:
  case 3:
    for (i = 0; (int) i < argc; i++) set_function(i, argv, F);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
  }
  return obj;
}

static VALUE rb_gsl_multimin_function_set_params(int argc, VALUE *argv, VALUE obj)
{
  gsl_multimin_function *F = NULL;
  VALUE ary, ary2;
  size_t i;
  if (argc == 0) return obj;
  Data_Get_Struct(obj, gsl_multimin_function, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  if (argc == 1) rb_ary_store(ary, 1, argv[0]);
  else {
    ary2 = rb_ary_new2(argc);
    for (i = 0; (int) i < argc; i++) rb_ary_store(ary2, i, argv[i]);
    rb_ary_store(ary, 1, ary2);
  }
  return obj;
}

static VALUE rb_gsl_multimin_function_params(VALUE obj)
{
  gsl_multimin_function *F = NULL;
  Data_Get_Struct(obj, gsl_multimin_function, F);
  return rb_ary_entry((VALUE) F->params, 1);
}

/*** multimin_function_fdf ***/
static void set_function_fdf(int argc, VALUE *argv, gsl_multimin_function_fdf *F);
static VALUE rb_gsl_multimin_function_fdf_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_multimin_function_fdf *F = NULL;
  VALUE ary;
  F = ALLOC(gsl_multimin_function_fdf);
  F->f = &rb_gsl_multimin_function_fdf_f;
  F->df = &rb_gsl_multimin_function_fdf_df;
  F->fdf = &rb_gsl_multimin_function_fdf_fdf;
  ary = rb_ary_new2(4);
  /*  (VALUE) F->params = ary;*/
  F->params = (void *) ary;
  rb_ary_store(ary, 2, Qnil);
  rb_ary_store(ary, 3, Qnil);
  set_function_fdf(argc, argv, F);
  return Data_Wrap_Struct(klass, gsl_multimin_function_fdf_mark, gsl_multimin_function_fdf_free, F);
}

static void gsl_multimin_function_fdf_free(gsl_multimin_function_fdf *f)
{
  free((gsl_multimin_function_fdf *) f);
}

static VALUE rb_gsl_multimin_function_fdf_n(VALUE obj)
{
  gsl_multimin_function_fdf *F = NULL;
  Data_Get_Struct(obj, gsl_multimin_function_fdf, F);
  return INT2FIX(F->n);
}

static void set_function_fdf(int argc, VALUE *argv, gsl_multimin_function_fdf *F)
{
  VALUE ary;
  if (F->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  switch (argc) {
  case 1:
    if (TYPE(argv[0]) != T_FIXNUM) rb_raise(rb_eTypeError, "Fixnum expected");
    F->n = FIX2INT(argv[0]);
    break;
  case 2:
    rb_ary_store(ary, 0, argv[0]);
    rb_ary_store(ary, 1, argv[1]);
    rb_ary_store(ary, 2, Qnil);
    break;
  case 3:
    rb_ary_store(ary, 0, argv[0]);
    rb_ary_store(ary, 1, argv[1]);
    if (TYPE(argv[2]) == T_FIXNUM) {
      F->n = FIX2INT(argv[2]);
      rb_ary_store(ary, 2, Qnil);
    } else {
			rb_ary_store(ary, 2, argv[2]);
		}
    break;
  case 4:
  case 5:
    if (TYPE(argv[0]) == T_FIXNUM) {
      F->n = FIX2INT(argv[0]);
      rb_ary_store(ary, 0, argv[1]);
      rb_ary_store(ary, 1, argv[2]);
      rb_ary_store(ary, 2, argv[3]);
    } else {
      rb_ary_store(ary, 0, argv[0]);
      rb_ary_store(ary, 1, argv[1]);
      rb_ary_store(ary, 2, argv[2]);
      F->n = FIX2INT(argv[3]);
    }
    if (argc == 5) rb_ary_store(ary, 3, argv[4]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (1, 3, or 4)");
  }
}

static VALUE rb_gsl_multimin_function_fdf_set_procs(int argc, VALUE *argv, VALUE obj)
{
  VALUE ary;
  gsl_multimin_function_fdf *F = NULL;
  Data_Get_Struct(obj, gsl_multimin_function_fdf, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  switch (argc) {
  case 2:
    rb_ary_store(ary, 0, argv[0]);
    rb_ary_store(ary, 1, argv[1]);
    rb_ary_store(ary, 2, Qnil);
    break;
  case 3:
    rb_ary_store(ary, 0, argv[0]);
    rb_ary_store(ary, 1, argv[1]);
    if (TYPE(argv[2]) == T_FIXNUM) {
      F->n = FIX2INT(argv[2]);
      rb_ary_store(ary, 2, Qnil);
    } else {
			rb_ary_store(ary, 2, argv[2]);
		}
    break;
  case 4:
  case 5:
    if (TYPE(argv[0]) == T_FIXNUM) {
      F->n = FIX2INT(argv[0]);
      rb_ary_store(ary, 0, argv[1]);
      rb_ary_store(ary, 1, argv[2]);
      rb_ary_store(ary, 2, argv[3]);
    } else {
      rb_ary_store(ary, 0, argv[0]);
      rb_ary_store(ary, 1, argv[1]);
      rb_ary_store(ary, 2, argv[2]);
      F->n = FIX2INT(argv[3]);
    }
    if (argc == 5) rb_ary_store(ary, 3, argv[4]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (2, 3, or 4)");
  }
  return obj;
}

static VALUE rb_gsl_multimin_function_fdf_set_params(int argc, VALUE *argv, VALUE obj)
{
  gsl_multimin_function_fdf *F = NULL;
  VALUE ary, ary2;
  size_t i;
  if (argc == 0) return obj;
  Data_Get_Struct(obj, gsl_multimin_function_fdf, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  if (argc == 1) {
		rb_ary_store(ary, 3, argv[0]);
  } else {
    ary2 = rb_ary_new2(argc);
    for (i = 0; (int) i < argc; i++) rb_ary_store(ary2, i, argv[i]);
    rb_ary_store(ary, 3, ary2);
  }
  return obj;
}

static VALUE rb_gsl_multimin_function_fdf_set(int argc, VALUE *argv, VALUE obj)
{
  gsl_multimin_function_fdf *F = NULL;
  Data_Get_Struct(obj, gsl_multimin_function_fdf, F);
  set_function_fdf(argc, argv, F);
  return obj;
}

double rb_gsl_multimin_function_fdf_f(const gsl_vector *x, void *p)
{
  VALUE vx, proc, vp, result, ary;
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  ary = (VALUE) p;
  proc = rb_ary_entry(ary, 0);
  vp = rb_ary_entry(ary, RARRAY_LEN(ary)-1);
  if (NIL_P(vp)) result = rb_funcall(proc, RBGSL_ID_call, 1, vx);
  else result = rb_funcall(proc, RBGSL_ID_call, 2, vx, vp);
  return NUM2DBL(result);
}

void rb_gsl_multimin_function_fdf_df(const gsl_vector *x, void *p, 
				     gsl_vector *g)
{
  VALUE vx, vg, proc, vp, ary;
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  vg = Data_Wrap_Struct(cgsl_vector, 0, NULL, g);
  ary = (VALUE) p;
  proc = rb_ary_entry(ary, 1);
  vp = rb_ary_entry(ary, RARRAY_LEN(ary)-1);
  if (NIL_P(vp)) {
    rb_funcall(proc, RBGSL_ID_call, 2, vx, vg);
  } else {
    rb_funcall(proc, RBGSL_ID_call, 3, vx, vp, vg);
  }
}

void rb_gsl_multimin_function_fdf_fdf(const gsl_vector *x, void *p, 
				      double *f, gsl_vector *g)
{
  VALUE vx, vg, proc_f, proc_df, vp, ary, result;
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  vg = Data_Wrap_Struct(cgsl_vector, 0, NULL, g);
  ary = (VALUE) p;
  proc_f = rb_ary_entry(ary, 0);
  proc_df = rb_ary_entry(ary, 1);
  vp = rb_ary_entry(ary, RARRAY_LEN(ary)-1);
  if (NIL_P(vp)) {
    result = rb_funcall(proc_f, RBGSL_ID_call, 1, vx);
    rb_funcall(proc_df, RBGSL_ID_call, 2, vx, vg);
  } else {
    result = rb_funcall(proc_f, RBGSL_ID_call, 2, vx, vp);
    rb_funcall(proc_df, RBGSL_ID_call, 3, vx, vp, vg);
  }
  *f = NUM2DBL(result);
}

static VALUE rb_gsl_multimin_function_fdf_params(VALUE obj)
{
  gsl_multimin_function_fdf *F = NULL;
  Data_Get_Struct(obj, gsl_multimin_function_fdf, F);
  return rb_ary_entry((VALUE) F->params, 3);
}

/****************/

static void define_const(VALUE klass1, VALUE klass2)
{
  rb_define_const(klass1, 
		  "CONJUGATE_FR", INT2FIX(GSL_FDFMINIMIZER_CONJUGATE_FR));
  rb_define_const(klass1, 
		  "CONJUGATE_PR", INT2FIX(GSL_FDFMINIMIZER_CONJUGATE_PR));
  rb_define_const(klass1, 
		  "VECTOR_BFGS", INT2FIX(GSL_FDFMINIMIZER_VECTOR_BFGS));  
  rb_define_const(klass1, 
		  "STEEPEST_DESCENT", INT2FIX(GSL_FDFMINIMIZER_STEEPEST_DESCENT));
#ifdef GSL_1_3_LATER
  rb_define_const(klass2, 
		  "NMSIMPLEX", INT2FIX(GSL_FMINIMIZER_NMSIMPLEX));
#endif
#ifdef GSL_1_9_LATER
  rb_define_const(klass1, 
		  "VECTOR_BFGS2", INT2FIX(GSL_FDFMINIMIZER_VECTOR_BFGS2));
#endif
#ifdef GSL_1_13_LATER
  rb_define_const(klass2, 
		  "NMSIMPLEX2RAND", INT2FIX(GSL_FMINIMIZER_NMSIMPLEX2RAND));
#endif
}

static const gsl_multimin_fdfminimizer_type* get_fdfminimizer_type(VALUE t)
{
  char name[64];
  switch (TYPE(t)) {
  case T_STRING:
    strcpy(name, STR2CSTR(t));
    if (str_tail_grep(name, "conjugate_fr") == 0) 
      return gsl_multimin_fdfminimizer_conjugate_fr;
    else if (str_tail_grep(name, "conjugate_pr") == 0) 
      return gsl_multimin_fdfminimizer_conjugate_pr;
    else if (str_tail_grep(name, "vector_bfgs") == 0) 
      return gsl_multimin_fdfminimizer_vector_bfgs;
    else if (str_tail_grep(name, "steepest_descent") == 0) 
      return gsl_multimin_fdfminimizer_steepest_descent;
#ifdef GSL_1_9_LATER
    else if (str_tail_grep(name, "vector_bfgs2") == 0) 
      return gsl_multimin_fdfminimizer_vector_bfgs2;
#endif
    else
      rb_raise(rb_eTypeError, "%s: unknown minimizer type", name);
    break;
  case T_FIXNUM:
    switch (FIX2INT(t)) {
    case GSL_FDFMINIMIZER_CONJUGATE_FR:
      return gsl_multimin_fdfminimizer_conjugate_fr; break;
    case GSL_FDFMINIMIZER_CONJUGATE_PR:
      return gsl_multimin_fdfminimizer_conjugate_pr; break;
    case GSL_FDFMINIMIZER_VECTOR_BFGS:
      return gsl_multimin_fdfminimizer_vector_bfgs; break;
    case GSL_FDFMINIMIZER_STEEPEST_DESCENT:
      return gsl_multimin_fdfminimizer_steepest_descent; break;
#ifdef GSL_1_9_LATER
    case GSL_FDFMINIMIZER_VECTOR_BFGS2:
      return gsl_multimin_fdfminimizer_vector_bfgs2; break;
#endif      
    default:
      rb_raise(rb_eTypeError, "%d: unknown type", FIX2INT(t));
      break;
    }
    break;
  default:
    rb_raise(rb_eTypeError, "type is given by a String or a Fixnum");
    break;
  }
}

static VALUE rb_gsl_fdfminimizer_new(VALUE klass, VALUE t, VALUE n)
{
  gsl_multimin_fdfminimizer *gmf = NULL;
  const gsl_multimin_fdfminimizer_type *T;
  T = get_fdfminimizer_type(t);
  gmf = gsl_multimin_fdfminimizer_alloc(T, FIX2INT(n));
  return Data_Wrap_Struct(klass, 0, gsl_multimin_fdfminimizer_free, gmf);
}

static VALUE rb_gsl_fdfminimizer_set(VALUE obj, VALUE ff, VALUE xx, VALUE ss, 
				     VALUE tt)
{
  gsl_multimin_fdfminimizer *gmf = NULL;
  gsl_multimin_function_fdf *F = NULL;
  gsl_vector *x;
  double stepsize, tol;
  int status;
  CHECK_MULTIMIN_FUNCTION_FDF(ff);
  Need_Float(ss); Need_Float(tt);
  Data_Get_Struct(obj, gsl_multimin_fdfminimizer, gmf);
  Data_Get_Struct(ff, gsl_multimin_function_fdf, F);
  Data_Get_Vector(xx, x);
  stepsize = NUM2DBL(ss);
  tol = NUM2DBL(tt);
  status = gsl_multimin_fdfminimizer_set(gmf, F, x, stepsize, tol);
  return INT2FIX(status);
}

static VALUE rb_gsl_fdfminimizer_name(VALUE obj)
{
  gsl_multimin_fdfminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fdfminimizer, gmf);
  return rb_str_new2(gsl_multimin_fdfminimizer_name(gmf));
}

static VALUE rb_gsl_fdfminimizer_iterate(VALUE obj)
{
  gsl_multimin_fdfminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fdfminimizer, gmf);
  return INT2FIX(gsl_multimin_fdfminimizer_iterate(gmf));
}

static VALUE rb_gsl_fdfminimizer_x(VALUE obj)
{
  gsl_multimin_fdfminimizer *gmf = NULL;
  gsl_vector *x = NULL;
  Data_Get_Struct(obj, gsl_multimin_fdfminimizer, gmf);
  x = gsl_multimin_fdfminimizer_x(gmf);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, x);
}

static VALUE rb_gsl_fdfminimizer_gradient(VALUE obj)
{
  gsl_multimin_fdfminimizer *gmf = NULL;
  gsl_vector *gradient = NULL;
  Data_Get_Struct(obj, gsl_multimin_fdfminimizer, gmf);
  gradient = gsl_multimin_fdfminimizer_gradient(gmf);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, gradient);
}

static VALUE rb_gsl_fdfminimizer_minimum(VALUE obj)
{
  gsl_multimin_fdfminimizer *gmf = NULL;
  double min;
  Data_Get_Struct(obj, gsl_multimin_fdfminimizer, gmf);
  min = gsl_multimin_fdfminimizer_minimum(gmf);
  return rb_float_new(min);
}

static VALUE rb_gsl_fdfminimizer_f(VALUE obj)
{
  gsl_multimin_fdfminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fdfminimizer, gmf);
  return rb_float_new(gmf->f);
}

static VALUE rb_gsl_fdfminimizer_restart(VALUE obj)
{
  gsl_multimin_fdfminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fdfminimizer, gmf);
  return INT2FIX(gsl_multimin_fdfminimizer_restart(gmf));
}

static VALUE rb_gsl_fdfminimizer_test_gradient(VALUE obj, VALUE ea)
{
  gsl_multimin_fdfminimizer *gmf = NULL;
  gsl_vector *g = NULL;
  Need_Float(ea);
  Data_Get_Struct(obj, gsl_multimin_fdfminimizer, gmf);
  g = gsl_multimin_fdfminimizer_gradient(gmf);
  return INT2FIX(gsl_multimin_test_gradient(g, NUM2DBL(ea)));
}

static VALUE rb_gsl_multimin_test_gradient(VALUE obj, VALUE gg, VALUE ea)
{
  gsl_vector *g = NULL;
  Need_Float(ea);
  Data_Get_Vector(gg, g);
  return INT2FIX(gsl_multimin_test_gradient(g, NUM2DBL(ea)));
}

/*****/
#ifdef GSL_1_3_LATER
static const gsl_multimin_fminimizer_type* get_fminimizer_type(VALUE t)
{
  char name[64];

  switch (TYPE(t)) {
  case T_STRING:
    strcpy(name, STR2CSTR(t));
    if (str_tail_grep(name, "nmsimplex") == 0) 
      return gsl_multimin_fminimizer_nmsimplex;
#ifdef GSL_1_13_LATER
    if (str_tail_grep(name, "nmsimplex2rand") == 0) 
      return gsl_multimin_fminimizer_nmsimplex2rand;
#endif
    else
      rb_raise(rb_eTypeError, "unknown type %s (nmsimplex and nmsimplex2rand supported)", name);
    break;
  case T_FIXNUM:
    switch (FIX2INT(t)) {
    case GSL_FMINIMIZER_NMSIMPLEX:
      return gsl_multimin_fminimizer_nmsimplex; break;
#ifdef GSL_1_13_LATER
    case GSL_FMINIMIZER_NMSIMPLEX2RAND:
      return gsl_multimin_fminimizer_nmsimplex2rand; break;
#endif
    default:
      rb_raise(rb_eTypeError, "%d: unknown type (not supported)", FIX2INT(t));
      break;
    }
    break;
  default:
    rb_raise(rb_eTypeError, "wrong argument type %s (String or Fixnum expected)",
	     rb_class2name(CLASS_OF(t)));
    break;
  }
}

static VALUE rb_gsl_fminimizer_new(VALUE klass, VALUE t, VALUE n)
{
  gsl_multimin_fminimizer *gmf = NULL;
  const gsl_multimin_fminimizer_type *T;
  CHECK_FIXNUM(n);
  T = get_fminimizer_type(t);
  gmf = gsl_multimin_fminimizer_alloc(T, FIX2INT(n));
  return Data_Wrap_Struct(klass, 0, gsl_multimin_fminimizer_free, gmf);
}

static VALUE rb_gsl_fminimizer_set(VALUE obj, VALUE ff, VALUE xx, VALUE ss)
{
  gsl_multimin_fminimizer *gmf = NULL;
  gsl_multimin_function *F = NULL;
  gsl_vector *x = NULL, *s = NULL;
  CHECK_MULTIMIN_FUNCTION(ff);
  Data_Get_Struct(obj, gsl_multimin_fminimizer, gmf);
  Data_Get_Struct(ff, gsl_multimin_function, F);
  Data_Get_Vector(xx, x);
  Data_Get_Vector(ss, s);
  return INT2FIX(gsl_multimin_fminimizer_set(gmf, F, x, s));
}

static VALUE rb_gsl_fminimizer_name(VALUE obj)
{
  gsl_multimin_fminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fminimizer, gmf);
  return rb_str_new2(gsl_multimin_fminimizer_name(gmf));
}

static VALUE rb_gsl_fminimizer_iterate(VALUE obj)
{
  gsl_multimin_fminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fminimizer, gmf);
  return INT2FIX(gsl_multimin_fminimizer_iterate(gmf));
}

static VALUE rb_gsl_fminimizer_x(VALUE obj)
{
  gsl_multimin_fminimizer *gmf = NULL;
  gsl_vector *x = NULL;
  Data_Get_Struct(obj, gsl_multimin_fminimizer, gmf);
  x = gsl_multimin_fminimizer_x(gmf);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, x);
}

static VALUE rb_gsl_fminimizer_minimum(VALUE obj)
{
  gsl_multimin_fminimizer *gmf = NULL;
  double min;
  Data_Get_Struct(obj, gsl_multimin_fminimizer, gmf);
  min = gsl_multimin_fminimizer_minimum(gmf);
  return rb_float_new(min);
}

static VALUE rb_gsl_fminimizer_size(VALUE obj)
{
  gsl_multimin_fminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fminimizer, gmf);
  return rb_float_new(gsl_multimin_fminimizer_size(gmf));
}

static VALUE rb_gsl_fminimizer_fval(VALUE obj)
{
  gsl_multimin_fminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fminimizer, gmf);
  return rb_float_new(gmf->fval);
}

static VALUE rb_gsl_fminimizer_test_size(VALUE obj, VALUE ea)
{
  gsl_multimin_fminimizer *gmf = NULL;
  Need_Float(ea);
  Data_Get_Struct(obj, gsl_multimin_fminimizer, gmf);
  return INT2FIX(gsl_multimin_test_size(gmf->size, NUM2DBL(ea)));
}
static VALUE rb_gsl_multimin_test_size(VALUE obj, VALUE ss, VALUE ea)
{
  Need_Float(ss); Need_Float(ea);
  return INT2FIX(gsl_multimin_test_size(NUM2DBL(ss), NUM2DBL(ea)));
}

#endif

#ifdef HAVE_GSL_GSL_MULTIMIN_FSDF_H
void Init_multimin_fsdf(VALUE module);
#endif

/*****/
void Init_gsl_multimin(VALUE module)
{
  VALUE mgsl_multimin;
  VALUE cgsl_multimin_fdfminimizer;
  VALUE cgsl_multimin_fminimizer;

  mgsl_multimin = rb_define_module_under(module, "MultiMin");
  rb_define_singleton_method(mgsl_multimin, "test_gradient", rb_gsl_multimin_test_gradient, 2);
#ifdef GSL_1_3_LATER
  rb_define_singleton_method(mgsl_multimin, "test_size", rb_gsl_multimin_test_size, 2);
#endif

  cgsl_multimin_fdfminimizer = rb_define_class_under(mgsl_multimin, "FdfMinimizer", cGSL_Object);
  cgsl_multimin_fminimizer = rb_define_class_under(mgsl_multimin, "FMinimizer", cGSL_Object);
  define_const(cgsl_multimin_fdfminimizer, cgsl_multimin_fminimizer);

  cgsl_multimin_function = rb_define_class_under(mgsl_multimin, "Function",
						  cgsl_function);
  rb_define_singleton_method(cgsl_multimin_function, "alloc",
			     rb_gsl_multimin_function_new, -1);
  rb_define_method(cgsl_multimin_function, "eval", rb_gsl_multimin_function_eval, 1);
  rb_define_alias(cgsl_multimin_function, "call", "eval");
  rb_define_method(cgsl_multimin_function, "set_proc", rb_gsl_multimin_function_set_f, -1);
  rb_define_alias(cgsl_multimin_function, "set_f", "set_proc");
  rb_define_method(cgsl_multimin_function, "set_params", rb_gsl_multimin_function_set_params, -1);
  rb_define_method(cgsl_multimin_function, "params", rb_gsl_multimin_function_params, 0);
  rb_define_method(cgsl_multimin_function, "n", rb_gsl_multimin_function_n, 0);

  cgsl_multimin_function_fdf = rb_define_class_under(mgsl_multimin, "Function_fdf",
						     cGSL_Object);
  rb_define_singleton_method(cgsl_multimin_function_fdf, "alloc",
			     rb_gsl_multimin_function_fdf_new, -1);

  rb_define_method(cgsl_multimin_function_fdf, "set", rb_gsl_multimin_function_fdf_set, -1);
  rb_define_method(cgsl_multimin_function_fdf, "set_params", rb_gsl_multimin_function_fdf_set_params, -1);
  rb_define_method(cgsl_multimin_function_fdf, "set_procs", rb_gsl_multimin_function_fdf_set_procs, -1);
  rb_define_method(cgsl_multimin_function_fdf, "params", rb_gsl_multimin_function_fdf_params, 0);
  rb_define_method(cgsl_multimin_function_fdf, "n", rb_gsl_multimin_function_fdf_n, 0);

  rb_define_singleton_method(cgsl_multimin_fdfminimizer, "alloc", rb_gsl_fdfminimizer_new, 2);

  rb_define_method(cgsl_multimin_fdfminimizer, "set", rb_gsl_fdfminimizer_set, 4);
  rb_define_method(cgsl_multimin_fdfminimizer, "name", rb_gsl_fdfminimizer_name, 0);
  rb_define_method(cgsl_multimin_fdfminimizer, "iterate", rb_gsl_fdfminimizer_iterate, 0);
  rb_define_method(cgsl_multimin_fdfminimizer, "x", rb_gsl_fdfminimizer_x, 0);
  rb_define_method(cgsl_multimin_fdfminimizer, "f", rb_gsl_fdfminimizer_f, 0);
  rb_define_method(cgsl_multimin_fdfminimizer, "gradient", rb_gsl_fdfminimizer_gradient, 0);
  rb_define_method(cgsl_multimin_fdfminimizer, "minimum", rb_gsl_fdfminimizer_minimum, 0);
  rb_define_method(cgsl_multimin_fdfminimizer, "restart", rb_gsl_fdfminimizer_restart, 0);
  rb_define_method(cgsl_multimin_fdfminimizer, "test_gradient", rb_gsl_fdfminimizer_test_gradient, 1);

  /*****/
#ifdef GSL_1_3_LATER
  rb_define_singleton_method(cgsl_multimin_fminimizer, "alloc", rb_gsl_fminimizer_new, 2);

  rb_define_method(cgsl_multimin_fminimizer, "set", rb_gsl_fminimizer_set, 3);
  rb_define_method(cgsl_multimin_fminimizer, "name", rb_gsl_fminimizer_name, 0);
  rb_define_method(cgsl_multimin_fminimizer, "iterate", rb_gsl_fminimizer_iterate, 0);
  rb_define_method(cgsl_multimin_fminimizer, "x", rb_gsl_fminimizer_x, 0);
  rb_define_method(cgsl_multimin_fminimizer, "fval", rb_gsl_fminimizer_fval, 0);
  rb_define_method(cgsl_multimin_fminimizer, "minimum", rb_gsl_fminimizer_minimum, 0);  
    rb_define_method(cgsl_multimin_fminimizer, "size", rb_gsl_fminimizer_size, 0);
  rb_define_method(cgsl_multimin_fminimizer, "test_size", rb_gsl_fminimizer_test_size, 1);
#endif


#ifdef HAVE_GSL_GSL_MULTIMIN_FSDF_H
	Init_multimin_fsdf(mgsl_multimin);
#endif
}
#ifdef CHECK_MULTIMIN_FUNCTION
#undef CHECK_MULTIMIN_FUNCTION
#endif
#ifdef CHECK_MULTIMIN_FUNCTION_FDF
#undef CHECK_MULTIMIN_FUNCTION_FDF
#endif

