/*
  function.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_function.h"
#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif

VALUE cgsl_function;
VALUE cgsl_function_fdf;

void gsl_function_free(gsl_function *f);
double rb_gsl_function_f(double x, void *p); 
ID RBGSL_ID_call, RBGSL_ID_arity;

static VALUE rb_gsl_function_set_f(int argc, VALUE *argv, VALUE obj)
{
  gsl_function *F = NULL;
  VALUE ary, ary2;
  size_t i;
  Data_Get_Struct(obj, gsl_function, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(2);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  rb_ary_store(ary, 1, Qnil);

  switch (argc) {
  case 0:
    break;
  case 1:
    CHECK_PROC(argv[0]);
    rb_ary_store(ary, 0, argv[0]);
    break;
  case 2:
    CHECK_PROC(argv[0]);
    rb_ary_store(ary, 0, argv[0]);
    rb_ary_store(ary, 1, argv[1]);
    break;
  default:
    CHECK_PROC(argv[0]);
    rb_ary_store(ary, 0, argv[0]);
    ary2 = rb_ary_new2(argc-1);
    for (i = 1; (int) i < argc; i++) rb_ary_store(ary2, i-1, argv[i]);
    rb_ary_store(ary, 1, ary2);
    break;
  }
  if (rb_block_given_p()) rb_ary_store(ary, 0, rb_block_proc());
  return obj;
}

void gsl_function_free(gsl_function *f)
{
  if (f) free((gsl_function *) f);
}

void gsl_function_mark(gsl_function *f)
{
  rb_gc_mark((VALUE) f->params);
}

/*
 * Create a Function object
 */
static VALUE rb_gsl_function_alloc(int argc, VALUE *argv, VALUE klass)
{
  gsl_function *f = NULL;
  VALUE obj;
  f = ALLOC(gsl_function);
  f->function = &rb_gsl_function_f;
  /*  (VALUE) f->params = rb_ary_new2(2);*/
  f->params = (void *) rb_ary_new2(2);
  rb_ary_store((VALUE) f->params, 1, Qnil);
  obj = Data_Wrap_Struct(klass, gsl_function_mark, gsl_function_free, f);
  rb_gsl_function_set_f(argc, argv, obj);
  return obj;
}			    

double rb_gsl_function_f(double x, void *p)
{
  VALUE result, ary, proc, params;
  ary = (VALUE) p;
  proc = rb_ary_entry(ary, 0);
  params = rb_ary_entry(ary, 1);
  if (NIL_P(params)) result = rb_funcall(proc, RBGSL_ID_call, 1, rb_float_new(x));
  else result = rb_funcall(proc, RBGSL_ID_call, 2, rb_float_new(x), params);
  return NUM2DBL(result);
}

/*
 * Calculates a function at x, and returns the rusult.
 */
static VALUE rb_gsl_function_eval(VALUE obj, VALUE x)
{
  gsl_function *F = NULL;
  VALUE ary, proc, params, result, arynew, x2;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  size_t i, j, n;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  Data_Get_Struct(obj, gsl_function, F);
  ary = (VALUE) F->params;
  proc = rb_ary_entry(ary, 0);
  params = rb_ary_entry(ary, 1);
  if (CLASS_OF(x) == rb_cRange) x = rb_gsl_range2ary(x);
  switch (TYPE(x)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    if (NIL_P(params)) result = rb_funcall(proc, RBGSL_ID_call, 1, x);
    else result = rb_funcall(proc, RBGSL_ID_call, 2, x, params);
    return result;
    break;
  case T_ARRAY:
    //    n = RARRAY(x)->len;
    n = RARRAY_LEN(x);
    arynew = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      x2 = rb_ary_entry(x, i);
      Need_Float(x2);
      if (NIL_P(params)) result = rb_funcall(proc, RBGSL_ID_call, 1, x2);
      else result = rb_funcall(proc, RBGSL_ID_call, 2, x2, params);
      rb_ary_store(arynew, i, result);
    }
    return arynew;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(x)) {
      GetNArray(x, na);
      ptr1 = (double *) na->ptr;
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(x));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) {
	x2 = rb_float_new(ptr1[i]);
	if (NIL_P(params)) result = rb_funcall(proc, RBGSL_ID_call, 1, x2);
	else result = rb_funcall(proc, RBGSL_ID_call, 2, x2, params);
	ptr2[i] = NUM2DBL(result);
      }
      return ary;
    }
#endif
    if (VECTOR_P(x)) {
      Data_Get_Struct(x, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	x2 = rb_float_new(gsl_vector_get(v, i));
	if (NIL_P(params)) result = rb_funcall(proc, RBGSL_ID_call, 1, x2);
	else result = rb_funcall(proc, RBGSL_ID_call, 2, x2, params);
	gsl_vector_set(vnew, i, NUM2DBL(result));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MATRIX_P(x)) {
      Data_Get_Struct(x, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  x2 = rb_float_new(gsl_matrix_get(m, i, j));
	  if (NIL_P(params)) result = rb_funcall(proc, RBGSL_ID_call, 1, x2);
	  else result = rb_funcall(proc, RBGSL_ID_call, 2, x2, params);
	  gsl_matrix_set(mnew, i, j, NUM2DBL(result));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  }
  /* never reach here */
  return Qnil;
}

static VALUE rb_gsl_function_arity(VALUE obj)
{
  gsl_function *F = NULL;
  VALUE proc;
  Data_Get_Struct(obj, gsl_function, F);
  proc = rb_ary_entry((VALUE) F->params, 0);
  return INT2FIX(rb_funcall(proc, RBGSL_ID_arity, 0));
}

static VALUE rb_gsl_function_proc(VALUE obj)
{
  gsl_function *F = NULL;
  Data_Get_Struct(obj, gsl_function, F);
  return rb_ary_entry((VALUE) F->params, 0);
}

static VALUE rb_gsl_function_params(VALUE obj)
{
  gsl_function *F = NULL;
  Data_Get_Struct(obj, gsl_function, F);
  return rb_ary_entry((VALUE) F->params, 1);
}

static VALUE rb_gsl_function_set_params(int argc, VALUE *argv, VALUE obj)
{
  gsl_function *F = NULL;
  VALUE ary, ary2;
  size_t i;
  if (argc == 0) return obj;
  Data_Get_Struct(obj, gsl_function, F);
  ary = (VALUE) F->params;
  if (argc == 1) {
    rb_ary_store(ary, 1, argv[0]);
  } else {
    ary2 = rb_ary_new2(argc);
    for (i = 0; (int) i < argc; i++) rb_ary_store(ary2, i, argv[i]);
    rb_ary_store(ary, 1, ary2);
  }
  return obj;
}

static VALUE rb_gsl_function_graph(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_GNU_GRAPH
  gsl_function *F = NULL;
  gsl_vector *v = NULL;
  double x, y;
  char opt[256] = "", command[1024];
  size_t i, n;
  int flag = 0;
  FILE *fp = NULL;
  VALUE ary, params, proc;
  switch (argc) {
  case 2:
    Check_Type(argv[1], T_STRING);
    strcpy(opt, STR2CSTR(argv[1]));
    /* no break, do next */
  case 1:
    if (CLASS_OF(argv[0]) == rb_cRange) argv[0] = rb_gsl_range2ary(argv[0]);
    if (TYPE(argv[0]) == T_ARRAY) {
      //      n = RARRAY(argv[0])->len;
      n = RARRAY_LEN(argv[0]);
      v = gsl_vector_alloc(n);
      flag = 1;
      for (i = 0; i < n; i++) 
	gsl_vector_set(v, i, NUM2DBL(rb_ary_entry(argv[0], i)));
    } else if (rb_obj_is_kind_of(argv[0], cgsl_vector)) {
      Data_Get_Struct(argv[0], gsl_vector, v);
      n = v->size;
      flag = 0;
    } else {
      rb_raise(rb_eTypeError, 
	       "wrong argument type %s (Array or GSL::Vector expected)", 
	       rb_class2name(CLASS_OF(argv[0])));
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
    break;
  }
  Data_Get_Struct(obj, gsl_function, F);
  ary = (VALUE) F->params;
  proc = rb_ary_entry(ary, 0);
  params = rb_ary_entry(ary, 1);
  sprintf(command, "graph -T X -g 3 %s", opt);
  fp = popen(command, "w");
  if (fp == NULL)
    rb_raise(rb_eIOError, "GNU graph not found.");
  for (i = 0; i < n; i++) {
    x = gsl_vector_get(v, i);
    if (NIL_P(params)) y = NUM2DBL(rb_funcall(proc, RBGSL_ID_call, 1, rb_float_new(x)));
    else y = NUM2DBL(rb_funcall(proc, RBGSL_ID_call, 2, rb_float_new(x), params));
    fprintf(fp, "%e %e\n", x, y);
  }
  fflush(fp);
 pclose(fp);
  fp = NULL;
  if (flag == 1) gsl_vector_free(v);
  return Qtrue;
#else
  rb_raise(rb_eNoMethodError, "not implemented");
  return Qfalse;
#endif
}


static double rb_gsl_function_fdf_f(double x, void *p);
static void gsl_function_fdf_free(gsl_function_fdf *f);

static double rb_gsl_function_fdf_f(double x, void *p);
static double rb_gsl_function_fdf_df(double x, void *p);
static void rb_gsl_function_fdf_fdf(double x, void *p, double *f, double *df);

static void setfunc(int i, VALUE *argv, gsl_function_fdf *F);
static void setfunc(int i, VALUE *argv, gsl_function_fdf *F)
{
  VALUE ary;
  if (F->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }

  if (rb_obj_is_kind_of(argv[i], rb_cProc)) {
    rb_ary_store(ary, i, argv[i]);
  } else if (TYPE(argv[i]) == T_ARRAY || rb_obj_is_kind_of(argv[i], cgsl_vector) 
	     || TYPE(argv[i]) == T_FIXNUM || TYPE(argv[i]) == T_FLOAT) {
    rb_ary_store(ary, 3, argv[i]);
  } else {
    rb_raise(rb_eArgError, 
	     "wrong type argument (Proc, Array, GSL::Vector or a number)");
  }
}

static void gsl_function_fdf_mark(gsl_function_fdf *f);
static VALUE rb_gsl_function_fdf_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_function_fdf *F = NULL;
  VALUE ary;
  size_t i;
  F = ALLOC(gsl_function_fdf);
  F->f = &rb_gsl_function_fdf_f;
  F->df = &rb_gsl_function_fdf_df;
  F->fdf = &rb_gsl_function_fdf_fdf;
  ary = rb_ary_new2(4);
  /*  (VALUE) F->params = ary;*/
  F->params = (void *) ary;
  rb_ary_store(ary, 2, Qnil);
  rb_ary_store(ary, 3, Qnil);
  for (i = 0; (int) i < argc; i++) setfunc(i, argv, F);
  return Data_Wrap_Struct(klass, gsl_function_fdf_mark, gsl_function_fdf_free, F);
}

static void gsl_function_fdf_free(gsl_function_fdf *f)
{
  free((gsl_function_fdf *) f);
}

static void gsl_function_fdf_mark(gsl_function_fdf *f)
{
  rb_gc_mark((VALUE) f->params);
}

static VALUE rb_gsl_function_fdf_set(int argc, VALUE *argv, VALUE obj)
{
  gsl_function_fdf *F = NULL;
  VALUE ary;
  size_t i;
  Data_Get_Struct(obj, gsl_function_fdf, F);
  ary = (VALUE) F->params;
  rb_ary_store(ary, 2, Qnil);
  rb_ary_store(ary, 3, Qnil);
  for (i = 0; (int) i < argc; i++) setfunc(i, argv, F);
  return obj;
}

static VALUE rb_gsl_function_fdf_set_f(VALUE obj, VALUE procf)
{
  gsl_function_fdf *F = NULL;
  VALUE ary;
  CHECK_PROC(procf);
  Data_Get_Struct(obj, gsl_function_fdf, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  rb_ary_store(ary, 0, procf);
  return obj;
}

static VALUE rb_gsl_function_fdf_set_df(VALUE obj, VALUE procdf)
{
  gsl_function_fdf *F = NULL;
  VALUE ary;
  CHECK_PROC(procdf);
  Data_Get_Struct(obj, gsl_function_fdf, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  rb_ary_store(ary, 1, procdf);
  return obj;
}

static VALUE rb_gsl_function_fdf_set_fdf(VALUE obj, VALUE procfdf)
{
  gsl_function_fdf *F = NULL;
  VALUE ary;
  CHECK_PROC(procfdf);
  Data_Get_Struct(obj, gsl_function_fdf, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  rb_ary_store(ary, 2, procfdf);
  return obj;
}

static VALUE rb_gsl_function_fdf_set_params(int argc, VALUE *argv, VALUE obj)
{
  gsl_function_fdf *F = NULL;
  VALUE ary, ary2;
  size_t i;
  Data_Get_Struct(obj, gsl_function_fdf, F);
  ary = (VALUE) F->params;
  if (argc == 0) return obj;
  if (argc == 1) {
    rb_ary_store(ary, 3, argv[0]);
  } else {
    ary2 = rb_ary_new2(argc);
    for (i = 0; (int) i < argc; i++) rb_ary_store(ary2, i, argv[i]);
    rb_ary_store(ary, 3, ary2);
  }
  return obj;
}

static double rb_gsl_function_fdf_f(double x, void *p)
{
  VALUE result, params, proc, ary;
  ary = (VALUE) p;
  proc = rb_ary_entry(ary, 0);
  params = rb_ary_entry(ary, 3);
  if (NIL_P(params)) result = rb_funcall(proc, RBGSL_ID_call, 1, rb_float_new(x));
  else result = rb_funcall(proc, RBGSL_ID_call, 2, rb_float_new(x), params);
  return NUM2DBL(result);
}

static double rb_gsl_function_fdf_df(double x, void *p)
{
  VALUE result, params, proc, ary;
  ary = (VALUE) p;
  proc = rb_ary_entry(ary, 1);
  params = rb_ary_entry(ary, 3);
  if (NIL_P(params)) result = rb_funcall(proc, RBGSL_ID_call, 1, rb_float_new(x));
  else result = rb_funcall(proc, RBGSL_ID_call, 2, rb_float_new(x), params);
  return NUM2DBL(result);
}

static void rb_gsl_function_fdf_fdf(double x, void *p, double *f, double *df)
{
  VALUE result, params, proc_f, proc_df, proc_fdf, ary;
  ary = (VALUE) p;
  proc_f = rb_ary_entry(ary, 0);
  proc_df = rb_ary_entry(ary, 1);
  proc_fdf = rb_ary_entry(ary, 2);
  params = rb_ary_entry(ary, 3);
  if (NIL_P(proc_fdf)) {
    if (NIL_P(params)) {
      result = rb_funcall(proc_f, RBGSL_ID_call, 1, rb_float_new(x));
      *f = NUM2DBL(result);
      result = rb_funcall(proc_df, RBGSL_ID_call, 1, rb_float_new(x));
      *df = NUM2DBL(result);
    } else {
      result = rb_funcall(proc_f, RBGSL_ID_call, 2, rb_float_new(x), params);
      *f = NUM2DBL(result);
      result = rb_funcall(proc_df, RBGSL_ID_call, 2, rb_float_new(x), params);
      *df = NUM2DBL(result);
    }
  } else {
    if (NIL_P(params)) result = rb_funcall(proc_fdf, RBGSL_ID_call, 1, rb_float_new(x));
    else result = rb_funcall(proc_fdf, RBGSL_ID_call, 2, rb_float_new(x), params);
    *f = NUM2DBL(rb_ary_entry(result, 0));
    *df = NUM2DBL(rb_ary_entry(result, 1));
  }
}

void Init_gsl_function(VALUE module)
{
  RBGSL_ID_call = rb_intern("call");
  RBGSL_ID_arity = rb_intern("arity");

  cgsl_function = rb_define_class_under(module, "Function", cGSL_Object);
  cgsl_function_fdf = rb_define_class_under(module, "Function_fdf", cGSL_Object);
  // This Fdf class seems superfluous.  Should probably be deleted?
  rb_define_class_under(cgsl_function_fdf, "Fdf", cgsl_function_fdf);

  /*  rb_define_singleton_method(cgsl_function, "new", rb_gsl_function_new, -1);*/
  rb_define_singleton_method(cgsl_function, "alloc", rb_gsl_function_alloc, -1);

  rb_define_method(cgsl_function, "eval", rb_gsl_function_eval, 1);
  rb_define_alias(cgsl_function, "call", "eval");
  rb_define_alias(cgsl_function, "[]", "eval");
  rb_define_alias(cgsl_function, "at", "eval");
  rb_define_method(cgsl_function, "arity", rb_gsl_function_arity, 0);
  rb_define_method(cgsl_function, "proc", rb_gsl_function_proc, 0);
  rb_define_alias(cgsl_function, "f", "proc");
  rb_define_method(cgsl_function, "params", rb_gsl_function_params, 0);
  rb_define_alias(cgsl_function, "param", "params");
  rb_define_method(cgsl_function, "set", rb_gsl_function_set_f, -1);
  rb_define_method(cgsl_function, "set_params", rb_gsl_function_set_params, -1);
  rb_define_alias(cgsl_function, "set_param", "set_params");
  rb_define_alias(cgsl_function, "params=", "set_params");
  rb_define_alias(cgsl_function, "param=", "set_params");

  rb_define_method(cgsl_function, "graph", rb_gsl_function_graph, -1);
  /*****/
  rb_define_singleton_method(cgsl_function_fdf, "new", rb_gsl_function_fdf_new, -1);
  rb_define_singleton_method(cgsl_function_fdf, "alloc", rb_gsl_function_fdf_new, -1);
  rb_define_method(cgsl_function_fdf, "set", rb_gsl_function_fdf_set, -1);
  rb_define_method(cgsl_function_fdf, "set_f", rb_gsl_function_fdf_set_f, 1);
  rb_define_method(cgsl_function_fdf, "set_df", rb_gsl_function_fdf_set_df, 1);
  rb_define_method(cgsl_function_fdf, "set_fdf", rb_gsl_function_fdf_set_fdf, 1);
  rb_define_method(cgsl_function_fdf, "set_params", rb_gsl_function_fdf_set_params, -1);

}
