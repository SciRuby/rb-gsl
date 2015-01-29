/*
  integration.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_array.h"
#include "include/rb_gsl_function.h"
#include "include/rb_gsl_integration.h"
#include "include/rb_gsl_common.h"

#ifndef CHECK_WORKSPACE
#define CHECK_WORKSPACE(x) if(CLASS_OF(x)!=cgsl_integration_workspace)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (Integration::Workspace expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

#define EPSABS_DEFAULT 0.0
#define EPSREL_DEFAULT 1e-10
#define LIMIT_DEFAULT 1000
#define KEY_DEFAULT GSL_INTEG_GAUSS61

static VALUE cgsl_integration_qaws_table, cgsl_integration_qawo_table;

static VALUE cgsl_integration_workspace;

#ifdef GSL_1_14_LATER
static VALUE cgsl_integration_glfixed_table;
#endif

static int get_a_b(int argc, VALUE *argv, int argstart, double *a, double *b);
static int get_epsabs_epsrel(int argc, VALUE *argv, int argstart, 
           double *epsabs, double *epsrel);
static int get_a_b_epsabs_epsrel(int argc, VALUE *argv, int argstart,
          double *a, double *b, double *epsabs,
          double *epsrel);
static int get_limit_key_workspace(int argc, VALUE *argv, int argstart,
           size_t *limit, int *key, 
           gsl_integration_workspace **w);
static int get_limit_workspace(int argc, VALUE *argv, int argstart,
           size_t *limit, 
             gsl_integration_workspace **w);
static int get_epsabs_epsrel_limit_workspace(int argc, VALUE *argv, int argstart,
               double *epsabs, double *epsrel,
               size_t *limit,
               gsl_integration_workspace **w);

static int get_a_b(int argc, VALUE *argv, int argstart, double *a, double *b)
{
  int itmp;
  VALUE aa, bb;
  if (argstart >= argc) return argstart;
  if (TYPE(argv[argstart]) == T_ARRAY) {
    aa = rb_ary_entry(argv[argstart], 0);
    bb = rb_ary_entry(argv[argstart], 1);
    Need_Float(aa); Need_Float(bb);
    //    *a = RFLOAT(aa)->value;
    //    *b = RFLOAT(bb)->value;
    *a = NUM2DBL(aa);
    *b = NUM2DBL(bb);
    itmp = argstart + 1;
  } else {
    Need_Float(argv[argstart]); Need_Float(argv[argstart+1]);
    *a = NUM2DBL(argv[argstart]);
    *b = NUM2DBL(argv[argstart+1]);
    itmp = argstart + 2;
  }
  return itmp;
}

static int get_epsabs_epsrel(int argc, VALUE *argv, int argstart, 
           double *epsabs, double *epsrel)
{
  int itmp;
  VALUE aa, bb;
  *epsabs = EPSABS_DEFAULT;
  *epsrel = EPSREL_DEFAULT;
  if (argstart >= argc) return argstart;
  if (TYPE(argv[argstart]) == T_ARRAY) {
    aa = rb_ary_entry(argv[argstart], 0);
    bb = rb_ary_entry(argv[argstart], 1);
    Need_Float(aa); Need_Float(bb);
    *epsabs = NUM2DBL(aa);
    *epsrel = NUM2DBL(bb);
    itmp = 1;
  } else {
    Need_Float(argv[argstart]); Need_Float(argv[argstart+1]);
    *epsabs = NUM2DBL(argv[argstart]);
    *epsrel = NUM2DBL(argv[argstart+1]);
    itmp = 2;
  } 
  return  argstart + itmp;
}

static int get_a_b_epsabs_epsrel(int argc, VALUE *argv, int argstart,
          double *a, double *b, double *epsabs,
          double *epsrel)
{
  int itmp;
  *epsabs = EPSABS_DEFAULT;
  *epsrel = EPSREL_DEFAULT;
  itmp = get_a_b(argc, argv, argstart, a, b);
  itmp = get_epsabs_epsrel(argc, argv, itmp, epsabs, epsrel);
  return itmp;
}

static int get_limit_key_workspace(int argc, VALUE *argv, int argstart,
           size_t *limit, int *key, 
           gsl_integration_workspace **w)
{
  int flag = 0;
  switch (argc-argstart) {
  case 3:
    CHECK_FIXNUM(argv[argstart]);
    CHECK_FIXNUM(argv[argstart+1]);
    CHECK_WORKSPACE(argv[argstart+2]);
    *limit = FIX2INT(argv[argstart]);
    *key = FIX2INT(argv[argstart+1]);
    Data_Get_Struct(argv[argstart+2], gsl_integration_workspace, *w);
    flag = 0;
    break;
  case 1:
    CHECK_FIXNUM(argv[argstart]);
    *key = FIX2INT(argv[argstart]);
    *limit = LIMIT_DEFAULT;
    *w = gsl_integration_workspace_alloc(*limit);
    flag = 1;
    break;
  case 2:
    if (TYPE(argv[argc-1]) == T_FIXNUM) {
      CHECK_FIXNUM(argv[argc-2]);
      *limit = FIX2INT(argv[argc-2]);
      *key = FIX2INT(argv[argc-1]);
      *w = gsl_integration_workspace_alloc(*limit);
      flag = 1;
    } else {
      CHECK_FIXNUM(argv[argc-2]);
      CHECK_WORKSPACE(argv[argc-1]);
      *key = FIX2INT(argv[argc-2]);
      Data_Get_Struct(argv[argc-1], gsl_integration_workspace, *w);
      *limit = (*w)->limit;
      flag = 0;
    }
    break;
  case 0:
    *key = KEY_DEFAULT;
    *limit = LIMIT_DEFAULT;
    *w = gsl_integration_workspace_alloc(*limit);
    flag = 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  if (*w == NULL) rb_raise(rb_eRuntimeError, "something wrong with workspace");
  return flag;
}

static int get_limit_workspace(int argc, VALUE *argv, int argstart,
           size_t *limit, 
           gsl_integration_workspace **w)
{
  int flag = 0;

  switch (argc-argstart) {
  case 2:
    CHECK_FIXNUM(argv[argstart]);
    *limit = FIX2INT(argv[argstart]);
    CHECK_WORKSPACE(argv[argstart+1]);
    Data_Get_Struct(argv[argstart+1], gsl_integration_workspace, *w);
    flag = 0;
    break;
  case 0:
    *limit = LIMIT_DEFAULT;
    *w = gsl_integration_workspace_alloc(*limit);
    flag = 1;
    break;
  case 1:
    switch (TYPE(argv[argstart])) {
    case T_FIXNUM:
    case T_BIGNUM:
      CHECK_FIXNUM(argv[argstart]);
      *limit = FIX2INT(argv[argstart]);
      *w = gsl_integration_workspace_alloc(*limit);
      flag = 1;
      break;
    default:
      CHECK_WORKSPACE(argv[argc-1]);
      Data_Get_Struct(argv[argc-1], gsl_integration_workspace, *w);
      *limit = (*w)->limit;
      flag = 0;
      break;
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  if (*w == NULL) rb_raise(rb_eRuntimeError, "something wrong with workspace");
  return flag;
}

static int get_epsabs_epsrel_limit_workspace(int argc, VALUE *argv, int argstart,
               double *epsabs, double *epsrel,
               size_t *limit,
               gsl_integration_workspace **w)
{
  int flag = 0, itmp;
  itmp = argstart;
  *epsabs = EPSABS_DEFAULT;
  *epsrel = EPSREL_DEFAULT;
  *limit = LIMIT_DEFAULT;
  switch (argc-itmp) {
  case 0:
    *w = gsl_integration_workspace_alloc(*limit);
    flag = 1;
    break;
  case 1:
    if (TYPE(argv[itmp]) == T_ARRAY) {
      get_epsabs_epsrel(argc, argv, itmp, epsabs, epsrel);
      *w = gsl_integration_workspace_alloc(*limit);
      flag = 1;
    } else {
      flag = get_limit_workspace(argc, argv, itmp, limit, w);
    }
    break;
  case 2:
  case 3:
    switch (TYPE(argv[itmp])) {
    case T_ARRAY:
      itmp = get_epsabs_epsrel(argc, argv, itmp, epsabs, epsrel);
      flag = get_limit_workspace(argc, argv, itmp, limit, w);
      break;
    case T_FLOAT:
      get_epsabs_epsrel(argc, argv, itmp, epsabs, epsrel);
      *w = gsl_integration_workspace_alloc(*limit);
      flag = 1;
      break;
    default:
      flag = get_limit_workspace(argc, argv, itmp, limit, w);
      break;
    }
    break;
  case 4:
    itmp = get_epsabs_epsrel(argc, argv, itmp, epsabs, epsrel);
    flag = get_limit_workspace(argc, argv, itmp, limit, w);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  if (*w == NULL) rb_raise(rb_eRuntimeError, "something wrong with workspace");
  return flag;
}

static VALUE rb_gsl_integration_qng(int argc, VALUE *argv, VALUE obj)
{
  double a, b, epsabs = EPSABS_DEFAULT, epsrel = EPSREL_DEFAULT;
  double result, abserr;
  size_t neval;
  gsl_function *F = NULL;
  int status;
  // local variable 'itmp' declared and set, but never used
  //int itmp;

  if (argc < 1) rb_raise(rb_eArgError, 
       "wrong number of arguments (%d for >= 1)", argc);

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    CHECK_FUNCTION(argv[0]);
    Data_Get_Struct(argv[0], gsl_function, F);
    /*itmp =*/ get_a_b_epsabs_epsrel(argc, argv, 1, &a, &b, &epsabs, &epsrel);
    break;
  default:
    /*itmp =*/ get_a_b_epsabs_epsrel(argc, argv, 0, &a, &b, &epsabs, &epsrel);
    Data_Get_Struct(obj, gsl_function, F);
    break;
  }
  status = gsl_integration_qng(F, a, b, epsabs, epsrel, 
             &result, &abserr, &neval);
  
  return rb_ary_new3(4, rb_float_new(result), rb_float_new(abserr),
         INT2FIX(neval), INT2FIX(status));
}          

static VALUE rb_gsl_integration_qag(int argc, VALUE *argv, VALUE obj)
{
  double a, b, epsabs = EPSABS_DEFAULT, epsrel = EPSREL_DEFAULT;
  double result, abserr;
  size_t limit = LIMIT_DEFAULT;
  gsl_function *F = NULL;
  gsl_integration_workspace *w = NULL;
  int key = KEY_DEFAULT, status, intervals, itmp, flag = 0;
  if (argc < 1) rb_raise(rb_eArgError, 
       "wrong number of arguments (%d for >= 1)", argc);

  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    CHECK_FUNCTION(argv[0]);
    Data_Get_Struct(argv[0], gsl_function, F);
    if (argc == 3) {
      CHECK_FIXNUM(argv[2]);
      get_a_b(argc, argv, 1, &a, &b);
      key = FIX2INT(argv[2]);
      w = gsl_integration_workspace_alloc(limit);
      flag = 1;
    } else if (argc == 4) {
      CHECK_FIXNUM(argv[3]);
      get_a_b(argc, argv, 1, &a, &b);
      key = FIX2INT(argv[3]);
      w = gsl_integration_workspace_alloc(limit);
      flag = 1;
    } else {
      itmp = get_a_b_epsabs_epsrel(argc, argv, 1, &a, &b, &epsabs, &epsrel);
      flag = get_limit_key_workspace(argc, argv, itmp, &limit, &key, &w);
    }
    break;
  default:
    if (argc == 2) {
      if (FIXNUM_P(argv[1])) {
  key = FIX2INT(argv[1]);
  w = gsl_integration_workspace_alloc(limit);
  flag = 1;
      } else if (rb_obj_is_kind_of(argv[1], cgsl_integration_workspace)) {
  Data_Get_Struct(argv[1], gsl_integration_workspace, w);
  flag = 0;
      } else {
  rb_raise(rb_eTypeError, "Key or Workspace expected");
      }
      itmp = get_a_b(argc, argv, 0, &a, &b);
    } else if (argc == 3) {
      if (FIXNUM_P(argv[2])) {
  key = FIX2INT(argv[2]);
  w = gsl_integration_workspace_alloc(limit);
  flag = 1;
      } else if (rb_obj_is_kind_of(argv[2], cgsl_integration_workspace)) {
  Data_Get_Struct(argv[2], gsl_integration_workspace, w);
  flag = 0;
      } else {
  rb_raise(rb_eTypeError, "Key or Workspace expected");
      }
      itmp = get_a_b(argc, argv, 0, &a, &b);
    } else {
      itmp = get_a_b_epsabs_epsrel(argc, argv, 0, &a, &b, &epsabs, &epsrel);
      flag = get_limit_key_workspace(argc, argv, itmp, &limit, &key, &w);
    }
    Data_Get_Struct(obj, gsl_function, F);
    break;
  }
  status = gsl_integration_qag(F, a, b, epsabs, epsrel, limit, key, w, 
             &result, &abserr);
  intervals = w->size;
  if (flag == 1) gsl_integration_workspace_free(w);

  return rb_ary_new3(4, rb_float_new(result), rb_float_new(abserr), 
         INT2FIX(intervals), INT2FIX(status));
}          

static VALUE rb_gsl_integration_qags(int argc, VALUE *argv, VALUE obj)
{
  double a, b, epsabs = EPSABS_DEFAULT, epsrel = EPSREL_DEFAULT;
  double result, abserr;
  size_t limit = LIMIT_DEFAULT;
  gsl_function *F = NULL;
  gsl_integration_workspace *w = NULL;
  int status, intervals, flag = 0, itmp;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    CHECK_FUNCTION(argv[0]);
    Data_Get_Struct(argv[0], gsl_function, F);
    itmp = get_a_b(argc, argv, 1, &a, &b);
    break;
  default:
    Data_Get_Struct(obj, gsl_function, F);
    itmp = get_a_b(argc, argv, 0, &a, &b);
    break;
  }
  flag = get_epsabs_epsrel_limit_workspace(argc, argv, itmp, &epsabs, &epsrel,
             &limit, &w);
 
  status = gsl_integration_qags(F, a, b, epsabs, epsrel, limit, w, 
        &result, &abserr);
  intervals = w->size;
  if (flag == 1) gsl_integration_workspace_free(w);

  return rb_ary_new3(4, rb_float_new(result), rb_float_new(abserr), 
      INT2FIX(intervals), INT2FIX(status));
}          

static VALUE rb_gsl_integration_qagp(int argc, VALUE *argv, VALUE obj)
{
  double epsabs, epsrel;
  double result, abserr;
  size_t limit;
  gsl_function *F = NULL;
  gsl_vector *v = NULL;
  gsl_integration_workspace *w = NULL;
  int status, intervals, flag = 0, flag2 = 0, itmp;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    CHECK_FUNCTION(argv[0]);
    Data_Get_Struct(argv[0], gsl_function, F);
    itmp = 1;
    break;
  default:
    Data_Get_Struct(obj, gsl_function, F);
    itmp = 0;
    break;
  }
  if (TYPE(argv[itmp]) == T_ARRAY) {
    v = make_cvector_from_rarray(argv[itmp]);
    flag2 = 1;
  } else {
    Data_Get_Vector(argv[itmp], v);
    flag2 = 0;
  }
  itmp += 1;
  flag = get_epsabs_epsrel_limit_workspace(argc, argv, itmp, &epsabs, &epsrel,
             &limit, &w);

  status = gsl_integration_qagp(F, v->data, v->size, epsabs, epsrel, limit, w, 
        &result, &abserr);
  intervals = w->size;
  if (flag == 1) gsl_integration_workspace_free(w);
  if (flag2 == 1) gsl_vector_free(v);

  return rb_ary_new3(4, rb_float_new(result), rb_float_new(abserr), 
         INT2FIX(intervals), INT2FIX(status));
}          

/* (-infty --- +infty) */
static VALUE rb_gsl_integration_qagi(int argc, VALUE *argv, VALUE obj)
{
  double epsabs, epsrel;
  double result, abserr;
  size_t limit;
  gsl_function *F = NULL;
  gsl_integration_workspace *w = NULL;
  int status, intervals, flag = 0, itmp;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    CHECK_FUNCTION(argv[0]);
    Data_Get_Struct(argv[0], gsl_function, F);
    itmp = 1;
    break;
  default:
    Data_Get_Struct(obj, gsl_function, F);
    itmp = 0;
    break;
  }
  flag = get_epsabs_epsrel_limit_workspace(argc, argv, itmp, &epsabs, &epsrel,
             &limit, &w);
  status = gsl_integration_qagi(F, epsabs, epsrel, limit, w, 
        &result, &abserr);
  intervals = w->size;
  if (flag == 1) gsl_integration_workspace_free(w);

  return rb_ary_new3(4, rb_float_new(result), rb_float_new(abserr), 
        INT2FIX(intervals), INT2FIX(status));
}          

/* (a --- +infty) */
static VALUE rb_gsl_integration_qagiu(int argc, VALUE *argv, VALUE obj)
{
  double a, epsabs, epsrel;
  double result, abserr;
  size_t limit;
  gsl_function *F = NULL;
  gsl_integration_workspace *w = NULL;
  int status, intervals, flag = 0, itmp;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    CHECK_FUNCTION(argv[0]);
    Data_Get_Struct(argv[0], gsl_function, F);
    itmp = 1;
    break;
  default:
    Data_Get_Struct(obj, gsl_function, F);
    itmp = 0;
    break;
  }
  Need_Float(argv[itmp]);
  a = NUM2DBL(argv[itmp]);
  itmp += 1;
  flag = get_epsabs_epsrel_limit_workspace(argc, argv, itmp, &epsabs, &epsrel,
             &limit, &w);
  status = gsl_integration_qagiu(F, a, epsabs, epsrel, limit, w, 
        &result, &abserr);
  intervals = w->size;
  if (flag == 1) gsl_integration_workspace_free(w);

  return rb_ary_new3(4, rb_float_new(result), rb_float_new(abserr), 
         INT2FIX(intervals), INT2FIX(status));
}          

/* (-infty --- b) */
static VALUE rb_gsl_integration_qagil(int argc, VALUE *argv, VALUE obj)
{
  double b, epsabs, epsrel;
  double result, abserr;
  size_t limit;
  gsl_function *F = NULL;
  gsl_integration_workspace *w = NULL;
  int status, intervals, flag = 0, itmp;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    CHECK_FUNCTION(argv[0]);
    Data_Get_Struct(argv[0], gsl_function, F);
    itmp = 1;
    break;
  default:
    Data_Get_Struct(obj, gsl_function, F);
    itmp = 0;
    break;
  }
  Need_Float(argv[itmp]);
  b = NUM2DBL(argv[itmp]);
  flag = get_epsabs_epsrel_limit_workspace(argc, argv, itmp+1, &epsabs, &epsrel,
             &limit, &w);
  Data_Get_Struct(obj, gsl_function, F);

  status = gsl_integration_qagil(F, b, epsabs, epsrel, limit, w, 
        &result, &abserr);
  intervals = w->size;
  if (flag == 1) gsl_integration_workspace_free(w);

  return rb_ary_new3(4, rb_float_new(result), rb_float_new(abserr), 
         INT2FIX(intervals), INT2FIX(status));
}          

static VALUE rb_gsl_integration_qawc(int argc, VALUE *argv, VALUE obj)
{
  double a, b, c, epsabs, epsrel;
  double result, abserr;
  size_t limit;
  gsl_function *F = NULL;
  gsl_integration_workspace *w = NULL;
  int status, intervals, itmp, flag = 0;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    CHECK_FUNCTION(argv[0]);
    Data_Get_Struct(argv[0], gsl_function, F);
    itmp = 1;
    break;
  default:
    Data_Get_Struct(obj, gsl_function, F);
    itmp = 0;
    break;
  }
  itmp = get_a_b(argc, argv, itmp, &a, &b);
  if (argc-itmp <= 0) rb_raise(rb_eArgError, "The pole is not given");
  Need_Float(argv[itmp]);
  c = NUM2DBL(argv[itmp]);
  flag = get_epsabs_epsrel_limit_workspace(argc, argv, itmp+1, &epsabs, &epsrel,
             &limit, &w);
  status = gsl_integration_qawc(F, a, b, c, epsabs, epsrel, limit, w, &result, &abserr);
  intervals = w->size;
  if (flag == 1) gsl_integration_workspace_free(w);

  return rb_ary_new3(4, rb_float_new(result), rb_float_new(abserr), INT2FIX(intervals),
         INT2FIX(status));
}          

VALUE rb_gsl_integration_qaws_table_alloc(int argc, VALUE *argv, VALUE klass)
{
  gsl_integration_qaws_table *t = NULL;
  VALUE alpha, beta, mu, nu;

  if (TYPE(argv[0]) == T_ARRAY) {
    alpha = rb_ary_entry(argv[0], 0);
    beta = rb_ary_entry(argv[0], 1);
    mu = rb_ary_entry(argv[0], 2);
    nu = rb_ary_entry(argv[0], 3);
  } else {
    Need_Float(argv[0]); Need_Float(argv[1]);
    CHECK_FIXNUM(argv[2]); CHECK_FIXNUM(argv[3]);
    alpha = argv[0];
    beta = argv[1];
    mu = argv[2];
    nu = argv[3];
  }
  t = gsl_integration_qaws_table_alloc(NUM2DBL(alpha), NUM2DBL(beta),
               FIX2INT(mu), FIX2INT(nu));
  return Data_Wrap_Struct(klass, 0, gsl_integration_qaws_table_free, t);
}

static VALUE rb_gsl_integration_qaws_table_set(int argc, VALUE *argv, VALUE obj)
{
  gsl_integration_qaws_table *t = NULL;
  double alpha, beta;
  int mu, nu, type;
  if (argc != 1 && argc != 4)
    rb_raise(rb_eArgError, "wrong number of argument (%d for 1 or 3)", argc);
  type = TYPE(argv[0]);
  Data_Get_Struct(obj, gsl_integration_qaws_table, t);

  if (type == T_FIXNUM || type == T_BIGNUM || type == T_FLOAT) {
    alpha = NUM2DBL(argv[0]);
    beta  = NUM2DBL(argv[1]);
    mu    = FIX2INT(argv[2]);
    nu    = FIX2INT(argv[3]);
  } else if (type == T_ARRAY) {
    alpha = NUM2DBL(rb_ary_entry(argv[0], 0));
    beta  = NUM2DBL(rb_ary_entry(argv[0], 1));
    mu    = FIX2INT(rb_ary_entry(argv[0], 2));
    nu    = FIX2INT(rb_ary_entry(argv[0], 3));
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s", rb_class2name(CLASS_OF(argv[0])));
  }

  gsl_integration_qaws_table_set(t, alpha, beta, mu, nu);
  return obj;
}

static VALUE rb_gsl_integration_qaws_table_to_a(VALUE obj)
{
  gsl_integration_qaws_table *t = NULL;
  VALUE ary;
  Data_Get_Struct(obj, gsl_integration_qaws_table, t);
  ary = rb_ary_new2(4);
  rb_ary_store(ary, 0, rb_float_new(t->alpha));
  rb_ary_store(ary, 1, rb_float_new(t->beta));
  rb_ary_store(ary, 2, INT2FIX(t->mu));
  rb_ary_store(ary, 3, INT2FIX(t->nu));
  return ary;
}

static gsl_integration_qaws_table* make_qaws_table(VALUE ary);
static VALUE rb_gsl_ary_to_integration_qaws_table(VALUE ary)
{
  gsl_integration_qaws_table *t = NULL;
  t = make_qaws_table(ary);
  return Data_Wrap_Struct(cgsl_integration_qaws_table, 
        0, gsl_integration_qaws_table_free, t);
}

static gsl_integration_qaws_table* make_qaws_table(VALUE ary)
{
  double alpha, beta;
  int mu, nu;
  alpha = NUM2DBL(rb_ary_entry(ary, 0));
  beta  = NUM2DBL(rb_ary_entry(ary, 1));
  mu    = FIX2INT(rb_ary_entry(ary, 2));
  nu    = FIX2INT(rb_ary_entry(ary, 3));
  return gsl_integration_qaws_table_alloc(alpha, beta, mu, nu);
}

static VALUE rb_gsl_integration_qaws(int argc, VALUE *argv, VALUE obj)
{
  double a, b, epsabs, epsrel;
  double result, abserr;
  size_t limit;
  gsl_function *F = NULL;
  gsl_integration_workspace *w = NULL;
  gsl_integration_qaws_table *t = NULL;
  int status, intervals, itmp, flag = 0, flagt = 0;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 2) rb_raise(rb_eArgError, "too few arguments");
    CHECK_FUNCTION(argv[0]);
    Data_Get_Struct(argv[0], gsl_function, F);
    itmp = 1;
    break;
  default:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments");
    Data_Get_Struct(obj, gsl_function, F);
    itmp = 0;
    break;
  }
  itmp = get_a_b(argc, argv, itmp, &a, &b);

  if (TYPE(argv[itmp]) == T_ARRAY) {
    flagt = 1;
    t = make_qaws_table(argv[itmp]);
  } else {
    flagt = 0;
    if (!rb_obj_is_kind_of(argv[itmp], cgsl_integration_qaws_table))
  rb_raise(rb_eTypeError, "Integration::QAWS_Table expected");

    Data_Get_Struct(argv[itmp], gsl_integration_qaws_table, t);
  }
  flag = get_epsabs_epsrel_limit_workspace(argc, argv, itmp+1, &epsabs, &epsrel,
             &limit, &w);
  status = gsl_integration_qaws(F, a, b, t, epsabs, epsrel, limit, w, &result, &abserr);
  intervals = w->size;
  if (flag == 1) gsl_integration_workspace_free(w);
  if (flagt == 1) gsl_integration_qaws_table_free(t);

  return rb_ary_new3(4, rb_float_new(result), rb_float_new(abserr), INT2FIX(intervals),
         INT2FIX(status));
}          

static gsl_integration_qawo_table* make_qawo_table(VALUE ary);

static VALUE rb_gsl_integration_qawo_table_alloc(int argc, VALUE *argv,
             VALUE klass)
{
  gsl_integration_qawo_table *t = NULL;
  double omega, L;
  enum gsl_integration_qawo_enum sine;
  size_t n;
  if (argc != 1 && argc != 4) 
    rb_raise(rb_eArgError, "wrong nubmer of arguments (%d for 1 or 4)", argc);

  if (TYPE(argv[0]) == T_ARRAY) {
    omega = NUM2DBL(rb_ary_entry(argv[0], 0));
    L = NUM2DBL(rb_ary_entry(argv[0], 1));
    sine = FIX2INT(rb_ary_entry(argv[0], 2));
    n = FIX2INT(rb_ary_entry(argv[0], 3));
  } else {
    omega = NUM2DBL(argv[0]);
    L = NUM2DBL(argv[1]);
    sine = FIX2INT(argv[2]);
    n = FIX2INT(argv[3]);
  }

  t = gsl_integration_qawo_table_alloc(omega, L, sine, n);
               
  return Data_Wrap_Struct(klass, 0, gsl_integration_qawo_table_free, t);
}

static VALUE rb_gsl_integration_qawo_table_to_a(VALUE obj)
{
  gsl_integration_qawo_table *t = NULL;
  VALUE ary;
  Data_Get_Struct(obj, gsl_integration_qawo_table, t);
  ary = rb_ary_new2(4);
  rb_ary_store(ary, 0, rb_float_new(t->omega));
  rb_ary_store(ary, 1, rb_float_new(t->L));
  rb_ary_store(ary, 2, INT2FIX(t->sine));
  rb_ary_store(ary, 3, INT2FIX(t->n));
  return ary;
}

static VALUE rb_gsl_ary_to_integration_qawo_table(VALUE ary)
{
  gsl_integration_qawo_table *t = NULL;
  t = make_qawo_table(ary);
  return Data_Wrap_Struct(cgsl_integration_qawo_table, 
        0, gsl_integration_qawo_table_free, t);
}

static gsl_integration_qawo_table* make_qawo_table(VALUE ary)
{
  double omega, L;
  enum gsl_integration_qawo_enum sine;
  size_t n;
  omega = NUM2DBL(rb_ary_entry(ary, 0));
  L     = NUM2DBL(rb_ary_entry(ary, 1));
  sine  = FIX2INT(rb_ary_entry(ary, 2));
  n     = FIX2INT(rb_ary_entry(ary, 3));
  return gsl_integration_qawo_table_alloc(omega, L, sine, n);
}

static VALUE rb_gsl_integration_qawo_table_set(int argc, VALUE *argv, VALUE obj)
{
  gsl_integration_qawo_table *t = NULL;
  double omega, L;
  enum gsl_integration_qawo_enum sine;
  int type;
  if (argc != 1 && argc != 3)
    rb_raise(rb_eArgError, "wrong number of argument (%d for 1 or 3)", argc);
  type = TYPE(argv[0]);
  Data_Get_Struct(obj, gsl_integration_qawo_table, t);
  if (type == T_FIXNUM || type == T_BIGNUM || type == T_FLOAT) {
    omega = NUM2DBL(argv[0]);
    L     = NUM2DBL(argv[1]);
    sine  = FIX2INT(argv[2]);
  } else if (type == T_ARRAY) {
    omega = NUM2DBL(rb_ary_entry(argv[0], 0));
    L     = NUM2DBL(rb_ary_entry(argv[0], 1));
    sine  = FIX2INT(rb_ary_entry(argv[0], 2));
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s", rb_class2name(CLASS_OF(argv[0])));
  }
  gsl_integration_qawo_table_set(t, omega, L, sine);
  return obj;
}

static VALUE rb_gsl_integration_qawo_table_set_length(VALUE obj, VALUE L)
{
  gsl_integration_qawo_table *t = NULL;
  Need_Float(L);
  Data_Get_Struct(obj, gsl_integration_qawo_table, t);
  gsl_integration_qawo_table_set_length(t, NUM2DBL(L));
  return obj;
}

static int get_qawo_table(VALUE tt, gsl_integration_qawo_table **t);

static VALUE rb_gsl_integration_qawo(int argc, VALUE *argv, VALUE obj)
{
  double a, epsabs, epsrel;
  double result, abserr;
  size_t limit;
  gsl_function *F = NULL;
  gsl_integration_workspace *w = NULL;
  gsl_integration_qawo_table *t = NULL;
  int status, intervals, itmp, flag = 0, flagt = 0;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
  if (argc < 2) rb_raise(rb_eArgError, "too few arguments");
    CHECK_FUNCTION(argv[0]);
    Data_Get_Struct(argv[0], gsl_function, F);
    itmp = 1;
    break;
  default:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments");
    Data_Get_Struct(obj, gsl_function, F);
    itmp = 0;
    break;
  }
  Need_Float(argv[itmp]);
  a = NUM2DBL(argv[itmp]);
  flagt = get_qawo_table(argv[argc-1], &t);
  flag = get_epsabs_epsrel_limit_workspace(argc-1, argv, itmp+1, &epsabs, &epsrel,
             &limit, &w);
  status = gsl_integration_qawo(F, a, epsabs, epsrel, limit, w, t, &result, &abserr);
  intervals = w->size;
  if (flag == 1) gsl_integration_workspace_free(w);
  if (flagt == 1) gsl_integration_qawo_table_free(t);

  return rb_ary_new3(4, rb_float_new(result), rb_float_new(abserr), INT2FIX(intervals),
         INT2FIX(status));
}          

static int get_qawo_table(VALUE tt,
        gsl_integration_qawo_table **t)
{
  int flagt;

  if (TYPE(tt) == T_ARRAY) {
    flagt = 1;
    *t = make_qawo_table(tt);
  } else {
    flagt = 0;
    if (!rb_obj_is_kind_of(tt, cgsl_integration_qawo_table))
      rb_raise(rb_eTypeError, "Integration::QAWO_Table expected");
    Data_Get_Struct(tt, gsl_integration_qawo_table, *t);
  }
  return flagt;
}

static VALUE rb_gsl_integration_qawf(int argc, VALUE *argv, VALUE obj)
{
  double a, epsabs = EPSREL_DEFAULT;
  double result, abserr;
  size_t limit = LIMIT_DEFAULT;
  gsl_function *F = NULL;
  gsl_integration_workspace *w = NULL, *cw = NULL;
  gsl_integration_qawo_table *t = NULL;
  int status, intervals, flag = 0, flagt = 0, itmp;
  VALUE *vtmp;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 2) rb_raise(rb_eArgError, "too few arguments");
    CHECK_FUNCTION(argv[0]);
    Data_Get_Struct(argv[0], gsl_function, F);
    itmp = 1;
    break;
  default:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments");
    Data_Get_Struct(obj, gsl_function, F);
    itmp = 0;
    break;
  }
  Need_Float(argv[itmp]);
  a = NUM2DBL(argv[itmp]);
  itmp += 1;
  if (TYPE(argv[itmp]) == T_FLOAT) {
    epsabs = NUM2DBL(argv[itmp]);
    itmp += 1;
  }
  vtmp = argv + itmp;
  flagt = get_qawo_table(argv[argc-1], &t);    

  switch (argc - 1 - itmp) {
  case 0:
    w = gsl_integration_workspace_alloc(limit);
    cw = gsl_integration_workspace_alloc(limit);
    flag = 1;
    break;
  case 1:
    CHECK_FIXNUM(vtmp[0]);
    limit = FIX2INT(vtmp[0]);
    w = gsl_integration_workspace_alloc(limit);
    cw = gsl_integration_workspace_alloc(limit);
    flag = 1;
    break;
  case 2:
    CHECK_WORKSPACE(vtmp[0]); CHECK_WORKSPACE(vtmp[1]);
    Data_Get_Struct(vtmp[0], gsl_integration_workspace, w);
    Data_Get_Struct(vtmp[1], gsl_integration_workspace, cw);
    flag = 0;
    break;
  case 3:
    CHECK_FIXNUM(vtmp[0]);
    CHECK_WORKSPACE(vtmp[1]); CHECK_WORKSPACE(vtmp[2]);
    limit = FIX2INT(vtmp[0]);
    Data_Get_Struct(vtmp[1], gsl_integration_workspace, w);
    Data_Get_Struct(vtmp[2], gsl_integration_workspace, cw);
    flag = 0;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }

  status = gsl_integration_qawf(F, a, epsabs, limit, w, cw, t, &result, &abserr);
  intervals = w->size;
  if (flag == 1) {
    gsl_integration_workspace_free(w);
    gsl_integration_workspace_free(cw);
  }
  if (flagt == 1) gsl_integration_qawo_table_free(t);

  return rb_ary_new3(4, rb_float_new(result), rb_float_new(abserr), 
         INT2FIX(intervals), INT2FIX(status));
}          


static void rb_gsl_integration_define_symbols(VALUE module)
{
  rb_define_const(module, "GAUSS15", INT2FIX(GSL_INTEG_GAUSS15));
  rb_define_const(module, "GAUSS21", INT2FIX(GSL_INTEG_GAUSS21));
  rb_define_const(module, "GAUSS31", INT2FIX(GSL_INTEG_GAUSS31));
  rb_define_const(module, "GAUSS41", INT2FIX(GSL_INTEG_GAUSS41));
  rb_define_const(module, "GAUSS51", INT2FIX(GSL_INTEG_GAUSS51));
  rb_define_const(module, "GAUSS61", INT2FIX(GSL_INTEG_GAUSS61));
  rb_define_const(module, "COSINE", INT2FIX(GSL_INTEG_COSINE));
  rb_define_const(module, "SINE", INT2FIX(GSL_INTEG_SINE));
}

static VALUE rb_gsl_integration_workspace_alloc(int argc, VALUE *argv,
            VALUE klass)
{
  size_t limit;
  if (argc == 1) limit = FIX2INT(argv[0]);
  else limit = LIMIT_DEFAULT;
  return Data_Wrap_Struct(klass, 0, 
        gsl_integration_workspace_free, 
        gsl_integration_workspace_alloc(limit));
}

static VALUE rb_gsl_integration_workspace_limit(VALUE obj)
{
  gsl_integration_workspace *w = NULL;
  Data_Get_Struct(obj, gsl_integration_workspace, w);
  return INT2FIX(w->limit);
}

static VALUE rb_gsl_integration_workspace_size(VALUE obj)
{
  gsl_integration_workspace *w = NULL;
  Data_Get_Struct(obj, gsl_integration_workspace, w);
  return INT2FIX(w->size);
}

static VALUE rb_gsl_integration_workspace_nrmax(VALUE obj)
{
  gsl_integration_workspace *w = NULL;
  Data_Get_Struct(obj, gsl_integration_workspace, w);
  return INT2FIX(w->nrmax);
}

static VALUE rb_gsl_integration_workspace_i(VALUE obj)
{
  gsl_integration_workspace *w = NULL;
  Data_Get_Struct(obj, gsl_integration_workspace, w);
  return INT2FIX(w->i);
}

static VALUE rb_gsl_integration_workspace_maximum_level(VALUE obj)
{
  gsl_integration_workspace *w = NULL;
  Data_Get_Struct(obj, gsl_integration_workspace, w);
  return INT2FIX(w->maximum_level);
}

static VALUE rb_gsl_integration_workspace_to_a(VALUE obj)
{
  gsl_integration_workspace *w = NULL;
  Data_Get_Struct(obj, gsl_integration_workspace, w);
  return rb_ary_new3(5, INT2FIX(w->limit), INT2FIX(w->size), INT2FIX(w->nrmax),
         INT2FIX(w->i), INT2FIX(w->maximum_level));
}

static VALUE rb_gsl_integration_workspace_alist(VALUE obj)
{
  gsl_integration_workspace *w = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_integration_workspace, w);
  v = rb_gsl_make_vector_view(w->alist, w->limit, 1);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, free, v);
}

static VALUE rb_gsl_integration_workspace_blist(VALUE obj)
{
  gsl_integration_workspace *w = NULL;
  gsl_vector_view *v = NULL;

  Data_Get_Struct(obj, gsl_integration_workspace, w);
  v = rb_gsl_make_vector_view(w->blist, w->limit, 1);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, free, v);
}

static VALUE rb_gsl_integration_workspace_rlist(VALUE obj)
{
  gsl_integration_workspace *w = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_integration_workspace, w);
  v = rb_gsl_make_vector_view(w->rlist, w->limit, 1);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, free, v);
}

static VALUE rb_gsl_integration_workspace_elist(VALUE obj)
{
  gsl_integration_workspace *w = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_integration_workspace, w);
  v = rb_gsl_make_vector_view(w->elist, w->limit, 1);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, free, v);
}

#ifdef GSL_1_14_LATER
static VALUE rb_gsl_integration_glfixed_table_alloc(VALUE klass, VALUE n)
{
  gsl_integration_glfixed_table *t;
  t = gsl_integration_glfixed_table_alloc(FIX2INT(n));
  return Data_Wrap_Struct(cgsl_integration_glfixed_table, 0, gsl_integration_glfixed_table_free, t);
}

static VALUE rb_gsl_integration_glfixed(VALUE obj, VALUE aa, VALUE bb, VALUE tt)
{
  gsl_function *f;
  double a, b;
  gsl_integration_glfixed_table *t;
  double res;
  if (!rb_obj_is_kind_of(tt, cgsl_integration_glfixed_table)) {
    rb_raise(rb_eTypeError, "Wrong arugment type (%s for GSL::Integration::Glfixed_table)",
       rb_class2name(CLASS_OF(tt)));
  }
  Data_Get_Struct(tt, gsl_integration_glfixed_table, t);
  a = NUM2DBL(aa);
  b = NUM2DBL(bb);
  Data_Get_Struct(obj, gsl_function, f);
  res = gsl_integration_glfixed(f, a, b, t);
  return rb_float_new(res);
}
#endif

void Init_gsl_integration(VALUE module)
{
  VALUE mgsl_integ;

  mgsl_integ = rb_define_module_under(module, "Integration");
  rb_gsl_integration_define_symbols(mgsl_integ);

  rb_define_method(cgsl_function, "integration_qng", rb_gsl_integration_qng, -1);
  rb_define_method(cgsl_function, "integration_qag", rb_gsl_integration_qag, -1);
  rb_define_method(cgsl_function, "integration_qags", rb_gsl_integration_qags, -1);
  rb_define_method(cgsl_function, "integration_qagp", rb_gsl_integration_qagp, -1);
  rb_define_method(cgsl_function, "integration_qagi", rb_gsl_integration_qagi, -1);
  rb_define_method(cgsl_function, "integration_qagiu", rb_gsl_integration_qagiu, -1);
  rb_define_method(cgsl_function, "integration_qagil", rb_gsl_integration_qagil, -1);
  rb_define_method(cgsl_function, "integration_qawc", rb_gsl_integration_qawc, -1);
  rb_define_alias(cgsl_function, "qng", "integration_qng");
  rb_define_alias(cgsl_function, "qag", "integration_qag");
  rb_define_alias(cgsl_function, "qags", "integration_qags");
  rb_define_alias(cgsl_function, "qagp", "integration_qagp");
  rb_define_alias(cgsl_function, "qagi", "integration_qagi");
  rb_define_alias(cgsl_function, "qagiu", "integration_qagiu");
  rb_define_alias(cgsl_function, "qagil", "integration_qagil");
  rb_define_alias(cgsl_function, "qawc", "integration_qawc");

  cgsl_integration_qaws_table = rb_define_class_under(mgsl_integ, "QAWS_Table", 
                  cGSL_Object);
  rb_define_singleton_method(cgsl_integration_qaws_table, "alloc",
           rb_gsl_integration_qaws_table_alloc, -1);
  /*  rb_define_singleton_method(cgsl_integration_qaws_table, "new",
      rb_gsl_integration_qaws_table_alloc, -1);*/
  rb_define_method(cgsl_integration_qaws_table, "to_a", 
       rb_gsl_integration_qaws_table_to_a, 0);
  rb_define_method(cgsl_integration_qaws_table, "set", 
       rb_gsl_integration_qaws_table_set, -1);
  rb_define_method(rb_cArray, "to_gsl_integration_qaws_table", 
       rb_gsl_ary_to_integration_qaws_table, 0);
  rb_define_alias(rb_cArray, "to_qaws_table", "to_gsl_integration_qaws_table");
  rb_define_method(cgsl_function, "integration_qaws", rb_gsl_integration_qaws, -1);
  rb_define_alias(cgsl_function, "qaws", "integration_qaws");

  cgsl_integration_qawo_table = rb_define_class_under(mgsl_integ, "QAWO_Table", 
                  cGSL_Object);
  rb_define_singleton_method(cgsl_integration_qawo_table, "alloc",
           rb_gsl_integration_qawo_table_alloc, -1);
  /*  rb_define_singleton_method(cgsl_integration_qawo_table, "new",
      rb_gsl_integration_qawo_table_alloc, -1);*/
  rb_define_method(cgsl_integration_qawo_table, "to_a", 
       rb_gsl_integration_qawo_table_to_a, 0);
  rb_define_method(rb_cArray, "to_gsl_integration_qawo_table", 
       rb_gsl_ary_to_integration_qawo_table, 0);
  rb_define_method(cgsl_integration_qawo_table, "set", 
       rb_gsl_integration_qawo_table_set, -1);
  rb_define_method(cgsl_integration_qawo_table, "set_length", 
       rb_gsl_integration_qawo_table_set_length, 1);
  rb_define_method(cgsl_function, "integration_qawo", rb_gsl_integration_qawo, -1);
  rb_define_method(cgsl_function, "integration_qawf", rb_gsl_integration_qawf, -1);
  rb_define_alias(cgsl_function, "qawo", "integration_qawo");
  rb_define_alias(cgsl_function, "qawf", "integration_qawf");

  cgsl_integration_workspace = rb_define_class_under(mgsl_integ, 
                 "Workspace", cGSL_Object);

  /*  rb_define_singleton_method(cgsl_integration_workspace, "new",
      rb_gsl_integration_workspace_alloc, -1);*/
  rb_define_singleton_method(cgsl_integration_workspace, "alloc",
           rb_gsl_integration_workspace_alloc, -1);

  rb_define_method(cgsl_integration_workspace, "limit", 
       rb_gsl_integration_workspace_limit, 0);
  rb_define_method(cgsl_integration_workspace, "size", 
       rb_gsl_integration_workspace_size, 0);
  rb_define_method(cgsl_integration_workspace, "nrmax", 
       rb_gsl_integration_workspace_nrmax, 0);
  rb_define_method(cgsl_integration_workspace, "i", 
       rb_gsl_integration_workspace_i, 0);
  rb_define_method(cgsl_integration_workspace, "maximum_level", 
       rb_gsl_integration_workspace_maximum_level, 0);
  rb_define_method(cgsl_integration_workspace, "to_a", 
       rb_gsl_integration_workspace_to_a, 0);
  rb_define_method(cgsl_integration_workspace, "alist", 
       rb_gsl_integration_workspace_alist, 0);
  rb_define_method(cgsl_integration_workspace, "blist", 
       rb_gsl_integration_workspace_blist, 0);
  rb_define_method(cgsl_integration_workspace, "rlist", 
       rb_gsl_integration_workspace_rlist, 0);
  rb_define_method(cgsl_integration_workspace, "elist", 
       rb_gsl_integration_workspace_elist, 0);

  /*****/
  rb_define_module_function(mgsl_integ, "qng", rb_gsl_integration_qng, -1);
  rb_define_module_function(mgsl_integ, "qag", rb_gsl_integration_qag, -1);
  rb_define_module_function(mgsl_integ, "qags", rb_gsl_integration_qags, -1);
  rb_define_module_function(mgsl_integ, "qagp", rb_gsl_integration_qagp, -1);
  rb_define_module_function(mgsl_integ, "qagi", rb_gsl_integration_qagi, -1);
  rb_define_module_function(mgsl_integ, "qagiu", rb_gsl_integration_qagiu, -1);
  rb_define_module_function(mgsl_integ, "qagil", rb_gsl_integration_qagil, -1);
  rb_define_module_function(mgsl_integ, "qawc", rb_gsl_integration_qawc, -1);
  rb_define_module_function(mgsl_integ, "qaws", rb_gsl_integration_qaws, -1);
  rb_define_module_function(mgsl_integ, "qawo", rb_gsl_integration_qawo, -1);
  rb_define_module_function(mgsl_integ, "qawf", rb_gsl_integration_qawf, -1);

#ifdef GSL_1_14_LATER
  cgsl_integration_glfixed_table = rb_define_class_under(mgsl_integ, "Glfixed_table", cGSL_Object);
  rb_define_singleton_method(cgsl_integration_glfixed_table, "alloc",
           rb_gsl_integration_glfixed_table_alloc, 1);
  rb_define_method(cgsl_function, "glfixed", rb_gsl_integration_glfixed, 3);
#endif

}

#undef EPSABS_DEFAULT
#undef EPSREL_DEFAULT
#undef LIMIT_DEFAULT
#undef KEY_DEFAULT

#ifdef CHECK_WORKSPACE
#undef CHECK_WORKSPACE
#endif

