/*
  multiroots.c
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
#include <gsl/gsl_multiroots.h>

#ifndef CHECK_MULTIROOT_FUNCTION
#define CHECK_MULTIROOT_FUNCTION(x) if(CLASS_OF(x)!=cgsl_multiroot_function)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (GSL::MultiRoot::Function expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

#ifndef CHECK_MULTIROOT_FUNCTION_FDF
#define CHECK_MULTIROOT_FUNCTION_FDF(x) if(CLASS_OF(x)!=cgsl_multiroot_function_fdf)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (GSL::MultiRoot::Function_fdf expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

static VALUE cgsl_multiroot_function;
static VALUE cgsl_multiroot_function_fdf;

enum {
  GSL_MULTIROOT_FDFSOLVER_HYBRIDSJ,
  GSL_MULTIROOT_FDFSOLVER_HYBRIDJ,
  GSL_MULTIROOT_FDFSOLVER_NEWTON,
  GSL_MULTIROOT_FDFSOLVER_GNEWTON,
  GSL_MULTIROOT_FSOLVER_HYBRIDS,
  GSL_MULTIROOT_FSOLVER_HYBRID,
  GSL_MULTIROOT_FSOLVER_DNEWTON,
  GSL_MULTIROOT_FSOLVER_BROYDEN,
};

static void gsl_multiroot_function_fdf_mark(gsl_multiroot_function_fdf *f);
static void gsl_multiroot_function_mark(gsl_multiroot_function *f);
static void gsl_multiroot_function_free(gsl_multiroot_function *f);
static int rb_gsl_multiroot_function_f(const gsl_vector *x, void *p, gsl_vector *f);
static void set_function(int i, VALUE *argv, gsl_multiroot_function *F);

static void gsl_multiroot_function_fdf_free(gsl_multiroot_function_fdf *f);
static int rb_gsl_multiroot_function_fdf_f(const gsl_vector *x, void *p,
             gsl_vector *f);
static int rb_gsl_multiroot_function_fdf_df(const gsl_vector *x, void *p,
              gsl_matrix *J);
static int rb_gsl_multiroot_function_fdf_fdf(const gsl_vector *x, void *p,
               gsl_vector *f, gsl_matrix *J);
static void set_function_fdf(int i, VALUE *argv, gsl_multiroot_function_fdf *F);
static const gsl_multiroot_fsolver_type* get_fsolver_type(VALUE t);
static const gsl_multiroot_fdfsolver_type* get_fdfsolver_type(VALUE t);

static VALUE rb_gsl_multiroot_function_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_multiroot_function *F = NULL;
  VALUE ary;
  size_t i;
  F = ALLOC(gsl_multiroot_function);
  F->f = &rb_gsl_multiroot_function_f;
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
    break;
  }
  return Data_Wrap_Struct(klass, gsl_multiroot_function_mark, gsl_multiroot_function_free, F);
}

static void gsl_multiroot_function_free(gsl_multiroot_function *f)
{
  free((gsl_multiroot_function *) f);
}

static void gsl_multiroot_function_mark(gsl_multiroot_function *f)
{
  size_t i;
  rb_gc_mark((VALUE) f->params);
  //  for (i = 0; i < RARRAY(f->params)->len; i++)
  for (i = 0; (int) i < RARRAY_LEN(f->params); i++)
    rb_gc_mark(rb_ary_entry((VALUE) f->params, i));
}

static int rb_gsl_multiroot_function_f(const gsl_vector *x, void *p, gsl_vector *f)
{
  VALUE vx, vf;
  VALUE vp, proc;
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  vf = Data_Wrap_Struct(cgsl_vector, 0, NULL, f);
  proc = rb_ary_entry((VALUE) p, 0);
  vp = rb_ary_entry((VALUE) p, 1);
  if (NIL_P(vp)) rb_funcall(proc, RBGSL_ID_call, 2, vx, vf);
  else rb_funcall(proc, RBGSL_ID_call, 3, vx, vp, vf);
  return GSL_SUCCESS;
}

static VALUE rb_gsl_multiroot_function_eval(VALUE obj, VALUE vx)
{
  gsl_multiroot_function *F = NULL;
  gsl_vector *f = NULL;
  VALUE vp, proc, vf, ary;
  Data_Get_Struct(obj, gsl_multiroot_function, F);
  ary = (VALUE) F->params;
  f = gsl_vector_alloc(F->n);
  vf = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, f);
  proc = rb_ary_entry(ary, 0);
  vp = rb_ary_entry(ary, 1);
  if (NIL_P(vp)) rb_funcall(proc, RBGSL_ID_call, 2, vx, vf);
  else rb_funcall(proc, RBGSL_ID_call, 3, vx, vp, vf);
  return vf;
}

static void set_function(int i, VALUE *argv, gsl_multiroot_function *F)
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

static VALUE rb_gsl_multiroot_function_set_f(int argc, VALUE *argv, VALUE obj)
{
  gsl_multiroot_function *F = NULL;
  VALUE ary;
  size_t i;
  Data_Get_Struct(obj, gsl_multiroot_function, F);
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
    break;
  }
  return obj;
}

static VALUE rb_gsl_multiroot_function_set_params(int argc, VALUE *argv, VALUE obj)
{
  gsl_multiroot_function *F = NULL;
  VALUE ary, ary2;
  size_t i;
  if (argc == 0) return obj;
  Data_Get_Struct(obj, gsl_multiroot_function, F);
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

static VALUE rb_gsl_multiroot_function_params(VALUE obj)
{
  gsl_multiroot_function *F = NULL;
  Data_Get_Struct(obj, gsl_multiroot_function, F);
  return rb_ary_entry((VALUE) F->params, 1);
}

static VALUE rb_gsl_multiroot_function_n(VALUE obj)
{
  gsl_multiroot_function *F = NULL;
  Data_Get_Struct(obj, gsl_multiroot_function, F);
  return INT2FIX(F->n);
}

/*** multiroot_function_fdf ***/
static void set_function_fdf(int argc, VALUE *argv, gsl_multiroot_function_fdf *F);
static VALUE rb_gsl_multiroot_function_fdf_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_multiroot_function_fdf *F = NULL;
  VALUE ary;
  F = ALLOC(gsl_multiroot_function_fdf);
  F->f = &rb_gsl_multiroot_function_fdf_f;
  F->df = &rb_gsl_multiroot_function_fdf_df;
  F->fdf = &rb_gsl_multiroot_function_fdf_fdf;
  ary = rb_ary_new2(4);
  /*  (VALUE) F->params = ary;*/
  F->params = (void *) ary;
  rb_ary_store(ary, 2, Qnil);
  rb_ary_store(ary, 3, Qnil);
  set_function_fdf(argc, argv, F);
  return Data_Wrap_Struct(klass, gsl_multiroot_function_fdf_mark, gsl_multiroot_function_fdf_free, F);
}

static void gsl_multiroot_function_fdf_free(gsl_multiroot_function_fdf *f)
{
  free((gsl_multiroot_function_fdf *) f);
}

static void gsl_multiroot_function_fdf_mark(gsl_multiroot_function_fdf *f)
{
  size_t i;
  rb_gc_mark((VALUE) f->params);
  //  for (i = 0; i < RARRAY(f->params)->len; i++)
  for (i = 0; (int) i < RARRAY_LEN(f->params); i++)
    rb_gc_mark(rb_ary_entry((VALUE) f->params, i));
}

static void set_function_fdf(int argc, VALUE *argv, gsl_multiroot_function_fdf *F)
{
  VALUE ary;
  if (F->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  rb_ary_store(ary, 2, Qnil);
  rb_ary_store(ary, 3, Qnil);
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
    rb_ary_store(ary, 0, argv[0]);
    rb_ary_store(ary, 1, argv[1]);
    if (TYPE(argv[2]) == T_FIXNUM) {
      F->n = FIX2INT(argv[2]);
      rb_ary_store(ary, 2, Qnil);
      rb_ary_store(ary, 3, argv[3]);
    } else {
      rb_ary_store(ary, 2, argv[2]);
      F->n = FIX2INT(argv[3]);
      rb_ary_store(ary, 3, Qnil);
    }
    break;
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
    rb_ary_store(ary, 3, argv[4]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (1, 3, or 4)");
    break;
  }
}

static VALUE rb_gsl_multiroot_function_fdf_set_params(int argc, VALUE *argv, VALUE obj)
{
  gsl_multiroot_function_fdf *F = NULL;
  VALUE ary, ary2;
  size_t i;
  if (argc == 0) return obj;
  Data_Get_Struct(obj, gsl_multiroot_function_fdf, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  if (argc == 1) rb_ary_store(ary, 3, argv[0]);
  else {
    ary2 = rb_ary_new2(argc);
    for (i = 0; (int) i < argc; i++) rb_ary_store(ary2, i, argv[i]);
    rb_ary_store(ary, 3, ary2);
  }
  return obj;
}

static VALUE rb_gsl_multiroot_function_fdf_set(int argc, VALUE *argv, VALUE obj)
{
  gsl_multiroot_function_fdf *F = NULL;
  Data_Get_Struct(obj, gsl_multiroot_function_fdf, F);
  set_function_fdf(argc, argv, F);
  return obj;
}

static int rb_gsl_multiroot_function_fdf_f(const gsl_vector *x, void *p,
             gsl_vector *f)
{
  VALUE vx, vf, ary;
  VALUE proc, vp;
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  vf = Data_Wrap_Struct(cgsl_vector, 0, NULL, f);
  ary = (VALUE) p;
  proc = rb_ary_entry(ary, 0);
  vp = rb_ary_entry(ary, 3);
  if (NIL_P(vp)) rb_funcall(proc, RBGSL_ID_call, 2, vx, vf);
  else rb_funcall(proc, RBGSL_ID_call, 3, vx, vp, vf);
  return GSL_SUCCESS;
}

static int rb_gsl_multiroot_function_fdf_df(const gsl_vector *x, void *p,
              gsl_matrix *J)
{
  VALUE vx, vJ, ary;
  VALUE proc, vp;
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  vJ = Data_Wrap_Struct(cgsl_matrix, 0, NULL, J);
  ary = (VALUE) p;
  proc = rb_ary_entry(ary, 1);
  vp = rb_ary_entry(ary, 3);
  if (NIL_P(vp)) rb_funcall(proc, RBGSL_ID_call, 2, vx, vJ);
  else rb_funcall(proc, RBGSL_ID_call, 3, vx, vp, vJ);
  return GSL_SUCCESS;
}

static int rb_gsl_multiroot_function_fdf_fdf(const gsl_vector *x, void *p,
               gsl_vector *f, gsl_matrix *J)
{
  VALUE vx, vf, vJ, ary;
  VALUE proc_f, proc_df, proc_fdf, vp;
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  vf = Data_Wrap_Struct(cgsl_vector, 0, NULL, f);
  vJ = Data_Wrap_Struct(cgsl_matrix, 0, NULL, J);
  ary = (VALUE) p;
  proc_f = rb_ary_entry(ary, 0);
  proc_df = rb_ary_entry(ary, 1);
  proc_fdf = rb_ary_entry(ary, 2);
  vp = rb_ary_entry(ary, 3);
  if (NIL_P(proc_fdf)) {
    if (NIL_P(vp)) {
      rb_funcall(proc_f, RBGSL_ID_call, 2, vx, vf);
      rb_funcall(proc_df, RBGSL_ID_call, 2, vx, vJ);
    } else {
      rb_funcall(proc_f, RBGSL_ID_call, 3, vx, vp, vf);
      rb_funcall(proc_df, RBGSL_ID_call, 3, vx, vp, vJ);
    }
  } else {
    if (NIL_P(vp)) rb_funcall(proc_fdf, RBGSL_ID_call, 3, vx,  vf, vJ);
    else rb_funcall(proc_fdf, RBGSL_ID_call, 4, vx, vp, vf, vJ);
  }
  return GSL_SUCCESS;
}

static VALUE rb_gsl_multiroot_function_fdf_params(VALUE obj)
{
  gsl_multiroot_function_fdf *F = NULL;
  Data_Get_Struct(obj, gsl_multiroot_function_fdf, F);
  return rb_ary_entry((VALUE) F->params, 3);
}

static VALUE rb_gsl_multiroot_function_fdf_n(VALUE obj)
{
  gsl_multiroot_function_fdf *F = NULL;
  Data_Get_Struct(obj, gsl_multiroot_function_fdf, F);
  return INT2FIX(F->n);
}

/**********/

static void multiroot_define_const(VALUE klass1, VALUE klass2);
static void multiroot_define_const(VALUE klass1, VALUE klass2)
{
  rb_define_const(klass1, "HYBRIDSJ", INT2FIX(GSL_MULTIROOT_FDFSOLVER_HYBRIDSJ));
  rb_define_const(klass1, "HYBRIDJ", INT2FIX(GSL_MULTIROOT_FDFSOLVER_HYBRIDJ));
  rb_define_const(klass1, "NEWTON", INT2FIX(GSL_MULTIROOT_FDFSOLVER_NEWTON));
  rb_define_const(klass1, "GNEWTON", INT2FIX(GSL_MULTIROOT_FDFSOLVER_GNEWTON));

  rb_define_const(klass2, "HYBRIDS", INT2FIX(GSL_MULTIROOT_FSOLVER_HYBRIDS));
  rb_define_const(klass2, "HYBRID", INT2FIX(GSL_MULTIROOT_FSOLVER_HYBRID));
  rb_define_const(klass2, "DNEWTON", INT2FIX(GSL_MULTIROOT_FSOLVER_DNEWTON));
  rb_define_const(klass2, "BROYDEN", INT2FIX(GSL_MULTIROOT_FSOLVER_BROYDEN));
}

#include <string.h>
static const gsl_multiroot_fsolver_type* get_fsolver_type(VALUE t)
{
  char name[32];
  switch (TYPE(t)) {
  case T_STRING:
    strcpy(name,STR2CSTR(t));
    if (str_tail_grep(name, "hybrids") == 0) return gsl_multiroot_fsolver_hybrids;
    else if (str_tail_grep(name, "hybrid") == 0) return gsl_multiroot_fsolver_hybrid;
    else if (str_tail_grep(name, "dnewton") == 0) return gsl_multiroot_fsolver_dnewton;
    else if (str_tail_grep(name, "broyden") == 0) return gsl_multiroot_fsolver_broyden;
    else rb_raise(rb_eTypeError, "%s: unknown algorithm", name);
    break;
  case T_FIXNUM:
    switch (FIX2INT(t)) {
    case GSL_MULTIROOT_FSOLVER_HYBRIDS: return gsl_multiroot_fsolver_hybrids; break;
    case GSL_MULTIROOT_FSOLVER_HYBRID: return gsl_multiroot_fsolver_hybrid; break;
    case GSL_MULTIROOT_FSOLVER_DNEWTON: return gsl_multiroot_fsolver_dnewton; break;
    case GSL_MULTIROOT_FSOLVER_BROYDEN: return gsl_multiroot_fsolver_broyden; break;
    default:
      rb_raise(rb_eTypeError, "%d: unknown algorithm", FIX2INT(t));
      break;
    }
    break;
  default:
    rb_raise(rb_eTypeError, "wrong type argument (Fixnum or String expected)");
    break;
  }
}

static const gsl_multiroot_fdfsolver_type* get_fdfsolver_type(VALUE t)
{
  char name[32];
  switch (TYPE(t)) {
  case T_STRING:
    strcpy(name,STR2CSTR(t));
    if (str_tail_grep(name, "hybridsj") == 0) return gsl_multiroot_fdfsolver_hybridsj;
    else if (str_tail_grep(name, "hybridj") == 0) return gsl_multiroot_fdfsolver_hybridj;
    else if (str_tail_grep(name, "gnewton") == 0) return gsl_multiroot_fdfsolver_gnewton;
    else if (str_tail_grep(name, "newton") == 0) return gsl_multiroot_fdfsolver_newton;
    else rb_raise(rb_eTypeError, "%s: unknown algorithm", name);
    break;
  case T_FIXNUM:
    switch (FIX2INT(t)) {
    case GSL_MULTIROOT_FDFSOLVER_HYBRIDSJ: return gsl_multiroot_fdfsolver_hybridsj; break;
    case GSL_MULTIROOT_FDFSOLVER_HYBRIDJ: return gsl_multiroot_fdfsolver_hybridj; break;
    case GSL_MULTIROOT_FDFSOLVER_NEWTON: return gsl_multiroot_fdfsolver_newton; break;
    case GSL_MULTIROOT_FDFSOLVER_GNEWTON: return gsl_multiroot_fdfsolver_gnewton; break;
    default:
      rb_raise(rb_eTypeError, "%d: unknown algorithm", FIX2INT(t));
      break;
    }
    break;
  default:
    rb_raise(rb_eTypeError, "wrong type argument (Fixnum or String expected)");
    break;
  }
}

static VALUE rb_gsl_multiroot_fsolver_new(VALUE klass, VALUE t, VALUE n)
{
  gsl_multiroot_fsolver *s = NULL;
  const gsl_multiroot_fsolver_type *T;
  CHECK_FIXNUM(n);
  T = get_fsolver_type(t);
  s = gsl_multiroot_fsolver_alloc(T, FIX2INT(n));
  return Data_Wrap_Struct(klass, 0, gsl_multiroot_fsolver_free, s);
}

static VALUE rb_gsl_multiroot_fsolver_set(VALUE obj, VALUE vf, VALUE vx)
{
  gsl_multiroot_fsolver *s = NULL;
  gsl_multiroot_function *f = NULL;
  gsl_vector *x = NULL;
  int flag = 0, status;
  CHECK_MULTIROOT_FUNCTION(vf);
  Data_Get_Struct(obj, gsl_multiroot_fsolver, s);
  Data_Get_Struct(vf, gsl_multiroot_function, f);
  if (TYPE(vx) == T_ARRAY) {
    x = gsl_vector_alloc(s->f->size);
    cvector_set_from_rarray(x, vx);
    flag = 1;
  } else {
    Data_Get_Vector(vx, x);
  }
  status = gsl_multiroot_fsolver_set(s, f, x);
  if (flag == 1) gsl_vector_free(x);
  return INT2FIX(status);
}

static VALUE rb_gsl_multiroot_fsolver_name(VALUE obj)
{
  gsl_multiroot_fsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fsolver, s);
  return rb_str_new2(gsl_multiroot_fsolver_name(s));
}

static VALUE rb_gsl_multiroot_fsolver_iterate(VALUE obj)
{
  gsl_multiroot_fsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fsolver, s);
  return INT2FIX(gsl_multiroot_fsolver_iterate(s));
}

static VALUE rb_gsl_multiroot_fsolver_root(VALUE obj)
{
  gsl_multiroot_fsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fsolver, s);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, gsl_multiroot_fsolver_root(s));
}

static VALUE rb_gsl_multiroot_fsolver_x(VALUE obj)
{
  gsl_multiroot_fsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fsolver, s);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, s->x);
}

static VALUE rb_gsl_multiroot_fsolver_dx(VALUE obj)
{
  gsl_multiroot_fsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fsolver, s);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, s->dx);
}

static VALUE rb_gsl_multiroot_fsolver_f(VALUE obj)
{
  gsl_multiroot_fsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fsolver, s);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, s->f);
}

static VALUE rb_gsl_multiroot_fsolver_test_delta(VALUE obj, VALUE ea, VALUE er)
{
  gsl_multiroot_fsolver *s = NULL;
  Need_Float(ea); Need_Float(er);
  Data_Get_Struct(obj, gsl_multiroot_fsolver, s);
  return INT2FIX(gsl_multiroot_test_delta(s->dx, s->x, NUM2DBL(ea), NUM2DBL(er)));
}

static VALUE rb_gsl_multiroot_fsolver_test_residual(VALUE obj, VALUE ea)
{
  gsl_multiroot_fsolver *s = NULL;
  Need_Float(ea);
  Data_Get_Struct(obj, gsl_multiroot_fsolver, s);
  return INT2FIX(gsl_multiroot_test_residual(s->f, NUM2DBL(ea)));
}

/***/

static VALUE rb_gsl_multiroot_fdfsolver_new(VALUE klass, VALUE t, VALUE n)
{
  gsl_multiroot_fdfsolver *s = NULL;
  const gsl_multiroot_fdfsolver_type *T;
  CHECK_FIXNUM(n);
  T = get_fdfsolver_type(t);
  s = gsl_multiroot_fdfsolver_alloc(T, FIX2INT(n));
  return Data_Wrap_Struct(klass, 0, gsl_multiroot_fdfsolver_free, s);
}

static VALUE rb_gsl_multiroot_fdfsolver_set(VALUE obj, VALUE vf, VALUE vx)
{
  gsl_multiroot_fdfsolver *s = NULL;
  gsl_multiroot_function_fdf *f = NULL;
  gsl_vector *x = NULL;
  int flag = 0, status;
  CHECK_MULTIROOT_FUNCTION_FDF(vf);
  Data_Get_Struct(obj, gsl_multiroot_fdfsolver, s);
  Data_Get_Struct(vf, gsl_multiroot_function_fdf, f);
  if (TYPE(vx) == T_ARRAY) {
    x = gsl_vector_alloc(s->f->size);
    cvector_set_from_rarray(x, vx);
    flag = 1;
  } else {
    Data_Get_Vector(vx, x);
  }
  status = gsl_multiroot_fdfsolver_set(s, f, x);
  if (flag == 0) gsl_vector_free(x);
  return INT2FIX(status);
}

static VALUE rb_gsl_multiroot_fdfsolver_name(VALUE obj)
{
  gsl_multiroot_fdfsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fdfsolver, s);
  return rb_str_new2(gsl_multiroot_fdfsolver_name(s));
}


static VALUE rb_gsl_multiroot_fdfsolver_iterate(VALUE obj)
{
  gsl_multiroot_fdfsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fdfsolver, s);
  return INT2FIX(gsl_multiroot_fdfsolver_iterate(s));
}

static VALUE rb_gsl_multiroot_fdfsolver_root(VALUE obj)
{
  gsl_multiroot_fdfsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fdfsolver, s);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, gsl_multiroot_fdfsolver_root(s));
}

static VALUE rb_gsl_multiroot_fdfsolver_x(VALUE obj)
{
  gsl_multiroot_fdfsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fdfsolver, s);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, s->x);
}

static VALUE rb_gsl_multiroot_fdfsolver_dx(VALUE obj)
{
  gsl_multiroot_fdfsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fdfsolver, s);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, s->dx);
}

static VALUE rb_gsl_multiroot_fdfsolver_f(VALUE obj)
{
  gsl_multiroot_fdfsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fdfsolver, s);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, s->f);
}

static VALUE rb_gsl_multiroot_fdfsolver_J(VALUE obj)
{
  gsl_multiroot_fdfsolver *s = NULL;
  Data_Get_Struct(obj, gsl_multiroot_fdfsolver, s);
  return Data_Wrap_Struct(cgsl_matrix_view_ro, 0, NULL, s->J);
}

static VALUE rb_gsl_multiroot_fdfsolver_test_delta(VALUE obj, VALUE ea, VALUE er)
{
  gsl_multiroot_fdfsolver *s = NULL;
  Need_Float(ea); Need_Float(er);
  Data_Get_Struct(obj, gsl_multiroot_fdfsolver, s);
  return INT2FIX(gsl_multiroot_test_delta(s->dx, s->x, NUM2DBL(ea), NUM2DBL(er)));
}

static VALUE rb_gsl_multiroot_fdfsolver_test_residual(VALUE obj, VALUE ea)
{
  gsl_multiroot_fdfsolver *s = NULL;
  Need_Float(ea);
  Data_Get_Struct(obj, gsl_multiroot_fdfsolver, s);
  return INT2FIX(gsl_multiroot_test_residual(s->f, NUM2DBL(ea)));
}

static VALUE rb_gsl_multiroot_test_delta(VALUE obj, VALUE vdx, VALUE vx,
           VALUE ea, VALUE er)
{
  gsl_vector *dx = NULL, *x = NULL;
  Need_Float(ea); Need_Float(er);
  Data_Get_Struct(vdx, gsl_vector, dx);
  Data_Get_Struct(vx, gsl_vector, x);
  return INT2FIX(gsl_multiroot_test_delta(dx, x, NUM2DBL(ea), NUM2DBL(er)));
}

static VALUE rb_gsl_multiroot_test_residual(VALUE obj, VALUE vf, VALUE ea)
{
  gsl_vector *f = NULL;
  Need_Float(ea);
  Data_Get_Struct(vf, gsl_vector, f);
  return INT2FIX(gsl_multiroot_test_residual(f, NUM2DBL(ea)));
}

static VALUE rb_gsl_multiroot_fsolver_fsolve(int argc, VALUE *argv, VALUE obj)
{
  gsl_multiroot_fsolver *s = NULL;
  int iter = 0, itmp = 0, i, status, max_iter = 1000;
  double eps = 1e-7;
  gsl_vector *xnew = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    Data_Get_Struct(argv[0], gsl_multiroot_fsolver, s);
    itmp = 1;
    break;
  default:
    Data_Get_Struct(obj, gsl_multiroot_fsolver, s);
    itmp = 0;
    break;
  }
  for (i = itmp; i < argc; i++) {
    switch (argv[i]) {
    case T_FIXNUM:
      max_iter = FIX2INT(argv[i]);
      break;
    case T_FLOAT:
      eps = NUM2DBL(argv[i]);
      break;
    default:
      rb_raise(rb_eTypeError, "wrong type of argument %s (Fixnum or Float expected)",
         rb_class2name(CLASS_OF(argv[i])));
      break;
    }
  }

  do {
    iter ++;
    status = gsl_multiroot_fsolver_iterate (s);
    if (status) break;
    status = gsl_multiroot_test_residual(s->f, eps);
  } while (status == GSL_CONTINUE && iter < max_iter);
  xnew = gsl_vector_alloc(s->x->size);
  gsl_vector_memcpy(xnew, gsl_multiroot_fsolver_root(s));
  return rb_ary_new3(3,
         Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, xnew),
         INT2FIX(iter), INT2FIX(status));
}

/* singleton */
static VALUE rb_gsl_multiroot_fdjacobian(int argc, VALUE *argv, VALUE obj)
{
  gsl_multiroot_function *F = NULL, func;
  gsl_multiroot_function_fdf *fdf = NULL;
  gsl_vector *x = NULL, *f = NULL;
  gsl_matrix *J = NULL;
  double eps;
  int status;
  if (argc != 4 && argc != 5)
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 4 or 5)", argc);

  if (rb_obj_is_kind_of(argv[0], cgsl_multiroot_function_fdf)) {
    Data_Get_Struct(argv[0], gsl_multiroot_function_fdf, fdf);
    func.f = fdf->f;
    func.n = fdf->n;
    func.params = fdf->params;
    F = &func;
  } else if (rb_obj_is_kind_of(argv[0], cgsl_multiroot_function)) {
    Data_Get_Struct(argv[0], gsl_multiroot_function, F);
  } else {
    rb_raise(rb_eArgError, "wrong argument type %s (MultiRoot::Function or MultiRoot::Function_fdf expected)", rb_class2name(CLASS_OF(argv[0])));
  }

  Need_Float(argv[3]);
  Data_Get_Vector(argv[1], x);
  Data_Get_Vector(argv[2], f);
  eps = NUM2DBL(argv[3]);
  if (argc == 4) {
    J = gsl_matrix_alloc(F->n, F->n);
    status = gsl_multiroot_fdjacobian(F, x, f, eps, J);
    return rb_ary_new3(2, Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, J),
           INT2FIX(status));
  } else {
    Data_Get_Struct(argv[4], gsl_matrix, J);
    status = gsl_multiroot_fdjacobian(F, x, f, eps, J);
    return rb_ary_new3(2, argv[4], INT2FIX(status));
  }
}

static VALUE rb_gsl_multiroot_function_get_f(VALUE obj)
{
  gsl_multiroot_function_fdf *F = NULL;
  Data_Get_Struct(obj, gsl_multiroot_function_fdf, F);
  return rb_ary_entry(((VALUE) F->params), 0);
}

static VALUE rb_gsl_multiroot_function_fdf_get_f(VALUE obj)
{
  gsl_multiroot_function_fdf *F = NULL;
  Data_Get_Struct(obj, gsl_multiroot_function_fdf, F);
  return rb_ary_entry(((VALUE) F->params), 0);
}

static VALUE rb_gsl_multiroot_function_fdf_get_df(VALUE obj)
{
  gsl_multiroot_function_fdf *F = NULL;
  Data_Get_Struct(obj, gsl_multiroot_function_fdf, F);
  return rb_ary_entry(((VALUE) F->params), 1);
}

static VALUE rb_gsl_multiroot_function_solve(int argc, VALUE *argv, VALUE obj)
{
  gsl_multiroot_function *F = NULL;
  gsl_vector *x0 = NULL, *xnew;
  int flag = 0;
  double epsabs = 1e-7;
  size_t max_iter = 10000, iter = 0, i;
  gsl_multiroot_fsolver_type *T
    = (gsl_multiroot_fsolver_type *) gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver *s = NULL;
  int status;
  if (argc < 1) rb_raise(rb_eArgError, "too few arguments (%d for >= 1)", argc);
  Data_Get_Struct(obj, gsl_multiroot_function, F);
  switch (argc) {
  case 4:
  case 3:
  case 2:
    for (i = 1; (int) i < argc; i++) {
      switch (TYPE(argv[i])) {
      case T_STRING:
  T = (gsl_multiroot_fsolver_type *) get_fsolver_type(argv[i]);
  break;
      case T_FLOAT:
  epsabs = NUM2DBL(argv[i]);
  break;
      case T_FIXNUM:
  max_iter = FIX2INT(argv[i]);
  break;
      }
    }
    /* no break */
  case 1:
    if (TYPE(argv[0]) == T_ARRAY) {
      //      if (RARRAY(argv[0])->len != F->n)
      if (RARRAY_LEN(argv[0]) != (int) F->n)
  rb_raise(rb_eRangeError, "array size are different.");
      x0 = gsl_vector_alloc(F->n);
      for (i = 0; i < x0->size; i++)
  gsl_vector_set(x0, i, NUM2DBL(rb_ary_entry(argv[0], i)));
      flag = 1;
    } else {
      Data_Get_Vector(argv[0], x0);
      flag = 0;
    }
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments (%d for 1 - 4)", argc);
    break;
  }
  s = gsl_multiroot_fsolver_alloc (T, F->n);
  gsl_multiroot_fsolver_set (s, F, x0);
  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate (s);
    if (status) break;
    status = gsl_multiroot_test_residual(s->f, epsabs);
  } while (status == GSL_CONTINUE && iter < max_iter);
  xnew = gsl_vector_alloc(x0->size);
  gsl_vector_memcpy(xnew, s->x);
  gsl_multiroot_fsolver_free (s);
  if (flag == 1) gsl_vector_free(x0);
  return rb_ary_new3(3, Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, xnew),
         INT2FIX(iter), INT2FIX(status));
}

void Init_gsl_multiroot(VALUE module)
{
  VALUE mgsl_multiroot;
  VALUE cgsl_multiroot_fdfsolver, cgsl_multiroot_fsolver;

  mgsl_multiroot = rb_define_module_under(module, "MultiRoot");

  rb_define_singleton_method(mgsl_multiroot, "test_delta",
           rb_gsl_multiroot_test_delta, 4);
  rb_define_singleton_method(mgsl_multiroot, "test_residual",
           rb_gsl_multiroot_test_residual, 2);

  rb_define_singleton_method(mgsl_multiroot, "fdjacobian",
           rb_gsl_multiroot_fdjacobian, -1);

  /* multiroot_function */
  cgsl_multiroot_function = rb_define_class_under(mgsl_multiroot, "Function",
              cgsl_function);
  rb_define_singleton_method(cgsl_multiroot_function, "alloc",
           rb_gsl_multiroot_function_new, -1);
  rb_define_method(cgsl_multiroot_function, "eval", rb_gsl_multiroot_function_eval, 1);
  rb_define_alias(cgsl_multiroot_function, "call", "eval");
  rb_define_method(cgsl_multiroot_function, "set", rb_gsl_multiroot_function_set_f, -1);
  rb_define_method(cgsl_multiroot_function, "set_params", rb_gsl_multiroot_function_set_params, -1);
  rb_define_method(cgsl_multiroot_function, "params", rb_gsl_multiroot_function_params, 0);
  rb_define_method(cgsl_multiroot_function, "n", rb_gsl_multiroot_function_n, 0);
  rb_define_method(cgsl_multiroot_function, "f", rb_gsl_multiroot_function_get_f, 0);

  /* multiroot_function_fdf */
  cgsl_multiroot_function_fdf = rb_define_class_under(mgsl_multiroot, "Function_fdf",
              cgsl_multiroot_function);
  rb_define_singleton_method(cgsl_multiroot_function_fdf, "alloc",
           rb_gsl_multiroot_function_fdf_new, -1);
  rb_define_method(cgsl_multiroot_function_fdf, "set", rb_gsl_multiroot_function_fdf_set, -1);
  rb_define_method(cgsl_multiroot_function_fdf, "set_params", rb_gsl_multiroot_function_fdf_set_params, -1);
  rb_define_method(cgsl_multiroot_function_fdf, "params", rb_gsl_multiroot_function_fdf_params, 0);
  rb_define_method(cgsl_multiroot_function_fdf, "n", rb_gsl_multiroot_function_fdf_n, 0);
  rb_define_method(cgsl_multiroot_function_fdf, "f", rb_gsl_multiroot_function_fdf_get_f, 0);
  rb_define_method(cgsl_multiroot_function_fdf, "df", rb_gsl_multiroot_function_fdf_get_df, 0);

  /* solver */
  cgsl_multiroot_fsolver = rb_define_class_under(mgsl_multiroot, "FSolver", cGSL_Object);
  cgsl_multiroot_fdfsolver = rb_define_class_under(mgsl_multiroot, "FdfSolver", cgsl_multiroot_fsolver);

  rb_define_singleton_method(cgsl_multiroot_fsolver, "alloc",
           rb_gsl_multiroot_fsolver_new, 2);
  rb_define_singleton_method(cgsl_multiroot_fdfsolver, "alloc",
           rb_gsl_multiroot_fdfsolver_new, 2);

  rb_define_method(cgsl_multiroot_fsolver, "set", rb_gsl_multiroot_fsolver_set, 2);
  rb_define_method(cgsl_multiroot_fsolver, "name", rb_gsl_multiroot_fsolver_name, 0);
  rb_define_method(cgsl_multiroot_fsolver, "iterate", rb_gsl_multiroot_fsolver_iterate, 0);
  rb_define_method(cgsl_multiroot_fsolver, "root", rb_gsl_multiroot_fsolver_root, 0);
  rb_define_method(cgsl_multiroot_fsolver, "x", rb_gsl_multiroot_fsolver_x, 0);
  rb_define_method(cgsl_multiroot_fsolver, "dx", rb_gsl_multiroot_fsolver_dx, 0);
  rb_define_method(cgsl_multiroot_fsolver, "f", rb_gsl_multiroot_fsolver_f, 0);
  rb_define_method(cgsl_multiroot_fsolver, "test_delta", rb_gsl_multiroot_fsolver_test_delta, 2);
  rb_define_method(cgsl_multiroot_fsolver, "test_residual", rb_gsl_multiroot_fsolver_test_residual, 1);

  rb_define_method(cgsl_multiroot_fdfsolver, "set", rb_gsl_multiroot_fdfsolver_set, 2);
  rb_define_method(cgsl_multiroot_fdfsolver, "name", rb_gsl_multiroot_fdfsolver_name, 0);
  rb_define_method(cgsl_multiroot_fdfsolver, "iterate", rb_gsl_multiroot_fdfsolver_iterate, 0);
  rb_define_method(cgsl_multiroot_fdfsolver, "root", rb_gsl_multiroot_fdfsolver_root, 0);
  rb_define_method(cgsl_multiroot_fdfsolver, "x", rb_gsl_multiroot_fdfsolver_x, 0);
  rb_define_method(cgsl_multiroot_fdfsolver, "dx", rb_gsl_multiroot_fdfsolver_dx, 0);
  rb_define_method(cgsl_multiroot_fdfsolver, "f", rb_gsl_multiroot_fdfsolver_f, 0);
  rb_define_method(cgsl_multiroot_fdfsolver, "J", rb_gsl_multiroot_fdfsolver_J, 0);
  rb_define_alias(cgsl_multiroot_fdfsolver, "jac", "J");
  rb_define_alias(cgsl_multiroot_fdfsolver, "jacobian", "J");

  rb_define_method(cgsl_multiroot_fdfsolver, "test_delta", rb_gsl_multiroot_fdfsolver_test_delta, 2);
  rb_define_method(cgsl_multiroot_fdfsolver, "test_residual", rb_gsl_multiroot_fdfsolver_test_residual, 1);


  multiroot_define_const(cgsl_multiroot_fdfsolver, cgsl_multiroot_fsolver);

  rb_define_method(cgsl_multiroot_fsolver, "fsolve", rb_gsl_multiroot_fsolver_fsolve, -1);
  rb_define_alias(cgsl_multiroot_fsolver, "solve", "fsolve");

  rb_define_singleton_method(cgsl_multiroot_fsolver, "fsolve", rb_gsl_multiroot_fsolver_fsolve, -1);
  rb_define_singleton_method(cgsl_multiroot_fsolver, "solve", rb_gsl_multiroot_fsolver_fsolve, -1);

  /*****/
  rb_define_method(cgsl_multiroot_function, "solve", rb_gsl_multiroot_function_solve, -1);
  rb_define_alias(cgsl_multiroot_function, "fsolve", "solve");
}
#ifdef CHECK_MULTIROOT_FUNCTION
#undef CHECK_MULTIROOT_FUNCTION
#endif
#ifdef CHECK_MULTIROOT_FUNCTION_FDF
#undef CHECK_MULTIROOT_FUNCTION_FDF
#endif
