/*
  ntuple.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl.h"
#include "include/rb_gsl_array.h"
#include "include/rb_gsl_function.h"
#include "include/rb_gsl_histogram.h"
#include <gsl/gsl_ntuple.h>

static VALUE cgsl_ntuple;
static VALUE cgsl_ntuple_select_fn;
static VALUE cgsl_ntuple_value_fn;

static VALUE rb_gsl_ntuple_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_vector *v = NULL;
  gsl_matrix *m = NULL;
  gsl_ntuple *n = NULL;
  double *data = NULL;
  size_t size;
  switch (argc) {
  case 2:
  case 3:
    if (VECTOR_P(argv[1])) {
      Data_Get_Struct(argv[1], gsl_vector, v);
      data = v->data;
      size = v->size;
    } else if (MATRIX_P(argv[1])) {
      Data_Get_Struct(argv[1], gsl_matrix, m);
      data = m->data;
      size = m->size1*m->size2;
    } else {
      rb_raise(rb_eTypeError, "Vector or Matrix expected");
    }
    if (argc == 3) size = FIX2INT(argv[2]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  n = gsl_ntuple_create(STR2CSTR(argv[0]), data, size*sizeof(double));
  return Data_Wrap_Struct(klass, 0, gsl_ntuple_close, n);
}

VALUE rb_gsl_ntuple_open(int argc, VALUE *argv, VALUE klass)
{
  gsl_vector *v = NULL;
  gsl_matrix *m = NULL;
  gsl_ntuple *n = NULL;
  double *data = NULL;
  size_t size;
  switch (argc) {
  case 2:
  case 3:
    if (VECTOR_P(argv[1])) {
      Data_Get_Struct(argv[1], gsl_vector, v);
      data = v->data;
      size = v->size;
    } else if (MATRIX_P(argv[1])) {
      Data_Get_Struct(argv[1], gsl_matrix, m);
      data = m->data;
      size = m->size1*m->size2;
    } else {
      rb_raise(rb_eTypeError, "Vector or Matrix expected");
    }
    if (argc == 3) size = FIX2INT(argv[2]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  n = gsl_ntuple_open(STR2CSTR(argv[0]), data, size*sizeof(double));
  return Data_Wrap_Struct(klass, 0, gsl_ntuple_close, n);
}

VALUE rb_gsl_ntuple_write(VALUE obj)
{
  gsl_ntuple *n = NULL;
  Data_Get_Struct(obj, gsl_ntuple, n);
  gsl_ntuple_write(n);
  return obj;
}

VALUE rb_gsl_ntuple_bookdata(VALUE obj)
{
  gsl_ntuple *n = NULL;
  Data_Get_Struct(obj, gsl_ntuple, n);
  gsl_ntuple_bookdata(n);
  return obj;
}

VALUE rb_gsl_ntuple_read(VALUE obj)
{
  gsl_ntuple *n = NULL;
  Data_Get_Struct(obj, gsl_ntuple, n);
  gsl_ntuple_read(n);
  return obj;
}

VALUE rb_gsl_ntuple_close(VALUE klass, VALUE obj)
{
  gsl_ntuple *n = NULL;
  Data_Get_Struct(obj, gsl_ntuple, n);
  gsl_ntuple_close(n);
  return Qnil;
}

VALUE rb_gsl_ntuple_size(VALUE klass, VALUE obj)
{
  gsl_ntuple *n = NULL;
  Data_Get_Struct(obj, gsl_ntuple, n);
  return INT2FIX(n->size);
}

VALUE rb_gsl_ntuple_data(VALUE obj)
{
  gsl_ntuple *n = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_ntuple, n);
  v = gsl_vector_view_alloc();
  v->vector.size = n->size;
  v->vector.data = n->ntuple_data;
  return Data_Wrap_Struct(cgsl_vector_view, 0, gsl_vector_view_free, v);
}

/***** select_fn *****/

static gsl_ntuple_select_fn* gsl_ntuple_select_fn_alloc();
static void gsl_ntuple_select_fn_free(gsl_ntuple_select_fn *ptr);
int rb_gsl_ntuple_select_fn_f(void *data, void *p);
static void gsl_ntuple_select_fn_mark(gsl_ntuple_select_fn *ptr);

static gsl_ntuple_select_fn* gsl_ntuple_select_fn_alloc()
{
  gsl_ntuple_select_fn *ptr = NULL;
  ptr = ALLOC(gsl_ntuple_select_fn);
  if (ptr == NULL) rb_raise(rb_eRuntimeError, "malloc failed");
  ptr->function = &rb_gsl_ntuple_select_fn_f;
  /*  (VALUE) ptr->params = rb_ary_new2(3);*/
  ptr->params = (void *) rb_ary_new2(3);
  rb_ary_store((VALUE) ptr->params, 1, Qnil);
  return ptr;
}

static void gsl_ntuple_select_fn_free(gsl_ntuple_select_fn *ptr)
{
  free((gsl_ntuple_select_fn *) ptr);
}

static void gsl_ntuple_select_fn_mark(gsl_ntuple_select_fn *ptr)
{
  rb_gc_mark((VALUE) ptr->params);
}

static VALUE rb_gsl_ntuple_select_fn_set_f(int argc, VALUE *argv, VALUE obj)
{
  gsl_ntuple_select_fn *F = NULL;
  VALUE ary, ary2;
  size_t i;
  Data_Get_Struct(obj, gsl_ntuple_select_fn, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(3);
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

int rb_gsl_ntuple_select_fn_f(void *data, void *p)
{
  VALUE result, ary, proc, params, vv;
  gsl_vector_view vtmp;
  size_t size;
  ary = (VALUE) p;
  proc = rb_ary_entry(ary, 0);
  params = rb_ary_entry(ary, 1);
  size = FIX2INT(rb_ary_entry(ary, 2));
  vtmp.vector.data = (double *) data;
  vtmp.vector.size = size;
  vtmp.vector.stride = 1;
  vv = Data_Wrap_Struct(cgsl_vector_view, 0, NULL, &vtmp);
  if (NIL_P(params)) result = rb_funcall(proc, RBGSL_ID_call, 1, vv);
  else result = rb_funcall(proc, RBGSL_ID_call, 2, vv, params);
  return FIX2INT(result);
}

static VALUE rb_gsl_ntuple_select_fn_params(VALUE obj)
{
  gsl_ntuple_select_fn *F = NULL;
  Data_Get_Struct(obj, gsl_ntuple_select_fn, F);
  return rb_ary_entry((VALUE) F->params, 1);
}

static VALUE rb_gsl_ntuple_select_fn_set_params(int argc, VALUE *argv, VALUE obj)
{
  gsl_ntuple_select_fn *F = NULL;
  VALUE ary, ary2;
  size_t i;
  if (argc == 0) return obj;
  Data_Get_Struct(obj, gsl_ntuple_select_fn, F);
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

static VALUE rb_gsl_ntuple_select_fn_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_ntuple_select_fn *F = NULL;
  VALUE ff;
  F = gsl_ntuple_select_fn_alloc();
  ff = Data_Wrap_Struct(klass, gsl_ntuple_select_fn_mark, gsl_ntuple_select_fn_free, F);
  rb_gsl_ntuple_select_fn_set_f(argc, argv, ff);
  return ff;
}

/***** value_fn *****/
static gsl_ntuple_value_fn* gsl_ntuple_value_fn_alloc();
static void gsl_ntuple_value_fn_free(gsl_ntuple_value_fn *ptr);
static double rb_gsl_ntuple_value_fn_f(void *data, void *p);
static void gsl_ntuple_value_fn_mark(gsl_ntuple_value_fn *ptr);

static gsl_ntuple_value_fn* gsl_ntuple_value_fn_alloc()
{
  gsl_ntuple_value_fn *ptr = NULL;
  ptr = ALLOC(gsl_ntuple_value_fn);
  if (ptr == NULL) rb_raise(rb_eRuntimeError, "malloc failed");
  ptr->function = &rb_gsl_ntuple_value_fn_f;
  /*  (VALUE) ptr->params = rb_ary_new2(3);*/
  ptr->params = (void *) rb_ary_new2(3);
  rb_ary_store((VALUE) ptr->params, 1, Qnil);
  return ptr;
}

static void gsl_ntuple_value_fn_mark(gsl_ntuple_value_fn *ptr)
{
  rb_gc_mark((VALUE) ptr->params);
}

static void gsl_ntuple_value_fn_free(gsl_ntuple_value_fn *ptr)
{
  free((gsl_ntuple_value_fn *) ptr);
}

static VALUE rb_gsl_ntuple_value_fn_set_f(int argc, VALUE *argv, VALUE obj)
{
  gsl_ntuple_value_fn *F = NULL;
  VALUE ary, ary2;
  size_t i;
  Data_Get_Struct(obj, gsl_ntuple_value_fn, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(3);
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

static double rb_gsl_ntuple_value_fn_f(void *data, void *p)
{
  VALUE result, ary, proc, params, vv;
  gsl_vector_view vtmp;
  ary = (VALUE) p;
  proc = rb_ary_entry(ary, 0);
  params = rb_ary_entry(ary, 1);
  vtmp.vector.data = (double *) data;
  vtmp.vector.size = FIX2INT(rb_ary_entry(ary, 2));
  vtmp.vector.stride = 1;
  vv = Data_Wrap_Struct(cgsl_vector_view, 0, NULL, &vtmp);
  if (NIL_P(params)) result = rb_funcall(proc, RBGSL_ID_call, 1, vv);
  else result = rb_funcall(proc, RBGSL_ID_call, 2, vv, params);
  return NUM2DBL(result);
}

static VALUE rb_gsl_ntuple_value_fn_params(VALUE obj)
{
  gsl_ntuple_value_fn *F = NULL;
  Data_Get_Struct(obj, gsl_ntuple_value_fn, F);
  return rb_ary_entry((VALUE) F->params, 1);
}

static VALUE rb_gsl_ntuple_value_fn_set_params(int argc, VALUE *argv, VALUE obj)
{
  gsl_ntuple_value_fn *F = NULL;
  VALUE ary, ary2;
  size_t i;
  if (argc == 0) return obj;
  Data_Get_Struct(obj, gsl_ntuple_value_fn, F);
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


static VALUE rb_gsl_ntuple_value_fn_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_ntuple_value_fn *F = NULL;
  VALUE ff;
  F = gsl_ntuple_value_fn_alloc();
  ff = Data_Wrap_Struct(klass, gsl_ntuple_value_fn_mark, gsl_ntuple_value_fn_free, F);
  rb_gsl_ntuple_value_fn_set_f(argc, argv, ff);
  return ff;
}

/* singleton method */
static VALUE rb_gsl_ntuple_project(VALUE obj, VALUE hh, VALUE nn,
				   VALUE vvfn, VALUE vsfn)
{
  gsl_histogram *h = NULL;
  gsl_ntuple *n = NULL;
  gsl_ntuple_value_fn *vfn = NULL;
  gsl_ntuple_select_fn *sfn = NULL;
  int status;
  size_t size;
  if (!rb_obj_is_kind_of(hh, cgsl_histogram))
    rb_raise(rb_eTypeError, "argument 1: Histogram expected");
  Data_Get_Struct(hh, gsl_histogram, h);
  if (!rb_obj_is_kind_of(nn, cgsl_ntuple))
    rb_raise(rb_eTypeError, "argument 2: Ntuple expected");
  Data_Get_Struct(nn, gsl_ntuple, n);
  if (!rb_obj_is_kind_of(vvfn, cgsl_ntuple_value_fn)) 
    rb_raise(rb_eTypeError, "argument 3: Ntuple::ValueFn expected");
  Data_Get_Struct(vvfn, gsl_ntuple_value_fn, vfn);
  if (!rb_obj_is_kind_of(vsfn, cgsl_ntuple_select_fn)) 
    rb_raise(rb_eTypeError, "argument 4: Ntuple::SelectFn expected");
  Data_Get_Struct(vsfn, gsl_ntuple_select_fn, sfn);

  size = n->size/sizeof(double);
  rb_ary_store((VALUE) vfn->params, 2, INT2FIX(size));
  rb_ary_store((VALUE) sfn->params, 2, INT2FIX(size));
  status = gsl_ntuple_project(h, n, vfn, sfn);
  return INT2FIX(status);
}

/* method */
static VALUE rb_gsl_ntuple_project2(VALUE obj, VALUE hh, VALUE vvfn, VALUE vsfn)
{
  gsl_histogram *h = NULL;
  gsl_ntuple *n = NULL;
  gsl_ntuple_value_fn *vfn = NULL;
  gsl_ntuple_select_fn *sfn = NULL;
  int status;
  size_t size;
  CHECK_HISTOGRAM(hh);
  Data_Get_Struct(obj, gsl_ntuple, n);
  Data_Get_Struct(hh, gsl_histogram, h);
  if (!rb_obj_is_kind_of(vvfn, cgsl_ntuple_value_fn)) 
    rb_raise(rb_eTypeError, "argument 2: Ntuple::ValueFn expected");
  Data_Get_Struct(vvfn, gsl_ntuple_value_fn, vfn);
  if (!rb_obj_is_kind_of(vsfn, cgsl_ntuple_select_fn)) 
    rb_raise(rb_eTypeError, "argument 3: Ntuple::SelectFn expected");
  Data_Get_Struct(vsfn, gsl_ntuple_select_fn, sfn);
  size = n->size/sizeof(double);
  rb_ary_store((VALUE) vfn->params, 2, INT2FIX(size));
  rb_ary_store((VALUE) sfn->params, 2, INT2FIX(size));
  status = gsl_ntuple_project(h, n, vfn, sfn);
  return INT2FIX(status);
}

void Init_gsl_ntuple(VALUE module)
{
  cgsl_ntuple = rb_define_class_under(module, "Ntuple", cGSL_Object);
  cgsl_ntuple_select_fn = rb_define_class_under(cgsl_ntuple, "SelectFn", cGSL_Object);
  cgsl_ntuple_value_fn = rb_define_class_under(cgsl_ntuple, "ValueFn", cGSL_Object);

  rb_define_singleton_method(cgsl_ntuple, "create", rb_gsl_ntuple_new, -1);
  rb_define_singleton_method(cgsl_ntuple, "alloc", rb_gsl_ntuple_new, -1);
  rb_define_singleton_method(cgsl_ntuple, "open", rb_gsl_ntuple_open, -1);
  rb_define_singleton_method(cgsl_ntuple, "close", rb_gsl_ntuple_close, 0);

  rb_define_method(cgsl_ntuple, "size", rb_gsl_ntuple_size, 0);
  rb_define_method(cgsl_ntuple, "write", rb_gsl_ntuple_write, 0);
  rb_define_method(cgsl_ntuple, "bookdata", rb_gsl_ntuple_bookdata, 0);

  rb_define_method(cgsl_ntuple, "read", rb_gsl_ntuple_read, 0);

  rb_define_method(cgsl_ntuple, "data", rb_gsl_ntuple_data, 0);
  rb_define_alias(cgsl_ntuple, "get_data", "data");
  rb_define_alias(cgsl_ntuple, "ntuple_data", "data");

  rb_define_singleton_method(cgsl_ntuple_select_fn, "alloc",
			     rb_gsl_ntuple_select_fn_new, -1);
  rb_define_method(cgsl_ntuple_select_fn, "set",
		   rb_gsl_ntuple_select_fn_set_f, -1);
  rb_define_method(cgsl_ntuple_select_fn, "set_params",
		   rb_gsl_ntuple_select_fn_set_params, -1);
  rb_define_method(cgsl_ntuple_select_fn, "params",
		   rb_gsl_ntuple_select_fn_params, 0);

  rb_define_singleton_method(cgsl_ntuple_value_fn, "alloc",
			     rb_gsl_ntuple_value_fn_new, -1);
  rb_define_method(cgsl_ntuple_value_fn, "set",
		   rb_gsl_ntuple_value_fn_set_f, -1);
  rb_define_method(cgsl_ntuple_value_fn, "set_params",
		   rb_gsl_ntuple_value_fn_set_params, -1);
  rb_define_method(cgsl_ntuple_value_fn, "params",
		   rb_gsl_ntuple_value_fn_params, 0);

  rb_define_singleton_method(cgsl_ntuple, "project", rb_gsl_ntuple_project, 4);
  rb_define_method(cgsl_ntuple, "project", rb_gsl_ntuple_project2, 3);
}
