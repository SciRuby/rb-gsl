/*
  monte.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl.h"
#include "rb_gsl_rng.h"
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#ifndef CHECK_MONTE_FUNCTION
#define CHECK_MONTE_FUNCTION(x) if(!rb_obj_is_kind_of(x,cgsl_monte_function))\
    rb_raise(rb_eTypeError, "wrong type (Function expected)");
#endif

static VALUE cgsl_monte_plain;
static VALUE cgsl_monte_miser;
static VALUE cgsl_monte_vegas;
static VALUE cgsl_monte_function;
#ifdef GSL_1_13_LATER
static VALUE cgsl_monte_miser_params, cgsl_monte_vegas_params;
#endif
EXTERN VALUE cgsl_vector;

enum {
  GSL_MONTE_PLAIN_STATE = 1,
  GSL_MONTE_MISER_STATE = 2,
  GSL_MONTE_VEGAS_STATE = 3,
};

static void gsl_monte_function_mark(gsl_monte_function *f);
static void gsl_monte_function_free(gsl_monte_function *f);
static double rb_gsl_monte_function_f(double *x, size_t dim, void *p);

static VALUE rb_gsl_monte_function_set_f(int argc, VALUE *argv, VALUE obj)
{
  gsl_monte_function *F = NULL;
  VALUE ary, ary2;
  size_t i;
  Data_Get_Struct(obj, gsl_monte_function, F);
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
    if (TYPE(argv[0]) == T_FIXNUM) {
      F->dim = FIX2INT(argv[0]);
    } else {
      rb_ary_store(ary, 0, argv[0]);
    }
    break;
  case 2:
    rb_ary_store(ary, 0, argv[0]);
    F->dim = FIX2INT(argv[1]);
    break;
  default:
    rb_ary_store(ary, 0, argv[0]);
    F->dim = FIX2INT(argv[1]);
    ary2 = rb_ary_new2(argc-2);
    for (i = 2; i < (size_t) argc; i++) rb_ary_store(ary2, i-2, argv[i]);
    rb_ary_store(ary, 1, ary2);
    break;
  }
  if (rb_block_given_p()) rb_ary_store(ary, 0, rb_block_proc());
  return obj;
}

static void gsl_monte_function_free(gsl_monte_function *f)
{
  free((gsl_monte_function *) f);
}

static void gsl_monte_function_mark(gsl_monte_function *f)
{
  rb_gc_mark((VALUE) f->params);
}

static VALUE rb_gsl_monte_function_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_monte_function *f;
  VALUE obj;
  f = ALLOC(gsl_monte_function);
  f->f = &rb_gsl_monte_function_f;
  /*  (VALUE) f->params = rb_ary_new2(2);*/
  f->params = (void *) rb_ary_new2(2);
  rb_ary_store((VALUE) f->params, 1, Qnil);
  obj = Data_Wrap_Struct(klass, gsl_monte_function_mark, gsl_monte_function_free, f);
  rb_gsl_monte_function_set_f(argc, argv, obj);
  return obj;
}			    

static double rb_gsl_monte_function_f(double *x, size_t dim, void *p)
{
  VALUE result, ary, proc, params;
  gsl_vector vtmp;
  VALUE vx;

  vtmp.data = x;
  vtmp.size = dim;
  vtmp.stride = 1;
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, &vtmp);
  ary = (VALUE) p;
  proc = rb_ary_entry(ary, 0);
  params = rb_ary_entry(ary, 1);
  if (NIL_P(params)) result = rb_funcall(proc, RBGSL_ID_call, 2, vx, INT2FIX(dim));
  else result = rb_funcall(proc, RBGSL_ID_call, 3,  vx, INT2FIX(dim), params);
  return NUM2DBL(result);
}


static VALUE rb_gsl_monte_function_eval(VALUE obj, VALUE vx)
{
  gsl_monte_function *F = NULL;
  VALUE result, ary, proc, params;

  Data_Get_Struct(obj, gsl_monte_function, F);
  ary = (VALUE) F->params;
  proc = rb_ary_entry(ary, 0);
  params = rb_ary_entry(ary, 1);
  if (NIL_P(params)) result = rb_funcall(proc, RBGSL_ID_call, 2, vx, INT2FIX(F->dim));
  else result = rb_funcall(proc, RBGSL_ID_call, 3,  vx, INT2FIX(F->dim), params);
  return result;
}

static VALUE rb_gsl_monte_function_proc(VALUE obj)
{
  gsl_monte_function *F = NULL;
  Data_Get_Struct(obj, gsl_monte_function, F);
  return rb_ary_entry((VALUE) F->params, 0);
}

static VALUE rb_gsl_monte_function_params(VALUE obj)
{
  gsl_monte_function *F = NULL;
  Data_Get_Struct(obj, gsl_monte_function, F);
  return rb_ary_entry((VALUE) F->params, 1);
}

static VALUE rb_gsl_monte_function_set_params(int argc, VALUE *argv, VALUE obj)
{
  gsl_monte_function *F = NULL;
  VALUE ary, ary2;
  size_t i;

  if (argc == 0) return obj;
  Data_Get_Struct(obj, gsl_monte_function, F);
  ary = (VALUE) F->params;
  if (argc == 1) {
    rb_ary_store(ary, 1, argv[0]);
  } else {
    ary2 = rb_ary_new2(argc);
    for (i = 0; i < (size_t) argc; i++) rb_ary_store(ary2, i, argv[i]);
    rb_ary_store(ary, 1, ary2);
  }
  return obj;
}

static VALUE rb_gsl_monte_plain_new(VALUE klass, VALUE d)
{
  gsl_monte_plain_state *s;
  size_t dim;
  CHECK_FIXNUM(d);
  dim = FIX2INT(d);
  s = gsl_monte_plain_alloc(dim);
  gsl_monte_plain_init(s);
  return Data_Wrap_Struct(klass, 0, gsl_monte_plain_free, s);
}

static VALUE rb_gsl_monte_plain_init(VALUE obj)
{
  gsl_monte_plain_state *s;
  Data_Get_Struct(obj, gsl_monte_plain_state, s);
  return INT2FIX(gsl_monte_plain_init(s));
}

static VALUE rb_gsl_monte_miser_new(VALUE klass, VALUE d)
{
  gsl_monte_miser_state *s;
  size_t dim;
  CHECK_FIXNUM(d);
  dim = FIX2INT(d);
  s = gsl_monte_miser_alloc(dim);
  gsl_monte_miser_init(s);
  return Data_Wrap_Struct(klass, 0, gsl_monte_miser_free, s);
}

static VALUE rb_gsl_monte_miser_init(VALUE obj)
{
  gsl_monte_miser_state *s;
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  return INT2FIX(gsl_monte_miser_init(s));
}

static VALUE rb_gsl_monte_vegas_new(VALUE klass, VALUE d)
{
  gsl_monte_vegas_state *s;
  size_t dim;
  CHECK_FIXNUM(d);
  dim = FIX2INT(d);
  s = gsl_monte_vegas_alloc(dim);
  gsl_monte_vegas_init(s);
  return Data_Wrap_Struct(klass, 0, gsl_monte_vegas_free, s);
}

static VALUE rb_gsl_monte_vegas_init(VALUE obj)
{
  gsl_monte_vegas_state *s;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  return INT2FIX(gsl_monte_vegas_init(s));
}

static int get_monte_type(VALUE vt);
static VALUE rb_gsl_monte_integrate(int argc, VALUE *argv, VALUE obj)
{
  gsl_monte_function *F = NULL;
  gsl_monte_plain_state *plain = NULL;
  gsl_monte_miser_state *miser = NULL;
  gsl_monte_vegas_state *vegas = NULL;
  gsl_vector *xl = NULL, *xu = NULL;
  gsl_rng *r = NULL;
  size_t dim, calls;
  int itmp = 0, flagr = 0, type;
  double result, abserr;

  if (argc < 4) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 4)", argc);
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (!rb_obj_is_kind_of(argv[0], cgsl_monte_function))
      rb_raise(rb_eTypeError, 
	       "wrong type argument %s (GSL::Monte::Function expected)", 
	       rb_class2name(CLASS_OF(argv[0])));
    Data_Get_Struct(argv[0], gsl_monte_function, F);
    itmp = 1;
    break;
  default:
    Data_Get_Struct(obj, gsl_monte_function, F);
    itmp = 0;
  }

  CHECK_VECTOR(argv[itmp]);
  CHECK_VECTOR(argv[itmp+1]);
  Data_Get_Struct(argv[itmp], gsl_vector, xl);
  Data_Get_Struct(argv[itmp+1], gsl_vector, xu);

  if (argc > itmp+3 && TYPE(argv[itmp+3]) == T_FIXNUM) {
    dim = FIX2INT(argv[itmp+2]);
    calls = FIX2INT(argv[itmp+3]);
  } else {
    dim = F->dim;
    calls = FIX2INT(argv[itmp+2]);
  }

  if (rb_obj_is_kind_of(argv[argc-2], cgsl_rng)) {
    Data_Get_Struct(argv[argc-2], gsl_rng, r);
  } else {
    r = gsl_rng_alloc(gsl_rng_default);
    flagr = 1;
  }

  type = get_monte_type(argv[argc-1]);

  switch (type) {
  case GSL_MONTE_PLAIN_STATE:
  case GSL_MONTE_PLAIN_STATE+100:
    if (type > 100) {
      plain = gsl_monte_plain_alloc(dim);
      gsl_monte_plain_init(plain);
    } else {
      if (!rb_obj_is_kind_of(argv[argc-1], cgsl_monte_plain))
	rb_raise(rb_eTypeError, "wrong argument type %s (Monte::Plain expected)",
		 rb_class2name(CLASS_OF(argv[argc-1])));
      Data_Get_Struct(argv[argc-1], gsl_monte_plain_state, plain);
    }
    gsl_monte_plain_integrate(F, xl->data, xu->data, dim, calls, r, plain, &result, &abserr);
    if (type > 100) gsl_monte_plain_free(plain);
    break;
  case GSL_MONTE_MISER_STATE:
  case GSL_MONTE_MISER_STATE+100:
    if (type > 100) {
      miser = gsl_monte_miser_alloc(dim);
      gsl_monte_miser_init(miser);
    } else {
      if (!rb_obj_is_kind_of(argv[argc-1], cgsl_monte_miser))
	rb_raise(rb_eTypeError, "wrong argument type %s (Monte::Miser expected)",
		 rb_class2name(CLASS_OF(argv[argc-1])));
      Data_Get_Struct(argv[argc-1], gsl_monte_miser_state, miser);
    }
    gsl_monte_miser_integrate(F, xl->data, xu->data, dim, calls, r, miser, &result, &abserr);
    if (type > 100) gsl_monte_miser_free(miser);
    break;
  case GSL_MONTE_VEGAS_STATE:
  case GSL_MONTE_VEGAS_STATE+100:
    if (type > 100) {
      vegas = gsl_monte_vegas_alloc(dim);
      gsl_monte_vegas_init(vegas);
    } else {      if (!rb_obj_is_kind_of(argv[argc-1], cgsl_monte_vegas))
	rb_raise(rb_eTypeError, "wrong argument type %s (Monte::Vegas expected)",
		 rb_class2name(CLASS_OF(argv[argc-1])));
      Data_Get_Struct(argv[argc-1], gsl_monte_vegas_state, vegas);
    }
    gsl_monte_vegas_integrate(F, xl->data, xu->data, dim, calls, r, vegas, &result, &abserr);
    if (type > 100) gsl_monte_vegas_free(vegas);
    break;
  }
  if (flagr == 1) gsl_rng_free(r);
  return rb_ary_new3(2, rb_float_new(result), rb_float_new(abserr));
}

static VALUE rb_gsl_monte_plain_integrate(int argc, VALUE *argv, VALUE obj)
{
  gsl_monte_function *F = NULL;
  gsl_monte_plain_state *plain = NULL;
  gsl_vector *xl = NULL, *xu = NULL;
  gsl_rng *r = NULL;
  size_t dim, calls;
  int flagr = 0;
  double result, abserr;

  if (argc < 4) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 4, 5 or 6)", argc);
  CHECK_MONTE_FUNCTION(argv[0]);
  CHECK_VECTOR(argv[1]);  CHECK_VECTOR(argv[2]);
  Data_Get_Struct(obj, gsl_monte_plain_state, plain);
  Data_Get_Struct(argv[0], gsl_monte_function, F);
  Data_Get_Struct(argv[1], gsl_vector, xl);
  Data_Get_Struct(argv[2], gsl_vector, xu);
  if (argc >= 5 && TYPE(argv[4]) == T_FIXNUM) {
    dim = FIX2INT(argv[3]);
    calls = FIX2INT(argv[4]);
  } else {
    dim = F->dim;
    calls = FIX2INT(argv[3]);
  }

  if (rb_obj_is_kind_of(argv[argc-1], cgsl_rng)) {
    Data_Get_Struct(argv[argc-1], gsl_rng, r);
  } else {
    r = gsl_rng_alloc(gsl_rng_default);
    flagr = 1;
  }
  gsl_monte_plain_integrate(F, xl->data, xu->data, dim, calls, r, plain, 
			    &result, &abserr);
  if (flagr == 1) gsl_rng_free(r);
  return rb_ary_new3(2, rb_float_new(result), rb_float_new(abserr));
}

static VALUE rb_gsl_monte_miser_integrate(int argc, VALUE *argv, VALUE obj)
{
  gsl_monte_function *F = NULL;
  gsl_monte_miser_state *miser = NULL;
  gsl_vector *xl = NULL, *xu = NULL;
  gsl_rng *r = NULL;
  size_t dim, calls;
  int flagr = 0;
  double result, abserr;

  if (argc < 4)
    rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 4)", argc);
  CHECK_MONTE_FUNCTION(argv[0]);
  CHECK_VECTOR(argv[1]);  CHECK_VECTOR(argv[2]);
  Data_Get_Struct(obj, gsl_monte_miser_state, miser);
  Data_Get_Struct(argv[0], gsl_monte_function, F);
  Data_Get_Struct(argv[1], gsl_vector, xl);
  Data_Get_Struct(argv[2], gsl_vector, xu);
  if (argc >= 5 && TYPE(argv[4]) == T_FIXNUM) {
    dim = FIX2INT(argv[3]);
    calls = FIX2INT(argv[4]);
  } else {
    dim = F->dim;
    calls = FIX2INT(argv[3]);
  }

  if (rb_obj_is_kind_of(argv[argc-1], cgsl_rng)) {
    Data_Get_Struct(argv[argc-1], gsl_rng, r);
  } else {
    r = gsl_rng_alloc(gsl_rng_default);
    flagr = 1;
  }
  gsl_monte_miser_integrate(F, xl->data, xu->data, dim, calls, r, miser, 
			    &result, &abserr);
  if (flagr == 1) gsl_rng_free(r);
  return rb_ary_new3(2, rb_float_new(result), rb_float_new(abserr));
}

static VALUE rb_gsl_monte_vegas_integrate(int argc, VALUE *argv, VALUE obj)
{
  gsl_monte_function *F = NULL;
  gsl_monte_vegas_state *vegas = NULL;
  gsl_vector *xl = NULL, *xu = NULL;
  gsl_rng *r = NULL;
  size_t dim, calls;
  int flagr = 0;
  double result, abserr;

  if (argc < 4) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 4)", argc);
  CHECK_MONTE_FUNCTION(argv[0]);
  CHECK_VECTOR(argv[1]);  CHECK_VECTOR(argv[2]);
  Data_Get_Struct(obj, gsl_monte_vegas_state, vegas);
  Data_Get_Struct(argv[0], gsl_monte_function, F);
  Data_Get_Struct(argv[1], gsl_vector, xl);
  Data_Get_Struct(argv[2], gsl_vector, xu);
  if (argc >= 5 && TYPE(argv[4]) == T_FIXNUM) {
    dim = FIX2INT(argv[3]);
    calls = FIX2INT(argv[4]);
  } else {
    dim = F->dim;
    calls = FIX2INT(argv[3]);
  }

  if (rb_obj_is_kind_of(argv[argc-1], cgsl_rng)) {
    Data_Get_Struct(argv[argc-1], gsl_rng, r);
  } else {
    r = gsl_rng_alloc(gsl_rng_default);
    flagr = 1;
  }
  gsl_monte_vegas_integrate(F, xl->data, xu->data, dim, calls, r, vegas, 
			    &result, &abserr);
  if (flagr == 1) gsl_rng_free(r);
  return rb_ary_new3(2, rb_float_new(result), rb_float_new(abserr));
}

static int get_monte_type(VALUE vt)
{
  char name[32];
  if (rb_obj_is_kind_of(vt, cgsl_monte_plain)) return GSL_MONTE_PLAIN_STATE;
  else if (rb_obj_is_kind_of(vt, cgsl_monte_miser)) return GSL_MONTE_MISER_STATE;
  else if (rb_obj_is_kind_of(vt, cgsl_monte_vegas)) return GSL_MONTE_VEGAS_STATE;
  else {
    /* do next */
  }

  switch(TYPE(vt)) {
  case T_STRING:
    strcpy(name, STR2CSTR(vt));
    if (str_tail_grep(name, "plain") == 0) {
      return GSL_MONTE_PLAIN_STATE + 100;
    } else if (str_tail_grep(name, "miser") == 0) {
      return GSL_MONTE_MISER_STATE + 100;
    } else if (str_tail_grep(name, "vegas") == 0) {
      return GSL_MONTE_VEGAS_STATE + 100;
    } else {
      rb_raise(rb_eArgError, "%s: unknown algorithm", name);
    }
    break;
  case T_FIXNUM:
    return FIX2INT(vt) + 100;
    break;
  default:
    rb_raise(rb_eTypeError, "String or Fixnum expected");
    break;
  }
  /* wrong argument if reach here */
  rb_raise(rb_eArgError, "wrong argument");
}

static VALUE rb_gsl_monte_miser_estimate_frac(VALUE obj)
{
  gsl_monte_miser_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  return rb_float_new(s->estimate_frac);
}

static VALUE rb_gsl_monte_miser_set_estimate_frac(VALUE obj, VALUE val)
{
  gsl_monte_miser_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  s->estimate_frac = NUM2DBL(val);
  return obj;
}

static VALUE rb_gsl_monte_miser_min_calls(VALUE obj)
{
  gsl_monte_miser_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  return INT2FIX(s->min_calls);
}

static VALUE rb_gsl_monte_miser_set_min_calls(VALUE obj, VALUE val)
{
  gsl_monte_miser_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  s->min_calls = FIX2INT(val);
  return obj;
}

static VALUE rb_gsl_monte_miser_min_calls_per_bisection(VALUE obj)
{
  gsl_monte_miser_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  return INT2FIX(s->min_calls_per_bisection);
}

static VALUE rb_gsl_monte_miser_set_min_calls_per_bisection(VALUE obj, VALUE val)
{
  gsl_monte_miser_state *s = NULL;
  CHECK_FIXNUM(val);
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  s->min_calls_per_bisection = FIX2INT(val);
  return obj;
}

static VALUE rb_gsl_monte_miser_alpha(VALUE obj)
{
  gsl_monte_miser_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  return rb_float_new(s->alpha);
}

static VALUE rb_gsl_monte_miser_set_alpha(VALUE obj, VALUE val)
{
  gsl_monte_miser_state *s = NULL;
  Need_Float(val);
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  s->alpha = NUM2DBL(val);
  return obj;
}

static VALUE rb_gsl_monte_miser_dither(VALUE obj)
{
  gsl_monte_miser_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  return rb_float_new(s->dither);
}

static VALUE rb_gsl_monte_miser_set_dither(VALUE obj, VALUE val)
{
  gsl_monte_miser_state *s = NULL;
  Need_Float(val);
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  s->dither = NUM2DBL(val);
  return obj;
}

static VALUE rb_gsl_monte_miser_state(VALUE obj)
{
  gsl_monte_miser_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  return rb_ary_new3(5, rb_float_new(s->estimate_frac), INT2FIX(s->min_calls),
		     INT2FIX(s->min_calls_per_bisection), rb_float_new(s->alpha),
		     rb_float_new(s->dither));
}

static VALUE rb_gsl_monte_vegas_result(VALUE obj)
{
  gsl_monte_vegas_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  return rb_float_new(s->result);
}

static VALUE rb_gsl_monte_vegas_sigma(VALUE obj)
{
  gsl_monte_vegas_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  return rb_float_new(s->sigma);
}
static VALUE rb_gsl_monte_vegas_chisq(VALUE obj)
{
  gsl_monte_vegas_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  return rb_float_new(s->chisq);
}
static VALUE rb_gsl_monte_vegas_alpha(VALUE obj)
{
  gsl_monte_vegas_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  return rb_float_new(s->alpha);
}

static VALUE rb_gsl_monte_vegas_set_alpha(VALUE obj, VALUE val)
{
  gsl_monte_vegas_state *s = NULL;
  Need_Float(val);
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  s->alpha = NUM2DBL(val);
  return obj;
}

static VALUE rb_gsl_monte_vegas_iterations(VALUE obj)
{
  gsl_monte_vegas_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  return INT2FIX(s->iterations);
}

static VALUE rb_gsl_monte_vegas_set_iterations(VALUE obj, VALUE val)
{
  gsl_monte_vegas_state *s = NULL;
  CHECK_FIXNUM(val);
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  s->iterations = FIX2INT(val);
  return obj;
}
static VALUE rb_gsl_monte_vegas_stage(VALUE obj)
{
  gsl_monte_vegas_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  return INT2FIX(s->stage);
}
static VALUE rb_gsl_monte_vegas_set_stage(VALUE obj, VALUE val)
{
  gsl_monte_vegas_state *s = NULL;
  CHECK_FIXNUM(val);
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  s->stage = FIX2INT(val);
  return obj;
}
static VALUE rb_gsl_monte_vegas_mode(VALUE obj)
{
  gsl_monte_vegas_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  return INT2FIX(s->mode);
}
static VALUE rb_gsl_monte_vegas_set_mode(VALUE obj, VALUE val)
{
  gsl_monte_vegas_state *s = NULL;
  CHECK_FIXNUM(val);
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  s->mode = FIX2INT(val);
  return obj;
}
static VALUE rb_gsl_monte_vegas_verbose(VALUE obj)
{
  gsl_monte_vegas_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  return INT2FIX(s->verbose);
}
static VALUE rb_gsl_monte_vegas_set_verbose(VALUE obj, VALUE val)
{
  gsl_monte_vegas_state *s = NULL;
  CHECK_FIXNUM(val);
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  s->verbose = FIX2INT(val);
  return obj;
}
static VALUE rb_gsl_monte_vegas_state(VALUE obj)
{
  gsl_monte_vegas_state *s = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  return rb_ary_new3(8, rb_float_new(s->result), rb_float_new(s->sigma), 
		     rb_float_new(s->chisq), rb_float_new(s->alpha), 
		     INT2FIX(s->iterations), INT2FIX(s->stage), 
		     INT2FIX(s->mode), INT2FIX(s->verbose));
}

#ifdef GSL_1_13_LATER
static VALUE rb_gsl_monte_miser_params_get(VALUE obj)
{
  gsl_monte_miser_state *s = NULL;
  gsl_monte_miser_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  p = (gsl_monte_miser_params *) malloc(sizeof(gsl_monte_miser_params));
  gsl_monte_miser_params_get(s, p);  
  return Data_Wrap_Struct(cgsl_monte_miser_params, 0, free, p);
}
static VALUE rb_gsl_monte_miser_params_set(VALUE obj, VALUE params)
{
  gsl_monte_miser_state *s = NULL;
  gsl_monte_miser_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_state, s);
  Data_Get_Struct(params, gsl_monte_miser_params, p);
  gsl_monte_miser_params_set(s, p);  
  return Qtrue;
}
static VALUE rb_gsl_monte_miser_params_get_estimate_frac(VALUE obj)
{
  gsl_monte_miser_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_params, p);
  return rb_float_new(p->estimate_frac);
}
static VALUE rb_gsl_monte_miser_params_set_estimate_frac(VALUE obj, VALUE val)
{
  gsl_monte_miser_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_params, p);
  p->estimate_frac = NUM2DBL(val);
  return val;
}
static VALUE rb_gsl_monte_miser_params_get_min_calls(VALUE obj)
{
  gsl_monte_miser_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_params, p);
  return INT2FIX(p->min_calls);
}
static VALUE rb_gsl_monte_miser_params_set_min_calls(VALUE obj, VALUE val)
{
  gsl_monte_miser_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_params, p);
  p->min_calls = (size_t) FIX2INT(val);
  return val;
}
static VALUE rb_gsl_monte_miser_params_get_min_calls_per_bisection(VALUE obj)
{
  gsl_monte_miser_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_params, p);
  return INT2FIX(p->min_calls_per_bisection);
}
static VALUE rb_gsl_monte_miser_params_set_min_calls_per_bisection(VALUE obj, VALUE val)
{
  gsl_monte_miser_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_params, p);
  p->min_calls_per_bisection = (size_t) FIX2INT(val);
  return val;
}
static VALUE rb_gsl_monte_miser_params_get_alpha(VALUE obj)
{
  gsl_monte_miser_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_params, p);
  return rb_float_new(p->alpha);
}
static VALUE rb_gsl_monte_miser_params_set_alpha(VALUE obj, VALUE val)
{
  gsl_monte_miser_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_params, p);
  p->alpha = NUM2DBL(val);
  return val;
}
static VALUE rb_gsl_monte_miser_params_get_dither(VALUE obj)
{
  gsl_monte_miser_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_params, p);
  return rb_float_new(p->dither);
}
static VALUE rb_gsl_monte_miser_params_set_dither(VALUE obj, VALUE val)
{
  gsl_monte_miser_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_miser_params, p);
  p->dither = NUM2DBL(val);
  return val;
}

static VALUE rb_gsl_monte_vegas_params_get(VALUE obj)
{
  gsl_monte_vegas_state *s = NULL;
  gsl_monte_vegas_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  p = (gsl_monte_vegas_params *) malloc(sizeof(gsl_monte_vegas_params));
  gsl_monte_vegas_params_get(s, p);  
  return Data_Wrap_Struct(cgsl_monte_vegas_params, 0, free, p);
}
static VALUE rb_gsl_monte_vegas_params_set(VALUE obj, VALUE params)
{
  gsl_monte_vegas_state *s = NULL;
  gsl_monte_vegas_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  Data_Get_Struct(params, gsl_monte_vegas_params, p);
  gsl_monte_vegas_params_set(s, p);  
  return Qtrue;
}
static VALUE rb_gsl_monte_vegas_params_get_alpha(VALUE obj)
{
  gsl_monte_vegas_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_params, p);
  return rb_float_new(p->alpha);
}
static VALUE rb_gsl_monte_vegas_params_set_alpha(VALUE obj, VALUE val)
{
  gsl_monte_vegas_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_params, p);
  p->alpha = NUM2DBL(val);
  return val;
}
static VALUE rb_gsl_monte_vegas_params_get_iterations(VALUE obj)
{
  gsl_monte_vegas_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_params, p);
  return INT2FIX(p->iterations);
}
static VALUE rb_gsl_monte_vegas_params_set_iterations(VALUE obj, VALUE val)
{
  gsl_monte_vegas_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_params, p);
  p->iterations = (size_t) FIX2INT(val);
  return val;
}
static VALUE rb_gsl_monte_vegas_params_get_stage(VALUE obj)
{
  gsl_monte_vegas_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_params, p);
  return INT2FIX(p->stage);
}
static VALUE rb_gsl_monte_vegas_params_set_stage(VALUE obj, VALUE val)
{
  gsl_monte_vegas_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_params, p);
  p->stage = FIX2INT(val);
  return val;
}
static VALUE rb_gsl_monte_vegas_params_get_mode(VALUE obj)
{
  gsl_monte_vegas_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_params, p);
  return INT2FIX(p->mode);
}
static VALUE rb_gsl_monte_vegas_params_set_mode(VALUE obj, VALUE val)
{
  gsl_monte_vegas_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_params, p);
  p->mode = FIX2INT(val);
  return val;
}
static VALUE rb_gsl_monte_vegas_params_get_verbose(VALUE obj)
{
  gsl_monte_vegas_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_params, p);
  return INT2FIX(p->verbose);
}
static VALUE rb_gsl_monte_vegas_params_set_verbose(VALUE obj, VALUE val)
{
  gsl_monte_vegas_params *p = NULL;
  Data_Get_Struct(obj, gsl_monte_vegas_params, p);
  p->verbose = FIX2INT(val);
  return val;
}
static VALUE rb_gsl_monte_vegas_runval(VALUE obj)
{
  gsl_monte_vegas_state *s = NULL;
  double res, sig;
  VALUE ary;
  Data_Get_Struct(obj, gsl_monte_vegas_state, s);
  gsl_monte_vegas_runval(s, &res, &sig);
  ary = rb_ary_new2(2);
  rb_ary_store(ary, 0, rb_float_new(res));
  rb_ary_store(ary, 1, rb_float_new(sig));
  return ary;
}
#endif

void Init_gsl_monte(VALUE module)
{
  VALUE mgsl_monte;


  mgsl_monte = rb_define_module_under(module, "Monte");

  rb_define_const(mgsl_monte, "PLAIN", INT2FIX(GSL_MONTE_PLAIN_STATE));
  rb_define_const(mgsl_monte, "MISER", INT2FIX(GSL_MONTE_MISER_STATE));
  rb_define_const(mgsl_monte, "VEGAS", INT2FIX(GSL_MONTE_VEGAS_STATE));

  cgsl_monte_function = rb_define_class_under(mgsl_monte, "Function", cGSL_Object);
  cgsl_monte_plain = rb_define_class_under(mgsl_monte, "Plain", cGSL_Object);
  cgsl_monte_miser = rb_define_class_under(mgsl_monte, "Miser", cGSL_Object);
  cgsl_monte_vegas = rb_define_class_under(mgsl_monte, "Vegas", cGSL_Object);

  rb_define_singleton_method(cgsl_monte_function, "new", rb_gsl_monte_function_new, -1);
  rb_define_singleton_method(cgsl_monte_function, "alloc", rb_gsl_monte_function_new, -1);

  rb_define_method(cgsl_monte_function, "proc", rb_gsl_monte_function_proc, 0);
  rb_define_method(cgsl_monte_function, "eval", rb_gsl_monte_function_eval, 0);
  rb_define_alias(cgsl_monte_function, "call", "eval");
  rb_define_method(cgsl_monte_function, "params", rb_gsl_monte_function_params, 0);
  rb_define_method(cgsl_monte_function, "set", rb_gsl_monte_function_set_f, -1);
  rb_define_alias(cgsl_monte_function, "set_proc", "set");
  rb_define_method(cgsl_monte_function, "set_params", rb_gsl_monte_function_set_params, -1);

  rb_define_method(cgsl_monte_function, "integrate", rb_gsl_monte_integrate, -1);

  /*****/
  rb_define_singleton_method(cgsl_monte_plain, "new", rb_gsl_monte_plain_new, 1);
  rb_define_singleton_method(cgsl_monte_plain, "alloc", rb_gsl_monte_plain_new, 1);
  rb_define_method(cgsl_monte_plain, "init", rb_gsl_monte_plain_init, 0);

  rb_define_singleton_method(cgsl_monte_miser, "new", rb_gsl_monte_miser_new, 1);
  rb_define_singleton_method(cgsl_monte_miser, "alloc", rb_gsl_monte_miser_new, 1);
  rb_define_method(cgsl_monte_miser, "init", rb_gsl_monte_miser_init, 0);
  rb_define_method(cgsl_monte_miser, "estimate_frac", rb_gsl_monte_miser_estimate_frac, 0);
  rb_define_method(cgsl_monte_miser, "min_calls", rb_gsl_monte_miser_min_calls, 0);
  rb_define_method(cgsl_monte_miser, "min_calls_per_bisection", rb_gsl_monte_miser_min_calls_per_bisection, 0);
  rb_define_method(cgsl_monte_miser, "alpha", rb_gsl_monte_miser_alpha, 0);
  rb_define_method(cgsl_monte_miser, "dither", rb_gsl_monte_miser_dither, 0);
  rb_define_method(cgsl_monte_miser, "state", rb_gsl_monte_miser_state, 0);

  rb_define_method(cgsl_monte_miser, "set_estimate_frac", rb_gsl_monte_miser_set_estimate_frac, 1);
  rb_define_alias(cgsl_monte_miser, "estimate_frac=", "set_estimate_frac");
  rb_define_method(cgsl_monte_miser, "set_min_calls", rb_gsl_monte_miser_set_min_calls, 1);
  rb_define_alias(cgsl_monte_miser, "min_calls=", "set_min_calls");
  rb_define_method(cgsl_monte_miser, "set_min_calls_per_bisection", rb_gsl_monte_miser_set_min_calls_per_bisection, 1);
  rb_define_alias(cgsl_monte_miser, "min_calls_per_bisection=", "set_min_calls_per_bisection");
  rb_define_method(cgsl_monte_miser, "set_alpha", rb_gsl_monte_miser_set_alpha, 1);
  rb_define_alias(cgsl_monte_miser, "alpha=", "set_alpha");
  rb_define_method(cgsl_monte_miser, "set_dither", rb_gsl_monte_miser_set_dither, 1);
  rb_define_alias(cgsl_monte_miser, "dither=", "set_dither");

  /*****/
  rb_define_singleton_method(cgsl_monte_vegas, "new", rb_gsl_monte_vegas_new, 1);
  rb_define_singleton_method(cgsl_monte_vegas, "alloc", rb_gsl_monte_vegas_new, 1);
  rb_define_method(cgsl_monte_vegas, "init", rb_gsl_monte_vegas_init, 0);

  rb_define_method(cgsl_monte_vegas, "result", rb_gsl_monte_vegas_result, 0);
  rb_define_method(cgsl_monte_vegas, "sigma", rb_gsl_monte_vegas_sigma, 0);
  rb_define_method(cgsl_monte_vegas, "chisq", rb_gsl_monte_vegas_chisq, 0);
  rb_define_method(cgsl_monte_vegas, "alpha", rb_gsl_monte_vegas_alpha, 0);
  rb_define_method(cgsl_monte_vegas, "iterations", rb_gsl_monte_vegas_iterations, 0);
  rb_define_method(cgsl_monte_vegas, "stage", rb_gsl_monte_vegas_stage, 0);
  rb_define_method(cgsl_monte_vegas, "mode", rb_gsl_monte_vegas_mode, 0);
  rb_define_method(cgsl_monte_vegas, "verbose", rb_gsl_monte_vegas_verbose, 0);
  rb_define_method(cgsl_monte_vegas, "state", rb_gsl_monte_vegas_state, 0);

  rb_define_method(cgsl_monte_vegas, "set_alpha", rb_gsl_monte_vegas_set_alpha, 1);
  rb_define_alias(cgsl_monte_vegas, "alpha=", "set_alpha");
  rb_define_method(cgsl_monte_vegas, "set_iterations", rb_gsl_monte_vegas_set_iterations, 1);
  rb_define_alias(cgsl_monte_vegas, "iterations=", "set_iterations");
  rb_define_method(cgsl_monte_vegas, "set_stage", rb_gsl_monte_vegas_set_stage, 1);
  rb_define_alias(cgsl_monte_vegas, "stage=", "set_stage");
  rb_define_method(cgsl_monte_vegas, "set_mode", rb_gsl_monte_vegas_set_mode, 1);
  rb_define_alias(cgsl_monte_vegas, "mode=", "set_mode");
  rb_define_method(cgsl_monte_vegas, "set_verbose", rb_gsl_monte_vegas_set_verbose, 1);
  rb_define_alias(cgsl_monte_vegas, "verbose=", "set_verbose");


  /*****/
  rb_define_singleton_method(cgsl_monte_plain, "integrate", 
			     rb_gsl_monte_integrate, -1);
  rb_define_method(cgsl_monte_plain, "integrate", 
		   rb_gsl_monte_plain_integrate, -1);
  rb_define_singleton_method(cgsl_monte_miser, "integrate", 
			     rb_gsl_monte_integrate, -1);
  rb_define_method(cgsl_monte_miser, "integrate", 
		   rb_gsl_monte_miser_integrate, -1);
  rb_define_singleton_method(cgsl_monte_vegas, "integrate", 
			     rb_gsl_monte_integrate, -1);
  rb_define_method(cgsl_monte_vegas, "integrate", 
		   rb_gsl_monte_vegas_integrate, -1);

#ifdef GSL_1_13_LATER
  cgsl_monte_miser_params = rb_define_class_under(cgsl_monte_miser, "Params", cGSL_Object);
  cgsl_monte_vegas_params = rb_define_class_under(cgsl_monte_vegas, "Params", cGSL_Object);

  rb_define_method(cgsl_monte_miser, "params_get", rb_gsl_monte_miser_params_get, 0);
  rb_define_method(cgsl_monte_miser, "params_set", rb_gsl_monte_miser_params_set, 1);
  rb_define_method(cgsl_monte_miser_params, "estimate_frac", rb_gsl_monte_miser_params_get_estimate_frac, 0);
  rb_define_method(cgsl_monte_miser_params, "set_estimate_frac", rb_gsl_monte_miser_params_set_estimate_frac, 1);
  rb_define_alias(cgsl_monte_miser_params, "estimate_frac=", "set_estimate_frac");
  rb_define_method(cgsl_monte_miser_params, "min_calls", rb_gsl_monte_miser_params_get_min_calls, 0);
  rb_define_method(cgsl_monte_miser_params, "set_min_calls", rb_gsl_monte_miser_params_set_min_calls, 1);
  rb_define_alias(cgsl_monte_miser_params, "min_calls=", "set_min_calls");
  rb_define_method(cgsl_monte_miser_params, "min_calls_per_bisection", rb_gsl_monte_miser_params_get_min_calls_per_bisection, 0);
  rb_define_method(cgsl_monte_miser_params, "set_min_calls_per_bisection", rb_gsl_monte_miser_params_set_min_calls_per_bisection, 1);
  rb_define_alias(cgsl_monte_miser_params, "min_calls_per_bisection=", "set_min_calls_per_bisection");
  rb_define_method(cgsl_monte_miser_params, "alpha", rb_gsl_monte_miser_params_get_alpha, 0);
  rb_define_method(cgsl_monte_miser_params, "set_alpha", rb_gsl_monte_miser_params_set_alpha, 1);
  rb_define_alias(cgsl_monte_miser_params, "alpha=", "set_alpha");
  rb_define_method(cgsl_monte_miser_params, "dither", rb_gsl_monte_miser_params_get_dither, 0);
  rb_define_method(cgsl_monte_miser_params, "set_dither", rb_gsl_monte_miser_params_set_dither, 1);
  rb_define_alias(cgsl_monte_miser_params, "dither=", "set_dither");

  rb_define_method(cgsl_monte_vegas, "params_get", rb_gsl_monte_vegas_params_get, 0);
  rb_define_method(cgsl_monte_vegas, "params_set", rb_gsl_monte_vegas_params_set, 1);
  rb_define_method(cgsl_monte_vegas_params, "alpha", rb_gsl_monte_vegas_params_get_alpha, 0);
  rb_define_method(cgsl_monte_vegas_params, "set_alpha", rb_gsl_monte_vegas_params_set_alpha, 1);
  rb_define_alias(cgsl_monte_vegas_params, "alpha=", "set_alpha");
  rb_define_method(cgsl_monte_vegas_params, "iterations", rb_gsl_monte_vegas_params_get_iterations, 0);
  rb_define_method(cgsl_monte_vegas_params, "set_iterations", rb_gsl_monte_vegas_params_set_iterations, 1);
  rb_define_alias(cgsl_monte_vegas_params, "iterations=", "set_iterations");
  rb_define_method(cgsl_monte_vegas_params, "stage", rb_gsl_monte_vegas_params_get_stage, 0);
  rb_define_method(cgsl_monte_vegas_params, "set_stage", rb_gsl_monte_vegas_params_set_stage, 1);
  rb_define_alias(cgsl_monte_vegas_params, "stage=", "set_stage");
  rb_define_method(cgsl_monte_vegas_params, "mode", rb_gsl_monte_vegas_params_get_mode, 0);
  rb_define_method(cgsl_monte_vegas_params, "set_mode", rb_gsl_monte_vegas_params_set_mode, 1);
  rb_define_alias(cgsl_monte_vegas_params, "mode=", "set_mode");
  rb_define_method(cgsl_monte_vegas_params, "verbose", rb_gsl_monte_vegas_params_get_verbose, 0);
  rb_define_method(cgsl_monte_vegas_params, "set_verbose", rb_gsl_monte_vegas_params_set_verbose, 1);
  rb_define_alias(cgsl_monte_vegas_params, "verbose=", "set_verbose");

  rb_define_method(cgsl_monte_vegas, "runval", rb_gsl_monte_vegas_runval, 0);

  rb_define_const(cgsl_monte_vegas, "MODE_IMPORTANCE", INT2FIX(GSL_VEGAS_MODE_IMPORTANCE));
  rb_define_const(cgsl_monte_vegas, "MODE_IMPORTANCE_ONLY", INT2FIX(GSL_VEGAS_MODE_IMPORTANCE_ONLY)); 
  rb_define_const(cgsl_monte_vegas, "MODE_STRATIFIED", INT2FIX(GSL_VEGAS_MODE_STRATIFIED));
#endif

}
#ifdef CHECK_MONTE_FUNCTION
#undef CHECK_MONTE_FUNCTION
#endif
