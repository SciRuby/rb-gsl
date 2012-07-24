/*
  min.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl.h"
#include "rb_gsl_function.h"
#include <gsl/gsl_min.h>

double rb_gsl_function_f(double x, void *p); 

enum {
  GSL_MIN_FMINIMIZER_GOLDENSECTION,
  GSL_MIN_FMINIMIZER_BRENT,
#ifdef GSL_1_13_LATER
  GSL_MIN_FMINIMIZER_QUAD_GOLDEN,
#endif
};
static const gsl_min_fminimizer_type* rb_gsl_min_fminimizer_type_get(VALUE t);

static const gsl_min_fminimizer_type* rb_gsl_min_fminimizer_type_get(VALUE t) 
{
  char name[32];
  switch (TYPE(t)) {
  case T_STRING:
    strcpy(name, STR2CSTR(t));
    if (str_tail_grep(name, "goldensection") == 0) 
      return gsl_min_fminimizer_goldensection;
    else if (str_tail_grep(name, "brent") == 0) 
      return gsl_min_fminimizer_brent;
#ifdef GSL_1_13_LATER
    else if (str_tail_grep(name, "quad_golden") == 0) 
      return gsl_min_fminimizer_quad_golden;
#endif
    else 
      rb_raise(rb_eTypeError, "unknown type %s (goldensection, brent or quad_golden expected)",
	       name);
    break;
  case T_FIXNUM:
    switch (FIX2INT(t)) {
    case GSL_MIN_FMINIMIZER_GOLDENSECTION: 
      return gsl_min_fminimizer_goldensection; 
      break;
    case GSL_MIN_FMINIMIZER_BRENT: 
      return gsl_min_fminimizer_brent; 
      break;
#ifdef GSL_1_13_LATER
    case GSL_MIN_FMINIMIZER_QUAD_GOLDEN: 
      return gsl_min_fminimizer_quad_golden; 
      break;
#endif
    default: 
      rb_raise(rb_eTypeError, "unknown type (GOLDENSECION or BRENT or QUAD_GOLDEN expected)"); 
      break;
    }
    break;
  default:
    rb_raise(rb_eTypeError, "wrong argument type %s (String of Fixnum)",
	     rb_class2name(CLASS_OF(t)));
    break;
  }
}

static VALUE rb_gsl_min_fminimizer_new(VALUE klass, VALUE t)
{
  gsl_min_fminimizer *gmf = NULL;
  const gsl_min_fminimizer_type *T;
  T = rb_gsl_min_fminimizer_type_get(t);
  gmf = gsl_min_fminimizer_alloc(T);
  return Data_Wrap_Struct(klass, 0, gsl_min_fminimizer_free, gmf);
}

static VALUE rb_gsl_min_fminimizer_name(VALUE obj)
{
  gsl_min_fminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_min_fminimizer, gmf);
  return rb_str_new2(gsl_min_fminimizer_name(gmf));
}

static VALUE rb_gsl_min_fminimizer_set(VALUE obj, VALUE ff, VALUE xmin, 
				       VALUE xl, VALUE xu)
{
  gsl_min_fminimizer *gmf = NULL;
  gsl_function *f = NULL;
  Need_Float(xmin); Need_Float(xl); Need_Float(xu); 
  CHECK_FUNCTION(ff);
  Data_Get_Struct(obj, gsl_min_fminimizer, gmf);
  Data_Get_Struct(ff, gsl_function, f);
  return INT2FIX(gsl_min_fminimizer_set(gmf, f, NUM2DBL(xmin), 
					NUM2DBL(xl), NUM2DBL(xu)));
}

static VALUE rb_gsl_min_fminimizer_set_with_values(VALUE obj, VALUE ff, 
						   VALUE xmin, VALUE fmin,
						   VALUE xl, VALUE fl, 
						   VALUE xu, VALUE fu)
{
  gsl_min_fminimizer *gmf = NULL;
  gsl_function *f = NULL;
  Need_Float(xmin); Need_Float(xl); Need_Float(xu); 
  Need_Float(fl); Need_Float(fu); 
  CHECK_FUNCTION(ff);
  Data_Get_Struct(obj, gsl_min_fminimizer, gmf);
  Data_Get_Struct(ff, gsl_function, f);
  return INT2FIX(gsl_min_fminimizer_set_with_values(gmf, f, NUM2DBL(xmin), 
						    NUM2DBL(fmin),
						    NUM2DBL(xl), NUM2DBL(fl),
						    NUM2DBL(xu), NUM2DBL(fu)));
}

static VALUE rb_gsl_min_fminimizer_iterate(VALUE obj)
{
  gsl_min_fminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_min_fminimizer, gmf);
  return INT2FIX(gsl_min_fminimizer_iterate(gmf));
}

static VALUE rb_gsl_min_fminimizer_x_lower(VALUE obj)
{
  gsl_min_fminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_min_fminimizer, gmf);
  return rb_float_new(gsl_min_fminimizer_x_lower(gmf));
}

static VALUE rb_gsl_min_fminimizer_x_upper(VALUE obj)
{
  gsl_min_fminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_min_fminimizer, gmf);
  return rb_float_new(gsl_min_fminimizer_x_upper(gmf));
}

#ifndef GSL_1_2_LATER
static double gsl_min_fminimizer_x_minimum(const gsl_min_fminimizer * s)
{  
  /*  return s->x_minimum;*/
  return s->minimum;
}

static double gsl_min_fminimizer_f_minimum(const gsl_min_fminimizer * s)
{  
  return s->f_minimum;
}

static double gsl_min_fminimizer_f_lower(const gsl_min_fminimizer * s)
{ 
  return s->f_lower;
}

static double gsl_min_fminimizer_f_upper(const gsl_min_fminimizer * s)
{ 
  return s->f_upper;
}
#endif

static VALUE rb_gsl_min_fminimizer_x_minimum(VALUE obj)
{
  gsl_min_fminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_min_fminimizer, gmf);
  return rb_float_new(gsl_min_fminimizer_x_minimum(gmf));
}

static VALUE rb_gsl_min_fminimizer_f_minimum(VALUE obj)
{
  gsl_min_fminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_min_fminimizer, gmf);
  return rb_float_new(gsl_min_fminimizer_f_minimum(gmf));
}

static VALUE rb_gsl_min_fminimizer_f_lower(VALUE obj)
{
  gsl_min_fminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_min_fminimizer, gmf);
  return rb_float_new(gsl_min_fminimizer_f_lower(gmf));
}

static VALUE rb_gsl_min_fminimizer_f_upper(VALUE obj)
{
  gsl_min_fminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_min_fminimizer, gmf);
  return rb_float_new(gsl_min_fminimizer_f_upper(gmf));
}

static VALUE rb_gsl_min_fminimizer_test_interval(VALUE obj, VALUE ea, VALUE er)
{
  gsl_min_fminimizer *gmf = NULL;
  double xl, xu;
  Need_Float(ea); Need_Float(er);
  Data_Get_Struct(obj, gsl_min_fminimizer, gmf);
  xl = gsl_min_fminimizer_x_lower(gmf);
  xu = gsl_min_fminimizer_x_upper(gmf);
  return INT2FIX(gsl_min_test_interval(xl, xu, NUM2DBL(ea), NUM2DBL(er)));
}

static VALUE rb_gsl_fminimizer_test_interval(VALUE obj, VALUE xl, VALUE xu,
					     VALUE ea, VALUE er)
{
  Need_Float(xl); Need_Float(xu);
  Need_Float(ea); Need_Float(er);
  return INT2FIX(gsl_min_test_interval(NUM2DBL(xl), NUM2DBL(xu),
				       NUM2DBL(ea), NUM2DBL(er)));
}

void Init_gsl_min(VALUE module)
{
  VALUE mgsl_min, cgsl_fminimizer;

  mgsl_min = rb_define_module_under(module, "Min");
  
  cgsl_fminimizer = rb_define_class_under(mgsl_min, "FMinimizer", cGSL_Object);

  rb_define_const(cgsl_fminimizer, "GOLDENSECTION", 
		  INT2FIX(GSL_MIN_FMINIMIZER_GOLDENSECTION));
  rb_define_const(cgsl_fminimizer, "Goldensection", 
		  INT2FIX(GSL_MIN_FMINIMIZER_GOLDENSECTION));
  rb_define_const(cgsl_fminimizer, "BRENT",
		  INT2FIX(GSL_MIN_FMINIMIZER_BRENT));
  rb_define_const(cgsl_fminimizer, "Brent",
		  INT2FIX(GSL_MIN_FMINIMIZER_BRENT));
#ifdef GSL_1_13_LATER
  rb_define_const(cgsl_fminimizer, "QUAD_GOLDEN",
		  INT2FIX(GSL_MIN_FMINIMIZER_QUAD_GOLDEN));
#endif

  rb_define_singleton_method(cgsl_fminimizer, "new", rb_gsl_min_fminimizer_new, 1);
  rb_define_singleton_method(cgsl_fminimizer, "alloc", rb_gsl_min_fminimizer_new, 1);

  rb_define_method(cgsl_fminimizer, "name", rb_gsl_min_fminimizer_name, 0);
  rb_define_method(cgsl_fminimizer, "set", rb_gsl_min_fminimizer_set, 4);
  rb_define_method(cgsl_fminimizer, "set_with_values", rb_gsl_min_fminimizer_set_with_values, 7);
  rb_define_method(cgsl_fminimizer, "iterate", rb_gsl_min_fminimizer_iterate, 0);

  rb_define_method(cgsl_fminimizer, "x_lower", rb_gsl_min_fminimizer_x_lower, 0);
  rb_define_method(cgsl_fminimizer, "x_upper", rb_gsl_min_fminimizer_x_upper, 0);
  rb_define_method(cgsl_fminimizer, "test_interval", rb_gsl_min_fminimizer_test_interval, 2);

  rb_define_singleton_method(mgsl_min, "test_interval", 
			     rb_gsl_fminimizer_test_interval, 4);

  rb_define_method(cgsl_fminimizer, "x_minimum", rb_gsl_min_fminimizer_x_minimum, 0);
  rb_define_method(cgsl_fminimizer, "f_minimum", rb_gsl_min_fminimizer_f_minimum, 0);
  rb_define_method(cgsl_fminimizer, "f_lower", rb_gsl_min_fminimizer_f_lower, 0);
  rb_define_method(cgsl_fminimizer, "f_upper", rb_gsl_min_fminimizer_f_upper, 0);

}
