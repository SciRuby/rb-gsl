#ifdef HAVE_GSL_GSL_MULTIMIN_FSDF_H
#include "include/rb_gsl.h"
#include "gsl/gsl_multimin_fsdf.h"

static VALUE cfsdf;
#ifndef CHECK_MULTIMIN_FUNCTION_FSDF
#define CHECK_MULTIMIN_FUNCTION_FSDF(x) if(CLASS_OF(x)!=cfsdf)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (GSL::MultiMin::Function_fsdf expected)",\
      rb_class2name(CLASS_OF(x)));
#endif
extern VALUE cgsl_multimin_function_fdf;

static const gsl_multimin_fsdfminimizer_type* get_fsdfminimizer_type(VALUE t)
{
  char name[64];
  switch (TYPE(t)) {
  case T_STRING:
    strcpy(name, STR2CSTR(t));
    if (strcmp(name, "bundle") == 0 || strcmp(name, "bundle_method") == 0) 
      return gsl_multimin_fsdfminimizer_bundle_method;
    else
      rb_raise(rb_eTypeError, "%s: unknown minimizer type", name);
    break;
  default:
    rb_raise(rb_eTypeError, "type is given by a String or a Fixnum");
    break;
  }
}

static VALUE rb_gsl_fsdfminimizer_alloc(VALUE klass, VALUE t, VALUE n)
{
  gsl_multimin_fsdfminimizer *gmf = NULL;
  const gsl_multimin_fsdfminimizer_type *T;
  T = get_fsdfminimizer_type(t);
  gmf = gsl_multimin_fsdfminimizer_alloc(T, FIX2INT(n));
  return Data_Wrap_Struct(klass, 0, gsl_multimin_fsdfminimizer_free, gmf);
}

static VALUE rb_gsl_fsdfminimizer_set(VALUE obj, VALUE ff, VALUE xx, VALUE ss)
{
  gsl_multimin_fsdfminimizer *gmf = NULL;
  gsl_multimin_function_fsdf *F = NULL;
  gsl_vector *x;
	size_t bundle_size;
  int status;
  CHECK_MULTIMIN_FUNCTION_FSDF(ff);
  Data_Get_Struct(obj, gsl_multimin_fsdfminimizer, gmf);
  Data_Get_Struct(ff, gsl_multimin_function_fsdf, F);
  Data_Get_Vector(xx, x);
	bundle_size = (size_t) FIX2INT(ss);
  status = gsl_multimin_fsdfminimizer_set(gmf, F, x, bundle_size);
  return INT2FIX(status);
}

static VALUE rb_gsl_fsdfminimizer_name(VALUE obj)
{
  gsl_multimin_fsdfminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fsdfminimizer, gmf);
  return rb_str_new2(gsl_multimin_fsdfminimizer_name(gmf));
}

static VALUE rb_gsl_fsdfminimizer_iterate(VALUE obj)
{
  gsl_multimin_fsdfminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fsdfminimizer, gmf);
  return INT2FIX(gsl_multimin_fsdfminimizer_iterate(gmf));
}

static VALUE rb_gsl_fsdfminimizer_x(VALUE obj)
{
  gsl_multimin_fsdfminimizer *gmf = NULL;
  gsl_vector *x = NULL;
  Data_Get_Struct(obj, gsl_multimin_fsdfminimizer, gmf);
  x = gsl_multimin_fsdfminimizer_x(gmf);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, x);
}

static VALUE rb_gsl_fsdfminimizer_subgradient(VALUE obj)
{
  gsl_multimin_fsdfminimizer *gmf = NULL;
  gsl_vector *gradient = NULL;
  Data_Get_Struct(obj, gsl_multimin_fsdfminimizer, gmf);
  gradient = gsl_multimin_fsdfminimizer_subgradient(gmf);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, gradient);
}

static VALUE rb_gsl_fsdfminimizer_minimum(VALUE obj)
{
  gsl_multimin_fsdfminimizer *gmf = NULL;
  double min;
  Data_Get_Struct(obj, gsl_multimin_fsdfminimizer, gmf);
  min = gsl_multimin_fsdfminimizer_minimum(gmf);
  return rb_float_new(min);
}

static VALUE rb_gsl_fsdfminimizer_f(VALUE obj)
{
  gsl_multimin_fsdfminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fsdfminimizer, gmf);
  return rb_float_new(gmf->f);
}

static VALUE rb_gsl_fsdfminimizer_restart(VALUE obj)
{
  gsl_multimin_fsdfminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fsdfminimizer, gmf);
  return INT2FIX(gsl_multimin_fsdfminimizer_restart(gmf));
}

static VALUE rb_gsl_fsdfminimizer_test_gradient(VALUE obj, VALUE ea)
{
  gsl_multimin_fsdfminimizer *gmf = NULL;
  gsl_vector *g = NULL;
  Need_Float(ea);
  Data_Get_Struct(obj, gsl_multimin_fsdfminimizer, gmf);
  g = gsl_multimin_fsdfminimizer_subgradient(gmf);
  return INT2FIX(gsl_multimin_test_gradient(g, NUM2DBL(ea)));
}

static VALUE rb_gsl_fsdfminimizer_test_convergence(VALUE obj, VALUE eps)
{
  gsl_multimin_fsdfminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fsdfminimizer, gmf);
  return INT2FIX(gsl_multimin_test_convergence(gmf, NUM2DBL(eps)));
}

static VALUE rb_gsl_fsdfminimizer_eps(VALUE obj)
{
  gsl_multimin_fsdfminimizer *gmf = NULL;
  Data_Get_Struct(obj, gsl_multimin_fsdfminimizer, gmf);
  return rb_float_new(gmf->eps);
}

void Init_multimin_fsdf(VALUE module)
{
	VALUE cmin;
	
	cmin = rb_define_class_under(module, "FsdfMinimizer",  cGSL_Object);
	cfsdf = rb_define_class_under(module, "Function_fsdf", cgsl_multimin_function_fdf);
	
	rb_define_singleton_method(cmin, "alloc", rb_gsl_fsdfminimizer_alloc, 2);
  rb_define_method(cmin, "set", rb_gsl_fsdfminimizer_set, 3);
  rb_define_method(cmin, "name", rb_gsl_fsdfminimizer_name, 0);
  rb_define_method(cmin, "iterate", rb_gsl_fsdfminimizer_iterate, 0);
  rb_define_method(cmin, "x", rb_gsl_fsdfminimizer_x, 0);
  rb_define_method(cmin, "f", rb_gsl_fsdfminimizer_f, 0);
  rb_define_method(cmin, "subgradient", rb_gsl_fsdfminimizer_subgradient, 0);
  rb_define_method(cmin, "minimum", rb_gsl_fsdfminimizer_minimum, 0);
  rb_define_method(cmin, "restart", rb_gsl_fsdfminimizer_restart, 0);
  rb_define_method(cmin, "test_gradient", rb_gsl_fsdfminimizer_test_gradient, 1);
  rb_define_method(cmin, "test_convergence", rb_gsl_fsdfminimizer_test_convergence, 1);  
  rb_define_method(cmin, "eps", rb_gsl_fsdfminimizer_eps, 0);	
}

#endif
