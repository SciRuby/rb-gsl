/*
  gsl.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl.h"
#include <gsl/gsl_machine.h>

ID rb_gsl_id_beg, rb_gsl_id_end, rb_gsl_id_excl, rb_gsl_id_to_a;
static ID rb_gsl_id_name, rb_gsl_id_size;
VALUE cGSL_Object;
static void rb_gsl_define_intern(VALUE module);
static void rb_gsl_define_const(VALUE module);
static void rb_gsl_define_methods(VALUE module);

static VALUE rb_gsl_object_inspect(VALUE obj)
{
  char buf[256];
  sprintf(buf, "%s", rb_class2name(CLASS_OF(obj)));
  return rb_str_new2(buf);
}

static VALUE rb_gsl_call_rescue(VALUE obj)
{
  return Qfalse;
}

static VALUE rb_gsl_call_name(VALUE obj)
{
  return rb_funcall(obj, rb_gsl_id_name, 0);
}

static VALUE rb_gsl_call_size(VALUE obj)
{
  return rb_funcall(obj, rb_gsl_id_size, 0);
}

static VALUE rb_gsl_object_info(VALUE obj)
{
  char buf[256];
  VALUE s;
  sprintf(buf, "Class:      %s\n", rb_class2name(CLASS_OF(obj)));
  sprintf(buf, "%sSuperClass: %s\n", buf, rb_class2name(RCLASS_SUPER(CLASS_OF(obj))));
  s = rb_rescue(rb_gsl_call_name, obj, rb_gsl_call_rescue, obj);
  if (s) sprintf(buf, "%sType:       %s\n", buf, STR2CSTR(s));
  s = rb_rescue(rb_gsl_call_size, obj, rb_gsl_call_rescue, obj);
  if (s) sprintf(buf, "%sSize:       %d\n", buf, (int) FIX2INT(s));
  return rb_str_new2(buf);
}

static VALUE rb_gsl_not_implemeted(VALUE obj)
{
  rb_raise(rb_eNotImpError, "%s#dup is not implemented", rb_class2name(CLASS_OF(obj)));
  return Qnil;
}

void Init_rb_gsl()
{
  VALUE mgsl;

  mgsl = rb_define_module("GSL");
  cGSL_Object = rb_define_class_under(mgsl, "Object", rb_cObject);
  rb_define_method(cGSL_Object, "inspect", rb_gsl_object_inspect, 0);
  rb_define_method(cGSL_Object, "info", rb_gsl_object_info, 0);
  // Override Object#dup to prevent invalid dup-ing of GSL_Objects.
  // Subclasses (e.g. GSL::Vector) must provide their own implementation of
  // #dup if desired.
  rb_define_method(cGSL_Object, "dup", rb_gsl_not_implemeted, 0);

  rb_gsl_define_intern(mgsl);

  Init_gsl_error(mgsl);

  Init_gsl_math(mgsl);
  Init_gsl_complex(mgsl);

  Init_gsl_array(mgsl);

  Init_gsl_blas(mgsl);

  Init_gsl_sort(mgsl);

  Init_gsl_poly(mgsl);  /* this must be called after the Vector class defined */
  Init_gsl_poly_int(mgsl);
  Init_gsl_poly2(mgsl);
  Init_gsl_rational(mgsl);

  Init_gsl_sf(mgsl);

  Init_gsl_linalg(mgsl); /*  Init_gsl_linalg_complex() is called in Init_gsl_linalg() */

  Init_gsl_eigen(mgsl);

  Init_gsl_fft(mgsl);
    Init_gsl_signal(mgsl);
  Init_gsl_function(mgsl);
  Init_gsl_integration(mgsl);

  Init_gsl_rng(mgsl);
  Init_gsl_qrng(mgsl);
  Init_gsl_ran(mgsl);
#ifdef GSL_1_4_LATER
  Init_gsl_cdf(mgsl);
#endif
  Init_gsl_stats(mgsl);

  Init_gsl_histogram(mgsl);
  Init_gsl_histogram2d(mgsl);
  Init_gsl_histogram3d(mgsl);
  Init_gsl_ntuple(mgsl);
  Init_gsl_monte(mgsl);
  Init_gsl_siman(mgsl);

  Init_gsl_odeiv(mgsl);
  Init_gsl_interp(mgsl);
  Init_gsl_spline(mgsl);
  Init_gsl_diff(mgsl);
#ifdef GSL_1_4_9_LATER
  Init_gsl_deriv(mgsl);
#endif

  Init_gsl_cheb(mgsl);
  Init_gsl_sum(mgsl);
  Init_gsl_dht(mgsl);

  Init_gsl_root(mgsl);
  Init_gsl_multiroot(mgsl);
  Init_gsl_min(mgsl);
  Init_gsl_multimin(mgsl);
  Init_gsl_fit(mgsl);
  Init_gsl_multifit(mgsl);

  Init_gsl_const(mgsl);

  Init_gsl_ieee(mgsl);

#ifdef HAVE_NARRAY_H
  Init_gsl_narray(mgsl);
#endif

  Init_wavelet(mgsl);

  rb_gsl_define_const(mgsl);

#ifdef HAVE_TENSOR_TENSOR_H
  Init_tensor_init(mgsl);
  Init_tensor_int_init(mgsl);
#endif

  Init_gsl_graph(mgsl);
  Init_gsl_dirac(mgsl);

#ifdef HAVE_TAMU_ANOVA_TAMU_ANOVA_H
  Init_tamu_anova(mgsl);
#endif

#ifdef HAVE_OOL_OOL_VERSION_H
	Init_ool(mgsl);
#endif

#ifdef HAVE_JACOBI_H
	Init_jacobi(mgsl);
#endif

#ifdef HAVE_GSL_GSL_CQP_H
	Init_cqp(mgsl);
#endif

	Init_fresnel(mgsl);
	
#ifdef GSL_1_9_LATER
	Init_bspline(mgsl);
#endif

#ifdef HAVE_ALF_ALF_H
	Init_alf(mgsl);
#endif

	Init_geometry(mgsl);

#ifdef GSL_1_14_LATER
	Init_multiset(mgsl);
#endif

  rb_gsl_define_methods(mgsl);
}

/**********/

static void rb_gsl_define_intern(VALUE module)
{
  rb_gsl_id_beg  = rb_intern("begin");
  rb_gsl_id_end  = rb_intern("end");
  rb_gsl_id_excl = rb_intern("exclude_end?");
  rb_gsl_id_to_a = rb_intern("to_a");
  rb_gsl_id_name  = rb_intern("name");
  rb_gsl_id_size  = rb_intern("size");
}

static void rb_gsl_define_const(VALUE module)
{
  rb_define_const(module, "MODE_DEFAULT", INT2FIX(GSL_MODE_DEFAULT));
  rb_define_const(module, "PREC_DOUBLE", INT2FIX(GSL_PREC_DOUBLE));
  rb_define_const(module, "PREC_SINGLE", INT2FIX(GSL_PREC_SINGLE));
  rb_define_const(module, "PREC_APPROX", INT2FIX(GSL_PREC_APPROX));
#ifdef GSL_VERSION
  rb_define_const(module, "VERSION", rb_str_new2(GSL_VERSION));
  rb_define_const(module, "GSL_VERSION", rb_str_new2(GSL_VERSION));
#endif

  rb_define_const(module, "DBL_EPSILON", rb_float_new(GSL_DBL_EPSILON));
  rb_define_const(module, "FLT_EPSILON", rb_float_new(GSL_FLT_EPSILON));
  rb_define_const(module, "MACH_EPS", rb_float_new(GSL_MACH_EPS));
  rb_define_const(module, "SQRT_DBL_EPSILON", rb_float_new(GSL_SQRT_DBL_EPSILON));
  rb_define_const(module, "ROOT3_DBL_EPSILON", rb_float_new(GSL_ROOT3_DBL_EPSILON));
  rb_define_const(module, "ROOT4_DBL_EPSILON", rb_float_new(GSL_ROOT4_DBL_EPSILON));
  rb_define_const(module, "ROOT5_DBL_EPSILON", rb_float_new(GSL_ROOT5_DBL_EPSILON));
  rb_define_const(module, "ROOT6_DBL_EPSILON", rb_float_new(GSL_ROOT6_DBL_EPSILON));
  rb_define_const(module, "LOG_DBL_EPSILON", rb_float_new(GSL_LOG_DBL_EPSILON));

  rb_define_const(module, "DBL_MAX", rb_float_new(GSL_DBL_MAX));
  rb_define_const(module, "SQRT_DBL_MAX", rb_float_new(GSL_SQRT_DBL_MAX));
  rb_define_const(module, "ROOT3_DBL_MAX", rb_float_new(GSL_ROOT3_DBL_MAX));
  rb_define_const(module, "ROOT4_DBL_MAX", rb_float_new(GSL_ROOT4_DBL_MAX));
  rb_define_const(module, "ROOT5_DBL_MAX", rb_float_new(GSL_ROOT5_DBL_MAX));
  rb_define_const(module, "ROOT6_DBL_MAX", rb_float_new(GSL_ROOT6_DBL_MAX));
  rb_define_const(module, "LOG_DBL_MAX", rb_float_new(GSL_LOG_DBL_MAX));

  rb_define_const(module, "DBL_MIN", rb_float_new(GSL_DBL_MIN));
  rb_define_const(module, "SQRT_DBL_MIN", rb_float_new(GSL_SQRT_DBL_MIN));
  rb_define_const(module, "ROOT3_DBL_MIN", rb_float_new(GSL_ROOT3_DBL_MIN));
  rb_define_const(module, "ROOT4_DBL_MIN", rb_float_new(GSL_ROOT4_DBL_MIN));
  rb_define_const(module, "ROOT5_DBL_MIN", rb_float_new(GSL_ROOT5_DBL_MIN));
  rb_define_const(module, "ROOT6_DBL_MIN", rb_float_new(GSL_ROOT6_DBL_MIN));
  rb_define_const(module, "LOG_DBL_MIN", rb_float_new(GSL_LOG_DBL_MIN));

#ifdef GSL_1_14_LATER
  rb_define_const(module, "MAJOR_VERSION", INT2FIX(GSL_MAJOR_VERSION));
  rb_define_const(module, "MINOR_VERSION", INT2FIX(GSL_MINOR_VERSION));
#endif

}

static VALUE rb_gsl_have_tensor(VALUE module)
{
#ifdef HAVE_TENSOR_TENSOR_H
  return Qtrue;
#else
  return Qfalse;
#endif
}

static VALUE rb_gsl_have_narray(VALUE module)
{
#ifdef HAVE_NARRAY_H
  return Qtrue;
#else
  return Qfalse;
#endif
}

static void rb_gsl_define_methods(VALUE module)
{
  rb_define_singleton_method(module, "have_tensor?", rb_gsl_have_tensor, 0);
  rb_define_singleton_method(module, "have_narray?", rb_gsl_have_narray, 0);
}
