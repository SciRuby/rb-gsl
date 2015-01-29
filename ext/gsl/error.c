/*
  error.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/
#include "include/rb_gsl.h"
#include <gsl/gsl_errno.h>
#include "include/rb_gsl_array.h"
#include "include/rb_gsl_common.h"
#include "include/rb_gsl_function.h"

static VALUE eHandler;
static VALUE cgsl_error[35];
static VALUE *pgsl_error;

static void Init_rb_gsl_define_GSL_CONST(VALUE module);
void rb_gsl_error_handler(const char *reason, const char *file,
        int line, int gsl_errno);
static void rb_gsl_my_error_handler(const char *reason, const char *file,
            int line, int gsl_errno);

void rb_gsl_error_handler(const char *reason, const char *file,
        int line, int gsl_errno)
{
  const char *emessage = gsl_strerror(gsl_errno);
  rb_raise(pgsl_error[gsl_errno],
     "Ruby/GSL error code %d, %s (file %s, line %d), %s",
     gsl_errno, reason, file, line, emessage);
}

static void rb_gsl_my_error_handler(const char *reason, const char *file,
            int line, int gsl_errno)
{
  VALUE vreason, vfile;
  VALUE vline, verrno;
  vreason = rb_str_new2(reason);
  vfile = rb_str_new2(file);
  vline = INT2FIX(line);
  verrno = INT2FIX(gsl_errno);
  rb_funcall(eHandler, RBGSL_ID_call, 4, vreason, vfile, vline, verrno);
}

static VALUE rb_gsl_set_error_handler(int argc, VALUE *argv, VALUE module)
{
  if (rb_block_given_p()) {
    eHandler = rb_block_proc();
    gsl_set_error_handler(&rb_gsl_my_error_handler);
    return Qtrue;
  }
  switch (argc) {
  case 0:
    gsl_set_error_handler(&rb_gsl_error_handler);
    return Qtrue;
    break;
  case 1:
    CHECK_PROC(argv[0]);
    eHandler = argv[0];
    gsl_set_error_handler(&rb_gsl_my_error_handler);
    return Qtrue;
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments (%d for 0 or 1 Proc)", argc);
    break;
  }
}

static VALUE rb_gsl_set_default_error_handler(VALUE module)
{
  gsl_set_error_handler(&rb_gsl_error_handler);
  return Qtrue;
}

static void rb_gsl_define_exceptions(VALUE module)
{
  VALUE mgsl_error;
  mgsl_error = rb_define_module_under(module, "ERROR");
  pgsl_error = &cgsl_error[2];
  pgsl_error[-2] = rb_define_class_under(mgsl_error, "CONTINUE", rb_cFixnum);
  pgsl_error[-1] = rb_define_class_under(mgsl_error, "FAILURE", rb_eRuntimeError);
  pgsl_error[0] = rb_define_class_under(mgsl_error, "SUCCESS", rb_cFixnum);
  pgsl_error[1] = rb_define_class_under(mgsl_error, "EDOM", rb_eRangeError);
  pgsl_error[2] = rb_define_class_under(mgsl_error, "ERANGE", rb_eRangeError);
  pgsl_error[3] = rb_define_class_under(mgsl_error, "EFAULT", rb_eRuntimeError);
  pgsl_error[4] = rb_define_class_under(mgsl_error, "EINVAL", rb_eIndexError);
  pgsl_error[5] = rb_define_class_under(mgsl_error, "EFAILED", rb_eRuntimeError);
  pgsl_error[6] = rb_define_class_under(mgsl_error, "EFACTOR", rb_eRuntimeError);
  pgsl_error[7] = rb_define_class_under(mgsl_error, "ESANITY", rb_eRuntimeError);
  pgsl_error[8] = rb_define_class_under(mgsl_error, "ENOMEM", rb_eNoMemError);
  pgsl_error[9] = rb_define_class_under(mgsl_error, "EBADFUNC", rb_eRuntimeError);
  pgsl_error[10] = rb_define_class_under(mgsl_error, "ERUNAWAY", rb_eRuntimeError);
  pgsl_error[11] = rb_define_class_under(mgsl_error, "EMAXITER", rb_eRuntimeError);
  pgsl_error[12] = rb_define_class_under(mgsl_error, "EZERODIV", rb_eZeroDivError);
  pgsl_error[13] = rb_define_class_under(mgsl_error, "EBADTOL", rb_eRuntimeError);
  pgsl_error[14] = rb_define_class_under(mgsl_error, "ETOL", rb_eRuntimeError);
  pgsl_error[15] = rb_define_class_under(mgsl_error, "EUNDRFLW", rb_eRangeError);
  pgsl_error[16] = rb_define_class_under(mgsl_error, "EOVRFLW", rb_eRangeError);
  pgsl_error[17] = rb_define_class_under(mgsl_error, "ELOSS", rb_eRuntimeError);
  pgsl_error[18] = rb_define_class_under(mgsl_error, "EROUND", rb_eRuntimeError);
  pgsl_error[19] = rb_define_class_under(mgsl_error, "EBADLEN", rb_eIndexError);
  pgsl_error[20] = rb_define_class_under(mgsl_error, "ENOTSQR", rb_eRuntimeError);
  pgsl_error[21] = rb_define_class_under(mgsl_error, "ESING", rb_eRuntimeError);
  pgsl_error[22] = rb_define_class_under(mgsl_error, "EDIVERGE", rb_eRuntimeError);
  pgsl_error[23] = rb_define_class_under(mgsl_error, "EUNSUP", rb_eRuntimeError);
  pgsl_error[24] = rb_define_class_under(mgsl_error, "EUNIMPL", rb_eNotImpError);
  pgsl_error[25] = rb_define_class_under(mgsl_error, "ECACHE", rb_eRuntimeError);
  pgsl_error[26] = rb_define_class_under(mgsl_error, "ETABLE", rb_eRuntimeError);
  pgsl_error[27] = rb_define_class_under(mgsl_error, "ENOPROG", rb_eRuntimeError);
  pgsl_error[28] = rb_define_class_under(mgsl_error, "ENOPROGJ", rb_eRuntimeError);
  pgsl_error[29] = rb_define_class_under(mgsl_error, "ETOLF", rb_eRuntimeError);
  pgsl_error[30] = rb_define_class_under(mgsl_error, "ETOLX", rb_eRuntimeError);
  pgsl_error[31] = rb_define_class_under(mgsl_error, "ETOLG", rb_eRuntimeError);
  pgsl_error[32] = rb_define_class_under(mgsl_error, "EOF", rb_eEOFError);
}

static void Init_rb_gsl_define_GSL_CONST(VALUE module)
{
  rb_define_const(module, "SUCCESS", INT2FIX(GSL_SUCCESS));
  rb_define_const(module, "FAILURE", INT2FIX(GSL_FAILURE));
  rb_define_const(module, "CONTINUE", INT2FIX(GSL_CONTINUE));
  rb_define_const(module, "EDOM", INT2FIX(GSL_EDOM));
  rb_define_const(module, "ERANGE", INT2FIX(GSL_ERANGE));
  rb_define_const(module, "EFAULT", INT2FIX(GSL_EFAULT));
  rb_define_const(module, "EINVAL", INT2FIX(GSL_EINVAL));
  rb_define_const(module, "EFAILED", INT2FIX(GSL_EFAILED));
  rb_define_const(module, "EFACTOR", INT2FIX(GSL_EFACTOR));
  rb_define_const(module, "ESANITY", INT2FIX(GSL_ESANITY));
  rb_define_const(module, "ENOMEM", INT2FIX(GSL_ENOMEM));
  rb_define_const(module, "EBADFUNC", INT2FIX(GSL_EBADFUNC));
  rb_define_const(module, "ERUNAWAY", INT2FIX(GSL_ERUNAWAY));
  rb_define_const(module, "EMAXITER", INT2FIX(GSL_EMAXITER));
  rb_define_const(module, "EZERODIV", INT2FIX(GSL_EZERODIV));
  rb_define_const(module, "EBADTOL", INT2FIX(GSL_EBADTOL));
  rb_define_const(module, "ETOL", INT2FIX(GSL_ETOL));
  rb_define_const(module, "EUNDRFLW", INT2FIX(GSL_EUNDRFLW));
  rb_define_const(module, "EOVRFLW", INT2FIX(GSL_EOVRFLW));
  rb_define_const(module, "ELOSS", INT2FIX(GSL_ELOSS));
  rb_define_const(module, "EROUND", INT2FIX(GSL_EROUND));
  rb_define_const(module, "EBADLEN", INT2FIX(GSL_EBADLEN));
  rb_define_const(module, "ENOTSQR", INT2FIX(GSL_ENOTSQR));
  rb_define_const(module, "ESING", INT2FIX(GSL_ESING));
  rb_define_const(module, "EDIVERGE", INT2FIX(GSL_EDIVERGE));
  rb_define_const(module, "EUNSUP", INT2FIX(GSL_EUNSUP));
  rb_define_const(module, "EUNIMPL", INT2FIX(GSL_EUNIMPL));
  rb_define_const(module, "ECACHE", INT2FIX(GSL_ECACHE));
  rb_define_const(module, "ETABLE", INT2FIX(GSL_ETABLE));
  rb_define_const(module, "ENOPROG", INT2FIX(GSL_ENOPROG));
  rb_define_const(module, "ENOPROGJ", INT2FIX(GSL_ENOPROGJ));
  rb_define_const(module, "ETOLF", INT2FIX(GSL_ETOLF));
  rb_define_const(module, "ETOLX", INT2FIX(GSL_ETOLX));
  rb_define_const(module, "ETOLG", INT2FIX(GSL_ETOLG));
  rb_define_const(module, "EOF", INT2FIX(GSL_EOF));
}

static VALUE rb_gsl_set_error_handler_off(VALUE module)
{
  gsl_set_error_handler_off();
  return Qtrue;
}

static void define_module_functions(VALUE module);
static VALUE rb_gsl_strerror(VALUE obj, VALUE errn);
static void define_module_functions(VALUE module)
{
  rb_define_module_function(module, "set_error_handler_off",
          rb_gsl_set_error_handler_off, 0);
  rb_define_module_function(module, "strerror",
          rb_gsl_strerror, 1);
  rb_define_module_function(module, "set_error_handler",
          rb_gsl_set_error_handler, -1);
  rb_define_module_function(module, "set_default_error_handler",
          rb_gsl_set_default_error_handler, 0);
}

static VALUE rb_gsl_strerror(VALUE obj, VALUE errn)
{
  int gsl_errno = FIX2INT(errn);
  return rb_str_new2(gsl_strerror(gsl_errno));
}

void Init_gsl_error(VALUE module)
{
  Init_rb_gsl_define_GSL_CONST(module);

  gsl_set_error_handler(&rb_gsl_error_handler);

  define_module_functions(module);
  rb_gsl_define_exceptions(module);
}
