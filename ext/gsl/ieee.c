/*
  ieee.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl.h"

static VALUE rb_gsl_ieee_env_setup(VALUE obj)
{
  gsl_ieee_env_setup();
  return obj;
}

static VALUE rb_gsl_ieee_fprintf_double(int argc, VALUE *argv, VALUE obj)
{
#ifdef HAVE_RUBY_IO_H
  rb_io_t *fptr = NULL;
#else
  OpenFile *fptr = NULL;
#endif
  FILE *fp = NULL;
  int flag = 0;
  VALUE vtmp;
  double ftmp;

  switch (argc) {
  case 2:
    switch (TYPE(argv[0])) {
    case T_STRING:
      fp = fopen(RSTRING_PTR(argv[0]), "w");
      flag = 1;
      break;
    case T_FILE:
      GetOpenFile(argv[0], fptr);
      rb_io_check_writable(fptr);
#ifdef HAVE_RUBY_IO_H
      fp = rb_io_stdio_file(fptr);
#else
      fp = GetWriteFile(fptr);
#endif
      break;
    default:
      rb_raise(rb_eTypeError, "wrong type argument %s (IO or String expected)",
         rb_class2name(CLASS_OF(argv[0])));
    }
    vtmp = argv[1];
    break;
  case 1:
    vtmp = argv[0];
    fp = stdout;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  }
  if (TYPE(vtmp) != T_FLOAT)
    rb_raise(rb_eTypeError, "wrong argument type %s (Float expected)",
       rb_class2name(CLASS_OF(vtmp)));
  ftmp = RFLOAT_VALUE(vtmp);
  gsl_ieee_fprintf_double(fp, &ftmp);
  if (fp == stdout) fprintf(stdout, "\n");
  if (flag == 1) fclose(fp);
  return obj;
}

static VALUE rb_gsl_ieee_printf_double(VALUE obj, VALUE xx)
{
  double x;
  x = NUM2DBL(xx);
  gsl_ieee_printf_double(&x);
  return xx;
}

void Init_gsl_ieee(VALUE module)
{
  VALUE mgsl_ieee;
  mgsl_ieee = rb_define_module_under(module, "IEEE");

  rb_define_singleton_method(mgsl_ieee, "env_setup",
           rb_gsl_ieee_env_setup, 0);
  rb_define_module_function(module, "ieee_env_setup", rb_gsl_ieee_env_setup, 0);
  rb_define_singleton_method(mgsl_ieee, "fprintf_double",
           rb_gsl_ieee_fprintf_double, -1);
  rb_define_singleton_method(mgsl_ieee, "fprintf",
           rb_gsl_ieee_fprintf_double, -1);
  rb_define_singleton_method(mgsl_ieee, "printf",
           rb_gsl_ieee_printf_double, -1);
  rb_define_singleton_method(mgsl_ieee, "printf_double",
           rb_gsl_ieee_printf_double, -1);

}
