/*
  sf_coupling.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "include/rb_gsl_sf.h"

static VALUE rb_gsl_sf_coupling_3j(VALUE obj, VALUE two_ja, VALUE two_jb, VALUE two_jc, VALUE two_ma, VALUE two_mb, VALUE two_mc) 
{
  CHECK_FIXNUM(two_ja); CHECK_FIXNUM(two_jb); CHECK_FIXNUM(two_jc);
  CHECK_FIXNUM(two_ma); CHECK_FIXNUM(two_mb); CHECK_FIXNUM(two_mc);
  return rb_float_new(gsl_sf_coupling_3j(FIX2INT(two_ja), FIX2INT(two_jb),
           FIX2INT(two_jc), FIX2INT(two_ma),
           FIX2INT(two_mb), FIX2INT(two_mc)));
}

static VALUE rb_gsl_sf_coupling_3j_e(VALUE obj, VALUE two_ja, VALUE two_jb, VALUE two_jc, VALUE two_ma, VALUE two_mb, VALUE two_mc) 
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  CHECK_FIXNUM(two_ja); CHECK_FIXNUM(two_jb); CHECK_FIXNUM(two_jc);
  CHECK_FIXNUM(two_ma); CHECK_FIXNUM(two_mb); CHECK_FIXNUM(two_mc);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_coupling_3j_e(FIX2INT(two_ja), FIX2INT(two_jb),
        FIX2INT(two_jc), FIX2INT(two_ma),
        FIX2INT(two_mb), FIX2INT(two_mc),
        rslt);

  return v;
}

static VALUE rb_gsl_sf_coupling_6j(VALUE obj, VALUE two_ja, VALUE two_jb, VALUE two_jc, VALUE two_jd, VALUE two_je, VALUE two_jf) 
{
  CHECK_FIXNUM(two_ja); CHECK_FIXNUM(two_jb); CHECK_FIXNUM(two_jc);
  CHECK_FIXNUM(two_jd); CHECK_FIXNUM(two_je); CHECK_FIXNUM(two_jf);
  return rb_float_new(gsl_sf_coupling_6j(FIX2INT(two_ja), FIX2INT(two_jb),
           FIX2INT(two_jc), FIX2INT(two_jd),
           FIX2INT(two_je), FIX2INT(two_jf)));
}

static VALUE rb_gsl_sf_coupling_6j_e(VALUE obj, VALUE two_ja, VALUE two_jb, VALUE two_jc, VALUE two_jd, VALUE two_je, VALUE two_jf) 
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  CHECK_FIXNUM(two_ja); CHECK_FIXNUM(two_jb); CHECK_FIXNUM(two_jc);
  CHECK_FIXNUM(two_jd); CHECK_FIXNUM(two_je); CHECK_FIXNUM(two_jf);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_coupling_6j_e(FIX2INT(two_ja), FIX2INT(two_jb),
        FIX2INT(two_jc), FIX2INT(two_jd),
        FIX2INT(two_je), FIX2INT(two_jf),
        rslt);
  return v;
}

static VALUE rb_gsl_sf_coupling_9j(VALUE obj, VALUE two_ja, VALUE two_jb,
           VALUE two_jc, VALUE two_jd, VALUE two_je,
           VALUE two_jf, VALUE two_jg, VALUE two_jh,
           VALUE two_ji) 
{
  CHECK_FIXNUM(two_ja); CHECK_FIXNUM(two_jb); CHECK_FIXNUM(two_jc);
  CHECK_FIXNUM(two_jd); CHECK_FIXNUM(two_je); CHECK_FIXNUM(two_jf);
  CHECK_FIXNUM(two_jg); CHECK_FIXNUM(two_jh); CHECK_FIXNUM(two_ji);
  return rb_float_new(gsl_sf_coupling_9j(FIX2INT(two_ja), FIX2INT(two_jb),
          FIX2INT(two_jc), FIX2INT(two_jd),
          FIX2INT(two_je), FIX2INT(two_jf),
          FIX2INT(two_jg), FIX2INT(two_jh),
          FIX2INT(two_ji)));
}

static VALUE rb_gsl_sf_coupling_9j_e(VALUE obj, VALUE two_ja, VALUE two_jb,
           VALUE two_jc, VALUE two_jd, VALUE two_je,
           VALUE two_jf, VALUE two_jg, VALUE two_jh,
           VALUE two_ji) 
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  CHECK_FIXNUM(two_ja); CHECK_FIXNUM(two_jb); CHECK_FIXNUM(two_jc);
  CHECK_FIXNUM(two_jd); CHECK_FIXNUM(two_je); CHECK_FIXNUM(two_jf);
  CHECK_FIXNUM(two_jg); CHECK_FIXNUM(two_jh); CHECK_FIXNUM(two_ji);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_coupling_9j_e(FIX2INT(two_ja), FIX2INT(two_jb),
         FIX2INT(two_jc), FIX2INT(two_jd),
         FIX2INT(two_je), FIX2INT(two_jf),
         FIX2INT(two_jg), FIX2INT(two_jh),
         FIX2INT(two_ji), rslt);
  return v;
}

void Init_gsl_sf_coupling(VALUE module)
{
  VALUE mgsl_sf_coupling;

  rb_define_module_function(module, "coupling_3j",  rb_gsl_sf_coupling_3j, 6);
  rb_define_module_function(module, "coupling_3j_e",  rb_gsl_sf_coupling_3j_e, 6);
  rb_define_module_function(module, "coupling_6j",  rb_gsl_sf_coupling_6j, 6);
  rb_define_module_function(module, "coupling_6j_e",  rb_gsl_sf_coupling_6j_e, 6);
  rb_define_module_function(module, "coupling_9j",  rb_gsl_sf_coupling_9j, 9);
  rb_define_module_function(module, "coupling_9j_e",  rb_gsl_sf_coupling_9j_e, 9);

  mgsl_sf_coupling = rb_define_module_under(module, "Coupling");
  
  rb_define_module_function(mgsl_sf_coupling, "3j",  rb_gsl_sf_coupling_3j, 6);
  rb_define_module_function(mgsl_sf_coupling, "3j_e",  rb_gsl_sf_coupling_3j_e, 6);
  rb_define_module_function(mgsl_sf_coupling, "6j",  rb_gsl_sf_coupling_6j, 6);
  rb_define_module_function(mgsl_sf_coupling, "6j_e",  rb_gsl_sf_coupling_6j_e, 6);
  rb_define_module_function(mgsl_sf_coupling, "9j",  rb_gsl_sf_coupling_9j, 9);
  rb_define_module_function(mgsl_sf_coupling, "9j_e",  rb_gsl_sf_coupling_9j_e, 9);
}
