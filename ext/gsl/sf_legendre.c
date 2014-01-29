/*
  legendre.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_sf.h"
EXTERN VALUE cgsl_vector;

static VALUE rb_gsl_sf_legendre_P1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_legendre_P1, x);
}

static VALUE rb_gsl_sf_legendre_P1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_legendre_P1_e, x);
}

static VALUE rb_gsl_sf_legendre_P2(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_legendre_P2, x);
}

static VALUE rb_gsl_sf_legendre_P2_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_legendre_P2_e, x);
}

static VALUE rb_gsl_sf_legendre_P3(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_legendre_P3, x);
}

static VALUE rb_gsl_sf_legendre_P3_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_legendre_P3_e, x);
}

static VALUE rb_gsl_sf_legendre_Pl(VALUE obj, VALUE l, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_legendre_Pl, l, x);
}

static VALUE rb_gsl_sf_legendre_Pl_e(VALUE obj, VALUE l, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_legendre_Pl_e, l, x);
}

static VALUE rb_gsl_sf_legendre_Pl_array(VALUE obj, VALUE lmax, VALUE x)
{
  gsl_vector *v = NULL;
  CHECK_FIXNUM(lmax);
  Need_Float(x);
  v = gsl_vector_alloc(FIX2INT(lmax) + 1);
  gsl_sf_legendre_Pl_array(FIX2INT(lmax), NUM2DBL(x), v->data);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

static VALUE rb_gsl_sf_legendre_Q0(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_legendre_Q0, x);
}

static VALUE rb_gsl_sf_legendre_Q0_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_legendre_Q0_e, x);
}

static VALUE rb_gsl_sf_legendre_Q1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_legendre_Q1, x);
}

static VALUE rb_gsl_sf_legendre_Q1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_legendre_Q1_e, x);
}

static VALUE rb_gsl_sf_legendre_Ql(VALUE obj, VALUE l, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_legendre_Ql, l, x);
}

static VALUE rb_gsl_sf_legendre_Ql_e(VALUE obj, VALUE l, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_legendre_Ql_e, l, x);
}

static VALUE rb_gsl_sf_legendre_Plm(VALUE obj, VALUE l, VALUE m, VALUE x)
{
  return rb_gsl_sf_eval_int_int_double(gsl_sf_legendre_Plm, l, m, x);
}

static VALUE rb_gsl_sf_legendre_Plm_e(VALUE obj, VALUE l, VALUE m, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  int status;
  CHECK_FIXNUM(l);
  CHECK_FIXNUM(m);
  Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  status = gsl_sf_legendre_Plm_e(FIX2INT(l), FIX2INT(m), NUM2DBL(x), rslt);
  return rb_ary_new3(2, v, INT2FIX(status));
}

static VALUE rb_gsl_sf_legendre_Plm_array(VALUE obj, VALUE lmax, VALUE m, VALUE x)
{
  gsl_vector *v = NULL;
  int size;
  int ll, mm;
  CHECK_FIXNUM(lmax); CHECK_FIXNUM(m);
  Need_Float(x);
  ll = FIX2INT(lmax);
  mm = FIX2INT(m);
  size = gsl_sf_legendre_array_size(ll, mm);
  v = gsl_vector_alloc(size);
  gsl_sf_legendre_Plm_array(ll, mm, NUM2DBL(x), v->data);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

static VALUE rb_gsl_sf_legendre_sphPlm(VALUE obj, VALUE l, VALUE m, VALUE x)
{
  return rb_gsl_sf_eval_int_int_double(gsl_sf_legendre_sphPlm, l, m, x);
}

static VALUE rb_gsl_sf_legendre_sphPlm_e(VALUE obj, VALUE l, VALUE m, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  int status;
  CHECK_FIXNUM(l); CHECK_FIXNUM(m);
  Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  status = gsl_sf_legendre_sphPlm_e(FIX2INT(l), FIX2INT(m), NUM2DBL(x), rslt);
  return rb_ary_new3(2, v, INT2FIX(status));
}

static VALUE rb_gsl_sf_legendre_sphPlm_array(VALUE obj, VALUE lmax, VALUE m, VALUE x)
{
  gsl_vector *v = NULL;
  int size;
  int ll, mm;
  CHECK_FIXNUM(lmax); CHECK_FIXNUM(m);
  Need_Float(x);
  ll = FIX2INT(lmax);
  mm = FIX2INT(m);
  size = gsl_sf_legendre_array_size(ll, mm);
  v = gsl_vector_alloc(size);
  gsl_sf_legendre_sphPlm_array(ll, mm, NUM2DBL(x), v->data);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

static VALUE rb_gsl_sf_legendre_array_size(VALUE obj, VALUE lmax, VALUE m)
{
  CHECK_FIXNUM(lmax);   CHECK_FIXNUM(m);
  return INT2FIX(gsl_sf_legendre_array_size(FIX2INT(lmax), FIX2INT(m)));
}

static VALUE rb_gsl_sf_conicalP_half(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_conicalP_half, lambda, x);
}

static VALUE rb_gsl_sf_conicalP_half_e(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_conicalP_half_e, lambda, x);
}

static VALUE rb_gsl_sf_conicalP_mhalf(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_conicalP_mhalf, lambda, x);
}

static VALUE rb_gsl_sf_conicalP_mhalf_e(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_conicalP_mhalf_e, lambda, x);
}

static VALUE rb_gsl_sf_conicalP_0(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_conicalP_0, lambda, x);
}

static VALUE rb_gsl_sf_conicalP_0_e(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_conicalP_0_e, lambda, x);
}

static VALUE rb_gsl_sf_conicalP_1(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_conicalP_1, lambda, x);
}

static VALUE rb_gsl_sf_conicalP_1_e(VALUE obj, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_conicalP_1_e, lambda, x);
}

static VALUE rb_gsl_sf_conicalP_sph_reg(VALUE obj, VALUE l, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_int_double_double(gsl_sf_conicalP_sph_reg, l, lambda, x);
}

static VALUE rb_gsl_sf_conicalP_sph_reg_e(VALUE obj, VALUE l, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double2(gsl_sf_conicalP_sph_reg_e, l, lambda, x);
}

static VALUE rb_gsl_sf_conicalP_cyl_reg(VALUE obj, VALUE m, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_int_double_double(gsl_sf_conicalP_cyl_reg, m, lambda, x);
}

static VALUE rb_gsl_sf_conicalP_cyl_reg_e(VALUE obj, VALUE m, VALUE lambda, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double2(gsl_sf_conicalP_cyl_reg_e, m, lambda, x);
}

static VALUE rb_gsl_sf_legendre_H3d_0(VALUE obj, VALUE lambda, VALUE eta)
{
  Need_Float(lambda); Need_Float(eta);
  return rb_float_new(gsl_sf_legendre_H3d_0(NUM2DBL(lambda), NUM2DBL(eta)));
}

static VALUE rb_gsl_sf_legendre_H3d_0_e(VALUE obj, VALUE lambda, VALUE eta)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_legendre_H3d_0_e, lambda, eta);
}

static VALUE rb_gsl_sf_legendre_H3d_1(VALUE obj, VALUE lambda, VALUE eta)
{
  Need_Float(lambda); Need_Float(eta);
  return rb_float_new(gsl_sf_legendre_H3d_1(NUM2DBL(lambda), NUM2DBL(eta)));
}

static VALUE rb_gsl_sf_legendre_H3d_1_e(VALUE obj, VALUE lambda, VALUE eta)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_legendre_H3d_1_e, lambda, eta);
}


static VALUE rb_gsl_sf_legendre_H3d(VALUE obj, VALUE l, VALUE lambda, VALUE eta)
{
  return rb_float_new(gsl_sf_legendre_H3d(FIX2INT(l), NUM2DBL(lambda), NUM2DBL(eta)));
}

static VALUE rb_gsl_sf_legendre_H3d_e(VALUE obj,VALUE l,  VALUE lambda, VALUE eta)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  CHECK_FIXNUM(l);
  Need_Float(lambda);
  Need_Float(eta);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_legendre_H3d_e(FIX2INT(l), NUM2DBL(lambda), NUM2DBL(eta), rslt);
  return v;
}

static VALUE rb_gsl_sf_legendre_H3d_array(VALUE obj, VALUE lmax, VALUE lambda, VALUE eta)
{
  gsl_vector *v = NULL;
  CHECK_FIXNUM(lmax);
  Need_Float(lambda);
  Need_Float(eta);
  v = gsl_vector_alloc(FIX2INT(lmax) + 1);
  gsl_sf_legendre_H3d_array(FIX2INT(lmax), NUM2DBL(lambda), NUM2DBL(eta), v->data);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

void Init_gsl_sf_legendre(VALUE module)
{
  VALUE mgsl_sf_leg;

  rb_define_module_function(module, "legendre_P1",  rb_gsl_sf_legendre_P1, 1);
  rb_define_module_function(module, "legendre_P1_e",  rb_gsl_sf_legendre_P1_e, 1);
  rb_define_module_function(module, "legendre_P2",  rb_gsl_sf_legendre_P2, 1);
  rb_define_module_function(module, "legendre_P2_e",  rb_gsl_sf_legendre_P2_e, 1);
  rb_define_module_function(module, "legendre_P3",  rb_gsl_sf_legendre_P3, 1);
  rb_define_module_function(module, "legendre_P3_e",  rb_gsl_sf_legendre_P3_e, 1);
  rb_define_module_function(module, "legendre_Pl",  rb_gsl_sf_legendre_Pl, 2);
  rb_define_module_function(module, "legendre_Pl_e",  rb_gsl_sf_legendre_Pl_e, 2);
  rb_define_module_function(module, "legendre_Pl_array",  rb_gsl_sf_legendre_Pl_array, 2);
  rb_define_module_function(module, "legendre_Q0",  rb_gsl_sf_legendre_Q0, 1);
  rb_define_module_function(module, "legendre_Q0_e",  rb_gsl_sf_legendre_Q0_e, 1);
  rb_define_module_function(module, "legendre_Q1",  rb_gsl_sf_legendre_Q1, 1);
  rb_define_module_function(module, "legendre_Q1_e",  rb_gsl_sf_legendre_Q1_e, 1);
  rb_define_module_function(module, "legendre_Ql",  rb_gsl_sf_legendre_Ql, 2);
  rb_define_module_function(module, "legendre_Ql_e",  rb_gsl_sf_legendre_Ql_e, 2);
  rb_define_module_function(module, "legendre_Plm",  rb_gsl_sf_legendre_Plm, 3);
  rb_define_module_function(module, "legendre_Plm_e",  rb_gsl_sf_legendre_Plm_e, 3);
  rb_define_module_function(module, "legendre_Plm_array",  rb_gsl_sf_legendre_Plm_array, 3);
  rb_define_module_function(module, "legendre_sphPlm",  rb_gsl_sf_legendre_sphPlm, 3);
  rb_define_module_function(module, "legendre_sphPlm_e",  rb_gsl_sf_legendre_sphPlm_e, 3);
  rb_define_module_function(module, "legendre_sphPlm_array",  rb_gsl_sf_legendre_sphPlm_array, 3);
  rb_define_module_function(module, "legendre_array_size",  rb_gsl_sf_legendre_array_size, 2);
  rb_define_module_function(module, "conicalP_half",  rb_gsl_sf_conicalP_half, 2);
  rb_define_module_function(module, "conicalP_half_e",  rb_gsl_sf_conicalP_half_e, 2);
  rb_define_module_function(module, "conicalP_mhalf",  rb_gsl_sf_conicalP_mhalf, 2);
  rb_define_module_function(module, "conicalP_mhalf_e",  rb_gsl_sf_conicalP_mhalf_e, 2);
  rb_define_module_function(module, "conicalP_0",  rb_gsl_sf_conicalP_0, 2);
  rb_define_module_function(module, "conicalP_0_e",  rb_gsl_sf_conicalP_0_e, 2);
  rb_define_module_function(module, "conicalP_1",  rb_gsl_sf_conicalP_1, 2);
  rb_define_module_function(module, "conicalP_1_e",  rb_gsl_sf_conicalP_1_e, 2);
  rb_define_module_function(module, "conicalP_sph_reg",  rb_gsl_sf_conicalP_sph_reg, 3);
  rb_define_module_function(module, "conicalP_sph_reg_e",  rb_gsl_sf_conicalP_sph_reg_e, 3);
  rb_define_module_function(module, "conicalP_cyl_reg",  rb_gsl_sf_conicalP_cyl_reg, 3);
  rb_define_module_function(module, "conicalP_cyl_reg_e",  rb_gsl_sf_conicalP_cyl_reg_e, 3);
  rb_define_module_function(module, "legendre_H3d_0",  rb_gsl_sf_legendre_H3d_0, 2);
  rb_define_module_function(module, "legendre_H3d_0_e",  rb_gsl_sf_legendre_H3d_0_e, 2);
  rb_define_module_function(module, "legendre_H3d_1",  rb_gsl_sf_legendre_H3d_1, 2);
  rb_define_module_function(module, "legendre_H3d_1_e",  rb_gsl_sf_legendre_H3d_1_e, 2);
  rb_define_module_function(module, "legendre_H3d",  rb_gsl_sf_legendre_H3d, 3);
  rb_define_module_function(module, "legendre_H3d_e",  rb_gsl_sf_legendre_H3d_e, 3);
  rb_define_module_function(module, "legendre_H3d_array",  rb_gsl_sf_legendre_H3d_array, 3);

  /*****/

  mgsl_sf_leg = rb_define_module_under(module, "Legendre");
  rb_define_module_function(mgsl_sf_leg, "P1",  rb_gsl_sf_legendre_P1, 1);
  rb_define_module_function(mgsl_sf_leg, "P1_e",  rb_gsl_sf_legendre_P1_e, 1);
  rb_define_module_function(mgsl_sf_leg, "P2",  rb_gsl_sf_legendre_P2, 1);
  rb_define_module_function(mgsl_sf_leg, "P2_e",  rb_gsl_sf_legendre_P2_e, 1);
  rb_define_module_function(mgsl_sf_leg, "P3",  rb_gsl_sf_legendre_P3, 1);
  rb_define_module_function(mgsl_sf_leg, "P3_e",  rb_gsl_sf_legendre_P3_e, 1);
  rb_define_module_function(mgsl_sf_leg, "Pl",  rb_gsl_sf_legendre_Pl, 2);
  rb_define_module_function(mgsl_sf_leg, "Pl_e",  rb_gsl_sf_legendre_Pl_e, 2);
  rb_define_module_function(mgsl_sf_leg, "Pl_array",  rb_gsl_sf_legendre_Pl_array, 2);
  rb_define_module_function(mgsl_sf_leg, "Q0",  rb_gsl_sf_legendre_Q0, 1);
  rb_define_module_function(mgsl_sf_leg, "Q0_e",  rb_gsl_sf_legendre_Q0_e, 1);
  rb_define_module_function(mgsl_sf_leg, "Q1",  rb_gsl_sf_legendre_Q1, 1);
  rb_define_module_function(mgsl_sf_leg, "Q1_e",  rb_gsl_sf_legendre_Q1_e, 1);
  rb_define_module_function(mgsl_sf_leg, "Plm",  rb_gsl_sf_legendre_Plm, 3);
  rb_define_module_function(mgsl_sf_leg, "Plm_e",  rb_gsl_sf_legendre_Plm_e, 3);
  rb_define_module_function(mgsl_sf_leg, "Plm_array",  rb_gsl_sf_legendre_Plm_array, 3);
  rb_define_module_function(mgsl_sf_leg, "sphPlm",  rb_gsl_sf_legendre_sphPlm, 3);
  rb_define_module_function(mgsl_sf_leg, "sphPlm_e",  rb_gsl_sf_legendre_sphPlm_e, 3);
  rb_define_module_function(mgsl_sf_leg, "sphPlm_array",  rb_gsl_sf_legendre_sphPlm_array, 3);
  rb_define_module_function(mgsl_sf_leg, "array_size",  rb_gsl_sf_legendre_array_size, 2);
  rb_define_module_function(mgsl_sf_leg, "conicalP_half",  rb_gsl_sf_conicalP_half, 2);
  rb_define_module_function(mgsl_sf_leg, "conicalP_half_e",  rb_gsl_sf_conicalP_half_e, 2);
  rb_define_module_function(mgsl_sf_leg, "conicalP_mhalf",  rb_gsl_sf_conicalP_mhalf, 2);
  rb_define_module_function(mgsl_sf_leg, "conicalP_mhalf_e",  rb_gsl_sf_conicalP_mhalf_e, 2);
  rb_define_module_function(mgsl_sf_leg, "conicalP_0",  rb_gsl_sf_conicalP_0, 2);
  rb_define_module_function(mgsl_sf_leg, "conicalP_0_e",  rb_gsl_sf_conicalP_0_e, 2);
  rb_define_module_function(mgsl_sf_leg, "conicalP_1",  rb_gsl_sf_conicalP_1, 2);
  rb_define_module_function(mgsl_sf_leg, "conicalP_1_e",  rb_gsl_sf_conicalP_1_e, 2);
  rb_define_module_function(mgsl_sf_leg, "conicalP_sph_reg",  rb_gsl_sf_conicalP_sph_reg, 3);
  rb_define_module_function(mgsl_sf_leg, "conicalP_sph_reg_e",  rb_gsl_sf_conicalP_sph_reg_e, 3);
  rb_define_module_function(mgsl_sf_leg, "conicalP_cyl_reg",  rb_gsl_sf_conicalP_cyl_reg, 3);
  rb_define_module_function(mgsl_sf_leg, "conicalP_cyl_reg_e",  rb_gsl_sf_conicalP_cyl_reg_e, 3);
  rb_define_module_function(mgsl_sf_leg, "H3d_0",  rb_gsl_sf_legendre_H3d_0, 2);
  rb_define_module_function(mgsl_sf_leg, "H3d_0_e",  rb_gsl_sf_legendre_H3d_0_e, 2);
  rb_define_module_function(mgsl_sf_leg, "H3d_1",  rb_gsl_sf_legendre_H3d_1, 2);
  rb_define_module_function(mgsl_sf_leg, "H3d_1_e",  rb_gsl_sf_legendre_H3d_1_e, 2);
  rb_define_module_function(mgsl_sf_leg, "H3d",  rb_gsl_sf_legendre_H3d, 3);
  rb_define_module_function(mgsl_sf_leg, "H3d_e",  rb_gsl_sf_legendre_H3d_e, 3);
  rb_define_module_function(mgsl_sf_leg, "H3d_array",  rb_gsl_sf_legendre_H3d_array, 3);

}
