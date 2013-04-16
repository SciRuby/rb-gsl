/*
  sf_coulomb.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "rb_gsl_sf.h"
EXTERN VALUE cgsl_vector;

static VALUE rb_gsl_sf_hydrogenicR_1(VALUE obj, VALUE Z, VALUE r)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_hydrogenicR_1, Z, r);
}

static VALUE rb_gsl_sf_hydrogenicR_1_e(VALUE obj,  VALUE Z, VALUE r)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_hydrogenicR_1_e, Z, r);
}

static VALUE rb_gsl_sf_hydrogenicR(VALUE obj, VALUE n, VALUE l, 
				   VALUE Z, VALUE r)
{
  return rb_float_new(gsl_sf_hydrogenicR(FIX2INT(n), FIX2INT(l),
					   NUM2DBL(Z), NUM2DBL(r)));
}

static VALUE rb_gsl_sf_hydrogenicR_e(VALUE obj, VALUE n, VALUE l, 
				     VALUE Z, VALUE r)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  // local variable "status" declared and set, but never used
  //int status;
  CHECK_FIXNUM(n); CHECK_FIXNUM(l);
  Need_Float(Z); Need_Float(r);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  /*status =*/ gsl_sf_hydrogenicR_e(FIX2INT(n), FIX2INT(l),
			    NUM2DBL(Z), NUM2DBL(r), rslt);
  return v;
}

static VALUE rb_gsl_sf_coulomb_wave_FG_e(VALUE obj, VALUE eta, VALUE x,
					 VALUE L_F, VALUE k)

{
  gsl_sf_result *F, *Fp, *G, *Gp;
  VALUE vF, vFp, vG, vGp;
  double exp_G, exp_F;
  int status;
  Need_Float(eta); Need_Float(x); Need_Float(L_F);
  CHECK_FIXNUM(k);
  vF = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, F);
  vFp = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, Fp);
  vG = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, G);
  vGp = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, Gp);
  status = gsl_sf_coulomb_wave_FG_e(NUM2DBL(eta), NUM2DBL(x), NUM2DBL(L_F),
				    FIX2INT(k), F, Fp, G, Gp, &exp_F, &exp_G);
  return rb_ary_new3(7, vF, vFp, vG, vGp,
		     rb_float_new(exp_F), rb_float_new(exp_G), INT2FIX(status));
}

static VALUE rb_gsl_sf_coulomb_wave_F_array(VALUE obj, VALUE Lmin, VALUE kmax,
					 VALUE eta, VALUE x)
{
  double F_exponent;
  int status;
  size_t size;
  gsl_vector *v = NULL;
  CHECK_FIXNUM(kmax);
  Need_Float(Lmin); Need_Float(eta); Need_Float(x);
  size = FIX2INT(kmax);
  v = gsl_vector_alloc(size);
  status = gsl_sf_coulomb_wave_F_array(NUM2DBL(Lmin), size, NUM2DBL(eta), 
				   NUM2DBL(x), v->data, &F_exponent);
  
  return rb_ary_new3(3, Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v), 
		     rb_float_new(F_exponent), INT2FIX(status));
}

static VALUE rb_gsl_sf_coulomb_wave_FG_array(VALUE obj, VALUE Lmin, VALUE kmax,
					 VALUE eta, VALUE x)
{
  double F_exponent, G_exponent;
  int status;
  size_t size;
  gsl_vector *vf = NULL, *vg = NULL;
  VALUE fary, gary;
  CHECK_FIXNUM(kmax);
  Need_Float(Lmin); Need_Float(eta); Need_Float(x);
  size = FIX2INT(kmax);
  vf = gsl_vector_alloc(size);
  vg = gsl_vector_alloc(size);

  status = gsl_sf_coulomb_wave_FG_array(NUM2DBL(Lmin), size, NUM2DBL(eta), 
				   NUM2DBL(x), vf->data, vg->data,
				    &F_exponent, &G_exponent);
  fary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vf);
  gary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vf);
  return rb_ary_new3(5, fary, gary,
		     rb_float_new(F_exponent), rb_float_new(G_exponent),
		     INT2FIX(status));
}

static VALUE rb_gsl_sf_coulomb_wave_FGp_array(VALUE obj, VALUE Lmin, VALUE kmax,
					 VALUE eta, VALUE x)
{
  double F_exponent, G_exponent;
  int status;
  size_t size;
  gsl_vector *vf = NULL, *vg = NULL, *vfp = NULL, *vgp = NULL;
  VALUE fary, gary, fpary, gpary;
  CHECK_FIXNUM(kmax);
  Need_Float(Lmin); Need_Float(eta); Need_Float(x);
  size = FIX2INT(kmax);
  vf = gsl_vector_alloc(size);
  vfp = gsl_vector_alloc(size);
  vg = gsl_vector_alloc(size);
  vgp = gsl_vector_alloc(size);

  status = gsl_sf_coulomb_wave_FGp_array(NUM2DBL(Lmin), size, NUM2DBL(eta), 
				     NUM2DBL(x), vf->data, vfp->data, 
				     vg->data, vgp->data,
				     &F_exponent, &G_exponent);
  fary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vf);
  fpary =Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vfp);
  gary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vg);
  gpary =Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vgp);
  return rb_ary_new3(7, fary, fpary, gary, gpary,
		     rb_float_new(F_exponent), rb_float_new(G_exponent),
		     INT2FIX(status));
}

static VALUE rb_gsl_sf_coulomb_wave_sphF_array(VALUE obj, VALUE Lmin, VALUE kmax,
					 VALUE eta, VALUE x)
{
  int status;
  size_t size;
  gsl_vector *v = NULL, *v2 = NULL;
  CHECK_FIXNUM(kmax);
  Need_Float(Lmin); Need_Float(eta); Need_Float(x);
  size = FIX2INT(kmax);
  v = gsl_vector_alloc(size);
  v2 = gsl_vector_alloc(size);
  status =  gsl_sf_coulomb_wave_sphF_array(NUM2DBL(Lmin), size, NUM2DBL(eta), 
				       NUM2DBL(x), v->data, v2->data);
  return rb_ary_new3(3, Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v), 
		     Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v2), 
		     INT2FIX(status));
}

static VALUE rb_gsl_sf_coulomb_CL_e(VALUE obj, VALUE L, VALUE eta)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_coulomb_CL_e, L, eta);
}

static VALUE rb_gsl_sf_coulomb_CL_array(VALUE obj, VALUE Lmin, VALUE kmax,
					VALUE eta)
{
  gsl_vector *v = NULL;
  size_t size;
  // local variable "status" declared and set, but never used
  //int status;
  CHECK_FIXNUM(kmax);
  Need_Float(Lmin); Need_Float(eta); 
  size = FIX2INT(kmax);
  v = gsl_vector_alloc(size);
  /*status =*/ gsl_sf_coulomb_CL_array(NUM2DBL(Lmin), size, NUM2DBL(eta), v->data);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

void Init_gsl_sf_coulomb(VALUE module)
{
  VALUE mgsl_sf_coulomb;

  rb_define_module_function(module, "hydrogenicR_1",  rb_gsl_sf_hydrogenicR_1, 2);
  rb_define_module_function(module, "hydrogenicR_1_e",  rb_gsl_sf_hydrogenicR_1_e, 2);
  rb_define_module_function(module, "hydrogenicR",  rb_gsl_sf_hydrogenicR, 4);
  rb_define_module_function(module, "hydrogenicR_e",  rb_gsl_sf_hydrogenicR_e, 4);
  rb_define_module_function(module, "coulomb_wave_FG_e",  rb_gsl_sf_coulomb_wave_FG_e, 4);
  rb_define_module_function(module, "coulomb_wave_F_array",  rb_gsl_sf_coulomb_wave_F_array, 4);
  rb_define_module_function(module, "coulomb_wave_FG_array",  rb_gsl_sf_coulomb_wave_FG_array, 4);
  rb_define_module_function(module, "coulomb_wave_FGp_array",  rb_gsl_sf_coulomb_wave_FGp_array, 4);
  rb_define_module_function(module, "coulomb_wave_sphF_array ",  rb_gsl_sf_coulomb_wave_sphF_array, 4);
  rb_define_module_function(module, "coulomb_CL_e",  rb_gsl_sf_coulomb_CL_e, 2);
  rb_define_module_function(module, "coulomb_CL_array",  rb_gsl_sf_coulomb_CL_array, 3);

  mgsl_sf_coulomb = rb_define_module_under(module, "Coulomb");
  
  rb_define_module_function(mgsl_sf_coulomb, "hydrogenicR_1",  rb_gsl_sf_hydrogenicR_1, 2);
  rb_define_module_function(mgsl_sf_coulomb, "hydrogenicR_1_e",  rb_gsl_sf_hydrogenicR_1_e, 2);
  rb_define_module_function(mgsl_sf_coulomb, "hydrogenicR",  rb_gsl_sf_hydrogenicR, 4);
  rb_define_module_function(mgsl_sf_coulomb, "hydrogenicR_e",  rb_gsl_sf_hydrogenicR_e, 4);
  rb_define_module_function(mgsl_sf_coulomb, "wave_FG_e",  rb_gsl_sf_coulomb_wave_FG_e, 4);
  rb_define_module_function(mgsl_sf_coulomb, "wave_F_array",  rb_gsl_sf_coulomb_wave_F_array, 4);
  rb_define_module_function(mgsl_sf_coulomb, "wave_FG_array",  rb_gsl_sf_coulomb_wave_FG_array, 4);
  rb_define_module_function(mgsl_sf_coulomb, "wave_FGp_array",  rb_gsl_sf_coulomb_wave_FGp_array, 4);
  rb_define_module_function(mgsl_sf_coulomb, "wave_sphF_array ",  rb_gsl_sf_coulomb_wave_sphF_array, 4);
  rb_define_module_function(mgsl_sf_coulomb, "CL_e",  rb_gsl_sf_coulomb_CL_e, 2);
  rb_define_module_function(mgsl_sf_coulomb, "CL_array",  rb_gsl_sf_coulomb_CL_array, 3);
}
