/*
  sf_ellint.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "include/rb_gsl_sf.h"

static VALUE rb_gsl_sf_ellint_Kcomp(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) return eval_sf(gsl_sf_ellint_Kcomp, argv[0]);
  else return rb_gsl_sf_eval_double_m(gsl_sf_ellint_Kcomp, argv[0], argv[1]);
}

static VALUE rb_gsl_sf_ellint_Kcomp_e(VALUE obj, VALUE k, VALUE m)
{
  return rb_gsl_sf_eval_e_m(gsl_sf_ellint_Kcomp_e, k, m);
}

static VALUE rb_gsl_sf_ellint_Ecomp(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 1) return eval_sf(gsl_sf_ellint_Ecomp, argv[0]);
  else return rb_gsl_sf_eval_double_m(gsl_sf_ellint_Ecomp, argv[0], argv[1]);
}

static VALUE rb_gsl_sf_ellint_Ecomp_e(VALUE obj, VALUE k, VALUE m)
{
  return rb_gsl_sf_eval_e_m(gsl_sf_ellint_Ecomp_e, k, m);
}

static VALUE rb_gsl_sf_ellint_F(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 2) 
    return rb_gsl_sf_eval_double2_m(gsl_sf_ellint_F, argv[0], argv[1], 
            INT2FIX(GSL_PREC_DOUBLE));
  else 
    return rb_gsl_sf_eval_double2_m(gsl_sf_ellint_F, argv[0], argv[1], argv[2]);
}

static VALUE rb_gsl_sf_ellint_F_e(VALUE obj, VALUE phi, VALUE k, VALUE m)
{
  return rb_gsl_sf_eval_e_double2_m(gsl_sf_ellint_F_e, phi, k, m);
}

static VALUE rb_gsl_sf_ellint_E(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 2) 
    return rb_gsl_sf_eval_double2_m(gsl_sf_ellint_E, argv[0], argv[1], 
            INT2FIX(GSL_PREC_DOUBLE));
  else 
    return rb_gsl_sf_eval_double2_m(gsl_sf_ellint_E, argv[0], argv[1], argv[2]);
}

static VALUE rb_gsl_sf_ellint_E_e(VALUE obj, VALUE phi, VALUE k, VALUE m)
{
  return rb_gsl_sf_eval_e_double2_m(gsl_sf_ellint_E_e, phi, k, m);
}

static VALUE rb_gsl_sf_ellint_P(int argc, VALUE *argv, VALUE obj)

{
  if (argc == 3) 
    return rb_gsl_sf_eval_double3_m(gsl_sf_ellint_P, argv[0], argv[1], argv[2],
            INT2FIX(GSL_PREC_DOUBLE));
  else 
    return rb_gsl_sf_eval_double3_m(gsl_sf_ellint_P, argv[0], argv[1], argv[2],
            argv[3]);
}

static VALUE rb_gsl_sf_ellint_P_e(VALUE obj, VALUE phi, VALUE k, 
          VALUE n, VALUE m)
{
  return rb_gsl_sf_eval_e_double3_m(gsl_sf_ellint_P_e, phi, k, n, m);
}

static VALUE rb_gsl_sf_ellint_D(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 3) 
    return rb_gsl_sf_eval_double3_m(gsl_sf_ellint_D, argv[0], argv[1], argv[2],
            INT2FIX(GSL_PREC_DOUBLE));
  else 
    return rb_gsl_sf_eval_double3_m(gsl_sf_ellint_D, argv[0], argv[1], argv[2],
            argv[3]);
}

static VALUE rb_gsl_sf_ellint_D_e(VALUE obj, VALUE phi, VALUE k, 
          VALUE n, VALUE m)
{
  return rb_gsl_sf_eval_e_double3_m(gsl_sf_ellint_D_e, phi, k, n, m);
}

static VALUE rb_gsl_sf_ellint_RC(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 2) 
    return rb_gsl_sf_eval_double2_m(gsl_sf_ellint_RC, argv[0], argv[1], 
            INT2FIX(GSL_PREC_DOUBLE));
  else 
    return rb_gsl_sf_eval_double2_m(gsl_sf_ellint_RC, argv[0], argv[1], argv[2]);
}

static VALUE rb_gsl_sf_ellint_RC_e(VALUE obj, VALUE x, VALUE y, VALUE m)
{
  return rb_gsl_sf_eval_e_double2_m(gsl_sf_ellint_RC_e, x, y, m);
}

static VALUE rb_gsl_sf_ellint_RD(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 3) 
    return rb_gsl_sf_eval_double3_m(gsl_sf_ellint_RD, argv[0], argv[1], argv[2],
            INT2FIX(GSL_PREC_DOUBLE));
  else 
    return rb_gsl_sf_eval_double3_m(gsl_sf_ellint_RD, argv[0], argv[1], argv[2],
            argv[3]);
}

static VALUE rb_gsl_sf_ellint_RD_e(VALUE obj, VALUE x, VALUE y, 
           VALUE z, VALUE m)
{
  return rb_gsl_sf_eval_e_double3_m(gsl_sf_ellint_RD_e, x, y, z, m);
}

static VALUE rb_gsl_sf_ellint_RF(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 3) 
    return rb_gsl_sf_eval_double3_m(gsl_sf_ellint_RF, argv[0], argv[1], argv[2],
            INT2FIX(GSL_PREC_DOUBLE));
  else 
    return rb_gsl_sf_eval_double3_m(gsl_sf_ellint_RF, argv[0], argv[1], argv[2],
            argv[3]);
}

static VALUE rb_gsl_sf_ellint_RF_e(VALUE obj, VALUE x, VALUE y, 
           VALUE z, VALUE m)
{
  return rb_gsl_sf_eval_e_double3_m(gsl_sf_ellint_RF_e, x, y, z, m);
}

static VALUE rb_gsl_sf_ellint_RJ(int argc, VALUE *argv, VALUE obj)
{
  if (argc == 4) 
    return rb_gsl_sf_eval_double4_m(gsl_sf_ellint_RJ, argv[0], argv[1], argv[2],
            argv[3], INT2FIX(GSL_PREC_DOUBLE));
  else 
    return rb_gsl_sf_eval_double4_m(gsl_sf_ellint_RJ, argv[0], argv[1], argv[2],
            argv[3], argv[4]);
}

static VALUE rb_gsl_sf_ellint_RJ_e(VALUE obj, VALUE x, VALUE y, 
           VALUE z, VALUE p, VALUE m)
{
  return rb_gsl_sf_eval_e_double4_m(gsl_sf_ellint_RJ_e, x, y, z, p, m);
}

void Init_gsl_sf_ellint(VALUE module)
{
  VALUE mgsl_sf_ellint;

  rb_define_module_function(module, "ellint_Kcomp",  rb_gsl_sf_ellint_Kcomp, -1);
  rb_define_module_function(module, "ellint_Kcomp_e",  rb_gsl_sf_ellint_Kcomp_e, 2);
  rb_define_module_function(module, "ellint_Ecomp",  rb_gsl_sf_ellint_Ecomp, -1);
  rb_define_module_function(module, "ellint_Ecomp_e",  rb_gsl_sf_ellint_Ecomp_e, 2);
  rb_define_module_function(module, "ellint_F",  rb_gsl_sf_ellint_F, -1);
  rb_define_module_function(module, "ellint_F_e",  rb_gsl_sf_ellint_F_e, 3);
  rb_define_module_function(module, "ellint_E",  rb_gsl_sf_ellint_E, -1);
  rb_define_module_function(module, "ellint_E_e",  rb_gsl_sf_ellint_E_e, 3);
  rb_define_module_function(module, "ellint_P",  rb_gsl_sf_ellint_P, -1);
  rb_define_module_function(module, "ellint_P_e",  rb_gsl_sf_ellint_P_e, 4);
  rb_define_module_function(module, "ellint_D",  rb_gsl_sf_ellint_D, -1);
  rb_define_module_function(module, "ellint_D_e",  rb_gsl_sf_ellint_D_e, 4);
  rb_define_module_function(module, "ellint_RC",  rb_gsl_sf_ellint_RC, -1);
  rb_define_module_function(module, "ellint_RC_e",  rb_gsl_sf_ellint_RC_e, 3);
  rb_define_module_function(module, "ellint_RD",  rb_gsl_sf_ellint_RD, -1);
  rb_define_module_function(module, "ellint_RD_e",  rb_gsl_sf_ellint_RD_e, 4);
  rb_define_module_function(module, "ellint_RF",  rb_gsl_sf_ellint_RF, -1);
  rb_define_module_function(module, "ellint_RF_e",  rb_gsl_sf_ellint_RF_e, 4);
  rb_define_module_function(module, "ellint_RJ",  rb_gsl_sf_ellint_RJ, -1);
  rb_define_module_function(module, "ellint_RJ_e",  rb_gsl_sf_ellint_RJ_e, 5);

  mgsl_sf_ellint = rb_define_module_under(module, "Ellint");
  rb_define_module_function(mgsl_sf_ellint, "Kcomp",  rb_gsl_sf_ellint_Kcomp, -1);
  rb_define_module_function(mgsl_sf_ellint, "Kcomp_e",  rb_gsl_sf_ellint_Kcomp_e, 2);
  rb_define_module_function(mgsl_sf_ellint, "Ecomp",  rb_gsl_sf_ellint_Ecomp, -1);
  rb_define_module_function(mgsl_sf_ellint, "Ecomp_e",  rb_gsl_sf_ellint_Ecomp_e, 2);
  rb_define_module_function(mgsl_sf_ellint, "F",  rb_gsl_sf_ellint_F, -1);
  rb_define_module_function(mgsl_sf_ellint, "F_e",  rb_gsl_sf_ellint_F_e, 3);
  rb_define_module_function(mgsl_sf_ellint, "E",  rb_gsl_sf_ellint_E, -1);
  rb_define_module_function(mgsl_sf_ellint, "E_e",  rb_gsl_sf_ellint_E_e, 3);
  rb_define_module_function(mgsl_sf_ellint, "P",  rb_gsl_sf_ellint_P, -1);
  rb_define_module_function(mgsl_sf_ellint, "P_e",  rb_gsl_sf_ellint_P_e, 4);
  rb_define_module_function(mgsl_sf_ellint, "D",  rb_gsl_sf_ellint_D, -1);
  rb_define_module_function(mgsl_sf_ellint, "D_e",  rb_gsl_sf_ellint_D_e, 4);
  rb_define_module_function(mgsl_sf_ellint, "RC",  rb_gsl_sf_ellint_RC, -1);
  rb_define_module_function(mgsl_sf_ellint, "RC_e",  rb_gsl_sf_ellint_RC_e, 3);
  rb_define_module_function(mgsl_sf_ellint, "RD",  rb_gsl_sf_ellint_RD, -1);
  rb_define_module_function(mgsl_sf_ellint, "RD_e",  rb_gsl_sf_ellint_RD_e, 4);
  rb_define_module_function(mgsl_sf_ellint, "RF",  rb_gsl_sf_ellint_RF, -1);
  rb_define_module_function(mgsl_sf_ellint, "RF_e",  rb_gsl_sf_ellint_RF_e, 4);
  rb_define_module_function(mgsl_sf_ellint, "RJ",  rb_gsl_sf_ellint_RJ, -1);
  rb_define_module_function(mgsl_sf_ellint, "RJ_e",  rb_gsl_sf_ellint_RJ_e, 5);
}
