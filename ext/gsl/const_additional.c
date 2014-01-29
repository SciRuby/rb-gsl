/*
  const_additional.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_const.h"

#ifndef GSL_1_4_9_LATER
#define MKS_STEFAN_BOLTZMANN_CONSTANT (5.6703993443638e-08)
#define MKS_THOMSON_CROSS_SECTION (6.65245852869236e-29)
#define CGS_STEFAN_BOLTZMANN_CONSTANT (5.6703993443638e-05)
#define CGS_THOMSON_CROSS_SECTION (6.65245852869236e-25)
#endif

#define MKS_CLASSICAL_ELECTRON_RADIUS (2.81794028310825e-15)
#define MKS_RADIATION_DENSITY_CONSTANT (7.56576650685962e-16)
#define MKS_RADIATION_NUMBER_DENSITY_CONSTANT (20.2869161066108e-6)
#define MKS_SOLAR_TIME (4.925490947e-6)
#define MKS_SOLAR_GM (1.3271243999e20)
#define MKS_PLANCK_MASS (2.17664610503472e-08)
#define MKS_PLANCK_LENGTH (1.61609744261067e-35)
#define MKS_PLANCK_TIME (5.39072081196475e-44)

#define CGS_CLASSICAL_ELECTRON_RADIUS (2.81794028310825e-13)
#define CGS_RADIATION_DENSITY_CONSTANT (7.56576650685962e-15)
#define CGS_RADIATION_NUMBER_DENSITY_CONSTANT (20.2869161066108)
#define CGS_SOLAR_TIME (4.925490947e-6)
#define CGS_SOLAR_GM (1.3271243999e26)
#define CGS_PLANCK_MASS (2.17664610503472e-05)
#define CGS_PLANCK_LENGTH (1.61609744261067e-33)
#define CGS_PLANCK_TIME (5.39072081196475e-44)

static void rb_gsl_const_mks(VALUE module);
static void rb_gsl_const_cgs(VALUE module);
static void rb_gsl_const_num(VALUE module);

static void rb_gsl_const_mks(VALUE module)
{
  rb_define_const(module, "RADIATION_DENSITY_CONSTANT",
		  rb_float_new(MKS_RADIATION_DENSITY_CONSTANT));
  rb_define_const(module, "RADIATION_NUMBER_DENSITY_CONSTANT", 
		  rb_float_new(MKS_RADIATION_NUMBER_DENSITY_CONSTANT));
  rb_define_const(module, "CLASSICAL_ELECTRON_RADIUS",
		  rb_float_new(MKS_CLASSICAL_ELECTRON_RADIUS));
  rb_define_const(module, "SOLAR_TIME", rb_float_new(MKS_SOLAR_TIME));
  rb_define_const(module, "SOLAR_GM", rb_float_new(MKS_SOLAR_GM));

  rb_define_const(module, "PLANCK_MASS", rb_float_new(MKS_PLANCK_MASS));
  rb_define_const(module, "PLANCK_LENGTH", rb_float_new(MKS_PLANCK_LENGTH));
  rb_define_const(module, "PLANCK_TIME", rb_float_new(MKS_PLANCK_TIME));

#ifndef GSL_1_4_9_LATER
  rb_define_const(module, "STEFAN_BOLTZMANN_CONSTANT", 
		  rb_float_new(MKS_STEFAN_BOLTZMANN_CONSTANT));
  rb_define_const(module, "THOMSON_CROSS_SECTION", 
		  rb_float_new(MKS_THOMSON_CROSS_SECTION));
#endif
}

static void rb_gsl_const_cgs(VALUE module)
{
  rb_define_const(module, "RADIATION_DENSITY_CONSTANT", 
		  rb_float_new(CGS_RADIATION_DENSITY_CONSTANT));
  rb_define_const(module, "RADIATION_NUMBER_DENSITY_CONSTANT", 
		  rb_float_new(CGS_RADIATION_NUMBER_DENSITY_CONSTANT));
  rb_define_const(module, "CLASSICAL_ELECTRON_RADIUS", 
		  rb_float_new(CGS_CLASSICAL_ELECTRON_RADIUS));
  rb_define_const(module, "SOLAR_TIME", rb_float_new(CGS_SOLAR_TIME));
  rb_define_const(module, "SOLAR_GM", rb_float_new(CGS_SOLAR_GM));

  rb_define_const(module, "PLANCK_MASS", rb_float_new(CGS_PLANCK_MASS));
  rb_define_const(module, "PLANCK_LENGTH", rb_float_new(CGS_PLANCK_LENGTH));
  rb_define_const(module, "PLANCK_TIME", rb_float_new(CGS_PLANCK_TIME));

#ifndef GSL_1_4_9_LATER
  rb_define_const(module, "STEFAN_BOLTZMANN_CONSTANT", 
		  rb_float_new(CGS_STEFAN_BOLTZMANN_CONSTANT));
  rb_define_const(module, "THOMSON_CROSS_SECTION", 
		  rb_float_new(CGS_THOMSON_CROSS_SECTION));
#endif
}

static void rb_gsl_const_num(VALUE module)
{

}

void Init_gsl_const_additional(VALUE mmks, VALUE mcgs, VALUE mnum)
{
  rb_gsl_const_mks(mmks);
  rb_gsl_const_cgs(mcgs);
  rb_gsl_const_num(mnum);
}

#undef MKS_CLASSICAL_ELECTRON_RADIUS
#undef MKS_STEFAN_BOLTZMANN_CONSTANT
#undef MKS_RADIATION_DENSITY_CONSTANT
#undef MKS_RADIATION_NUMBER_DENSITY_CONSTANT
#undef CGS_CLASSICAL_ELECTRON_RADIUS 
#undef CGS_STEFAN_BOLTZMANN_CONSTANT 
#undef CGS_RADIATION_DENSITY_CONSTANT
#undef CGS_RADIATION_NUMBER_DENSITY_CONSTANT
#undef CGS_THOMSON_CROSS_SECTION
#undef MKS_THOMSON_CROSS_SECTION
#undef MKS_SOLAR_TIME
#undef CGS_SOLAR_TIME
#undef MKS_SOLAR_GM
#undef CGS_SOLAR_GM
#undef MKS_PLANCK_MASS
#undef MKS_PLANCK_LENGTH
#undef MKS_PLANCK_TIME
#undef CGS_PLANCK_MASS
#undef CGS_PLANCK_LENGTH
#undef CGS_PLANCK_TIME
