/*
  const.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_const.h"

static void rb_gsl_const_mks(VALUE module);
static void rb_gsl_const_cgs(VALUE module);
static void rb_gsl_const_num(VALUE module);
void Init_gsl_const_additional(VALUE mmks, VALUE mcgs, VALUE mnum);

#ifndef GSL_1_4_LATER
static void rb_gsl_const_mks(VALUE module)
{
  rb_define_const(module, "SPEED_OF_LIGHT",
      rb_float_new(GSL_CONST_MKS_SPEED_OF_LIGHT));
  rb_define_const(module, "GRAVITATIONAL_CONSTANT",
      rb_float_new(GSL_CONST_MKS_GRAVITATIONAL_CONSTANT));
  rb_define_const(module, "PLANCKS_CONSTANT_H",
      rb_float_new(GSL_CONST_MKS_PLANCKS_CONSTANT_H));
  rb_define_const(module, "PLANCKS_CONSTANT_HBAR",
      rb_float_new(GSL_CONST_MKS_PLANCKS_CONSTANT_HBAR));
  rb_define_const(module, "VACUUM_PERMEABILITY",
      rb_float_new(GSL_CONST_MKS_VACUUM_PERMEABILITY));
  rb_define_const(module, "ASTRONOMICAL_UNIT",
      rb_float_new(GSL_CONST_MKS_ASTRONOMICAL_UNIT));
  rb_define_const(module, "LIGHT_YEAR", rb_float_new(GSL_CONST_MKS_LIGHT_YEAR));
  rb_define_const(module, "PARSEC", rb_float_new(GSL_CONST_MKS_PARSEC));
  rb_define_const(module, "GRAV_ACCEL", rb_float_new(GSL_CONST_MKS_GRAV_ACCEL));
  rb_define_const(module, "ELECTRON_VOLT",
      rb_float_new(GSL_CONST_MKS_ELECTRON_VOLT));
  rb_define_const(module, "MASS_ELECTRON",
      rb_float_new(GSL_CONST_MKS_MASS_ELECTRON));
  rb_define_const(module, "MASS_MUON", rb_float_new(GSL_CONST_MKS_MASS_MUON));
  rb_define_const(module, "MASS_PROTON", rb_float_new(GSL_CONST_MKS_MASS_PROTON));
  rb_define_const(module, "MASS_NEUTRON", rb_float_new(GSL_CONST_MKS_MASS_NEUTRON));
  rb_define_const(module, "RYDBERG", rb_float_new(GSL_CONST_MKS_RYDBERG));
  rb_define_const(module, "BOLTZMANN", rb_float_new(GSL_CONST_MKS_BOLTZMANN));
  rb_define_const(module, "MOLAR_GAS", rb_float_new(GSL_CONST_MKS_MOLAR_GAS));
  rb_define_const(module, "BOHR_MAGNETON",
      rb_float_new(GSL_CONST_MKS_BOHR_MAGNETON));
  rb_define_const(module, "NUCLEAR_MAGNETON",
      rb_float_new(GSL_CONST_MKS_NUCLEAR_MAGNETON));
  rb_define_const(module, "ELECTRON_MAGNETIC_MOMENT",
      rb_float_new(GSL_CONST_MKS_ELECTRON_MAGNETIC_MOMENT));
  rb_define_const(module, "PROTON_MAGNETIC_MOMENT",
      rb_float_new(GSL_CONST_MKS_PROTON_MAGNETIC_MOMENT));
  rb_define_const(module, "STANDARD_GAS_VOLUME",
      rb_float_new(GSL_CONST_MKS_STANDARD_GAS_VOLUME));

  rb_define_const(module, "MINUTE", rb_float_new(GSL_CONST_MKS_MINUTE));
  rb_define_const(module, "HOUR", rb_float_new(GSL_CONST_MKS_HOUR));
  rb_define_const(module, "DAY", rb_float_new(GSL_CONST_MKS_DAY));
  rb_define_const(module, "WEEK", rb_float_new(GSL_CONST_MKS_WEEK));
  rb_define_const(module, "INCH", rb_float_new(GSL_CONST_MKS_INCH));
  rb_define_const(module, "FOOT", rb_float_new(GSL_CONST_MKS_FOOT));
  rb_define_const(module, "YARD", rb_float_new(GSL_CONST_MKS_YARD));
  rb_define_const(module, "MILE", rb_float_new(GSL_CONST_MKS_MILE));
  rb_define_const(module, "NAUTICAL_MILE", rb_float_new(GSL_CONST_MKS_NAUTICAL_MILE));
  rb_define_const(module, "FATHOM", rb_float_new(GSL_CONST_MKS_FATHOM));
  rb_define_const(module, "MIL", rb_float_new(GSL_CONST_MKS_MIL));
  rb_define_const(module, "POINT", rb_float_new(GSL_CONST_MKS_POINT));
  rb_define_const(module, "TEXPOINT", rb_float_new(GSL_CONST_MKS_TEXPOINT));
  rb_define_const(module, "MICRON", rb_float_new(GSL_CONST_MKS_MICRON));
  rb_define_const(module, "ANGSTROM", rb_float_new(GSL_CONST_MKS_ANGSTROM));
  rb_define_const(module, "HECTARE", rb_float_new(GSL_CONST_MKS_HECTARE));
  rb_define_const(module, "ACRE", rb_float_new(GSL_CONST_MKS_ACRE));
#ifdef GSL_0_9_4_LATER
  rb_define_const(module, "BARN", rb_float_new(GSL_CONST_MKS_BARN));
  rb_define_const(module, "BTU", rb_float_new(GSL_CONST_MKS_BTU));
  rb_define_const(module, "SOLAR_MASS", rb_float_new(GSL_CONST_MKS_SOLAR_MASS));
#else
  rb_define_const(module, "BARN", rb_float_new(1.0e-28));
  rb_define_const(module, "BTU", rb_float_new(1.05505585262e3));
#endif
  rb_define_const(module, "LITER", rb_float_new(GSL_CONST_MKS_LITER));
  rb_define_const(module, "US_GALLON", rb_float_new(GSL_CONST_MKS_US_GALLON));
  rb_define_const(module, "QUART", rb_float_new(GSL_CONST_MKS_QUART));
  rb_define_const(module, "PINT", rb_float_new(GSL_CONST_MKS_PINT));
  rb_define_const(module, "CUP", rb_float_new(GSL_CONST_MKS_CUP));
  rb_define_const(module, "FLUID_OUNCE", rb_float_new(GSL_CONST_MKS_FLUID_OUNCE));
  rb_define_const(module, "TABLESPOON", rb_float_new(GSL_CONST_MKS_TABLESPOON));
  rb_define_const(module, "CANADIAN_GALLON",
      rb_float_new(GSL_CONST_MKS_CANADIAN_GALLON));

  rb_define_const(module, "UK_GALLON", rb_float_new(GSL_CONST_MKS_UK_GALLON));
  rb_define_const(module, "KILOMETERS_PER_HOUR",
      rb_float_new(GSL_CONST_MKS_MILES_PER_HOUR));
  rb_define_const(module, "MILES_PER_HOUR",
      rb_float_new(GSL_CONST_MKS_KILOMETERS_PER_HOUR));
  rb_define_const(module, "KNOT", rb_float_new(GSL_CONST_MKS_KNOT));
  rb_define_const(module, "POUND_MASS", rb_float_new(GSL_CONST_MKS_POUND_MASS));
  rb_define_const(module, "POUND_OUNCE", rb_float_new(GSL_CONST_MKS_OUNCE_MASS));
  rb_define_const(module, "POUND_TON", rb_float_new(GSL_CONST_MKS_TON));
  rb_define_const(module, "POUND_METRIC_TON",
      rb_float_new(GSL_CONST_MKS_METRIC_TON));
  rb_define_const(module, "POUND_UK_TON", rb_float_new(GSL_CONST_MKS_UK_TON));
  rb_define_const(module, "POUND_TROY_OUNCE",
      rb_float_new(GSL_CONST_MKS_TROY_OUNCE));
  rb_define_const(module, "CARAT", rb_float_new(GSL_CONST_MKS_CARAT));
  rb_define_const(module, "UNIFIED_ATOMIC_MASS",
      rb_float_new(GSL_CONST_MKS_UNIFIED_ATOMIC_MASS));
  rb_define_const(module, "GRAM_FORCE", rb_float_new(GSL_CONST_MKS_GRAM_FORCE));
  rb_define_const(module, "POUND_FORCE", rb_float_new(GSL_CONST_MKS_POUND_FORCE));
  rb_define_const(module, "KILOPOUND_FORCE",
      rb_float_new(GSL_CONST_MKS_KILOPOUND_FORCE));
  rb_define_const(module, "POUNDAL", rb_float_new(GSL_CONST_MKS_POUNDAL));
  rb_define_const(module, "CALORIE", rb_float_new(GSL_CONST_MKS_CALORIE));
  rb_define_const(module, "THERM", rb_float_new(GSL_CONST_MKS_THERM));
  rb_define_const(module, "HORSEPOWER", rb_float_new(GSL_CONST_MKS_HORSEPOWER));
  rb_define_const(module, "BAR", rb_float_new(GSL_CONST_MKS_BAR));
  rb_define_const(module, "STD_ATMOSPHERE",
      rb_float_new(GSL_CONST_MKS_STD_ATMOSPHERE));
  rb_define_const(module, "TORR", rb_float_new(GSL_CONST_MKS_TORR));
  rb_define_const(module, "METER_OF_MERCURY",
      rb_float_new(GSL_CONST_MKS_METER_OF_MERCURY));
  rb_define_const(module, "INCH_OF_MERCURY",
      rb_float_new(GSL_CONST_MKS_INCH_OF_MERCURY));
  rb_define_const(module, "INCH_OF_WATER",
      rb_float_new(GSL_CONST_MKS_INCH_OF_WATER));
  rb_define_const(module, "PSI", rb_float_new(GSL_CONST_MKS_PSI));
  rb_define_const(module, "POISE", rb_float_new(GSL_CONST_MKS_POISE));
  rb_define_const(module, "STOKES", rb_float_new(GSL_CONST_MKS_STOKES));
  rb_define_const(module, "FARADAY", rb_float_new(GSL_CONST_MKS_FARADAY));
  rb_define_const(module, "ELECTRON_CHARGE",
      rb_float_new(GSL_CONST_MKS_ELECTRON_CHARGE));
  rb_define_const(module, "GAUSS", rb_float_new(GSL_CONST_MKS_GAUSS));
  rb_define_const(module, "STILB", rb_float_new(GSL_CONST_MKS_STILB));
  rb_define_const(module, "LUMEN", rb_float_new(GSL_CONST_MKS_LUMEN));
  rb_define_const(module, "LUX", rb_float_new(GSL_CONST_MKS_LUX));
  rb_define_const(module, "PHOT", rb_float_new(GSL_CONST_MKS_PHOT));
  rb_define_const(module, "FOOTCANDLE", rb_float_new(GSL_CONST_MKS_FOOTCANDLE));
  rb_define_const(module, "LAMBERT", rb_float_new(GSL_CONST_MKS_LAMBERT));
  rb_define_const(module, "CURIE", rb_float_new(GSL_CONST_MKS_CURIE));
  rb_define_const(module, "ROENTGEN", rb_float_new(GSL_CONST_MKS_ROENTGEN));
  rb_define_const(module, "RAD", rb_float_new(GSL_CONST_MKS_RAD));

#ifdef GSL_1_1_LATER
  rb_define_const(module, "VACUUM_PERMITTIVITY",
      rb_float_new(GSL_CONST_MKS_VACUUM_PERMITTIVITY));
  rb_define_const(module, "BOHR_RADIUS", rb_float_new(GSL_CONST_MKS_BOHR_RADIUS));
#endif

#ifdef GSL_1_2_LATER
  rb_define_const(module, "NEWTON", rb_float_new(GSL_CONST_MKS_NEWTON));
  rb_define_const(module, "DYNE", rb_float_new(GSL_CONST_MKS_DYNE));
  rb_define_const(module, "JOULE", rb_float_new(GSL_CONST_MKS_JOULE));
  rb_define_const(module, "ERG", rb_float_new(GSL_CONST_MKS_ERG));
#endif

#ifdef GSL_1_8_LATER
  rb_define_const(module, "DEBYE", rb_float_new(GSL_CONST_MKS_DEBYE));
#endif

}

static void rb_gsl_const_cgs(VALUE module)
{
  rb_define_const(module, "SPEED_OF_LIGHT",
      rb_float_new(GSL_CONST_CGS_SPEED_OF_LIGHT));
  rb_define_const(module, "GRAVITATIONAL_CONSTANT",
      rb_float_new(GSL_CONST_CGS_GRAVITATIONAL_CONSTANT));
  rb_define_const(module, "PLANCKS_CONSTANT_H",
      rb_float_new(GSL_CONST_CGS_PLANCKS_CONSTANT_H));
  rb_define_const(module, "PLANCKS_CONSTANT_HBAR",
      rb_float_new(GSL_CONST_CGS_PLANCKS_CONSTANT_HBAR));
  rb_define_const(module, "ASTRONOMICAL_UNIT",
      rb_float_new(GSL_CONST_CGS_ASTRONOMICAL_UNIT));
  rb_define_const(module, "LIGHT_YEAR", rb_float_new(GSL_CONST_CGS_LIGHT_YEAR));
  rb_define_const(module, "PARSEC", rb_float_new(GSL_CONST_CGS_PARSEC));
  rb_define_const(module, "GRAV_ACCEL", rb_float_new(GSL_CONST_CGS_GRAV_ACCEL));
  rb_define_const(module, "ELECTRON_VOLT",
      rb_float_new(GSL_CONST_CGS_ELECTRON_VOLT));
  rb_define_const(module, "MASS_ELECTRON",
      rb_float_new(GSL_CONST_CGS_MASS_ELECTRON));
  rb_define_const(module, "MASS_MUON", rb_float_new(GSL_CONST_CGS_MASS_MUON));
  rb_define_const(module, "MASS_PROTON",
      rb_float_new(GSL_CONST_CGS_MASS_PROTON));
  rb_define_const(module, "MASS_NEUTRON", rb_float_new(GSL_CONST_CGS_MASS_NEUTRON));
  rb_define_const(module, "RYDBERG", rb_float_new(GSL_CONST_CGS_RYDBERG));
  rb_define_const(module, "BOHR_MAGNETON",
      rb_float_new(GSL_CONST_CGS_BOHR_MAGNETON));

  rb_define_const(module, "NUCLEAR_MAGNETON",
      rb_float_new(GSL_CONST_CGS_NUCLEAR_MAGNETON));
  rb_define_const(module, "ELECTRON_MAGNETIC_MOMENT",
      rb_float_new(GSL_CONST_CGS_ELECTRON_MAGNETIC_MOMENT));
  rb_define_const(module, "PROTON_MAGNETIC_MOMENT",
      rb_float_new(GSL_CONST_CGS_PROTON_MAGNETIC_MOMENT));
  rb_define_const(module, "STANDARD_GAS_VOLUME",
      rb_float_new(GSL_CONST_CGS_STANDARD_GAS_VOLUME));

  rb_define_const(module, "MINUTE", rb_float_new(GSL_CONST_CGS_MINUTE));
  rb_define_const(module, "HOUR", rb_float_new(GSL_CONST_CGS_HOUR));
  rb_define_const(module, "DAY", rb_float_new(GSL_CONST_CGS_DAY));
  rb_define_const(module, "WEEK", rb_float_new(GSL_CONST_CGS_WEEK));
  rb_define_const(module, "INCH", rb_float_new(GSL_CONST_CGS_INCH));
  rb_define_const(module, "FOOT", rb_float_new(GSL_CONST_CGS_FOOT));
  rb_define_const(module, "YARD", rb_float_new(GSL_CONST_CGS_YARD));
  rb_define_const(module, "MILE", rb_float_new(GSL_CONST_CGS_MILE));
  rb_define_const(module, "NAUTICAL_MILE",
      rb_float_new(GSL_CONST_CGS_NAUTICAL_MILE));
  rb_define_const(module, "FATHOM", rb_float_new(GSL_CONST_CGS_FATHOM));
  rb_define_const(module, "MIL", rb_float_new(GSL_CONST_CGS_MIL));
  rb_define_const(module, "POINT", rb_float_new(GSL_CONST_CGS_POINT));
  rb_define_const(module, "TEXPOINT", rb_float_new(GSL_CONST_CGS_TEXPOINT));
  rb_define_const(module, "MICRON", rb_float_new(GSL_CONST_CGS_MICRON));
  rb_define_const(module, "ANGSTROM", rb_float_new(GSL_CONST_CGS_ANGSTROM));
  rb_define_const(module, "HECTARE", rb_float_new(GSL_CONST_CGS_HECTARE));
  rb_define_const(module, "ACRE", rb_float_new(GSL_CONST_CGS_ACRE));
#ifdef GSL_0_9_4_LATER
  rb_define_const(module, "BARN", rb_float_new(GSL_CONST_CGS_BARN));
  rb_define_const(module, "BTU", rb_float_new(GSL_CONST_CGS_BTU));
  rb_define_const(module, "SOLAR_MASS", rb_float_new(GSL_CONST_CGS_SOLAR_MASS));
#else
  rb_define_const(module, "BARN", rb_float_new(1e-24));
  rb_define_const(module, "BTU", rb_float_new(1.05505585262e10));
#endif
  rb_define_const(module, "LITER", rb_float_new(GSL_CONST_CGS_LITER));
  rb_define_const(module, "US_GALLON", rb_float_new(GSL_CONST_CGS_US_GALLON));
  rb_define_const(module, "QUART", rb_float_new(GSL_CONST_CGS_QUART));
  rb_define_const(module, "PINT", rb_float_new(GSL_CONST_CGS_PINT));
  rb_define_const(module, "CUP", rb_float_new(GSL_CONST_CGS_CUP));
  rb_define_const(module, "FLUID_OUNCE", rb_float_new(GSL_CONST_CGS_FLUID_OUNCE));
  rb_define_const(module, "TABLESPOON", rb_float_new(GSL_CONST_CGS_TABLESPOON));
  rb_define_const(module, "CANADIAN_GALLON",
      rb_float_new(GSL_CONST_CGS_CANADIAN_GALLON));

  rb_define_const(module, "UK_GALLON", rb_float_new(GSL_CONST_CGS_UK_GALLON));
  rb_define_const(module, "KILOMETERS_PER_HOUR",
      rb_float_new(GSL_CONST_CGS_MILES_PER_HOUR));
  rb_define_const(module, "MILES_PER_HOUR",
      rb_float_new(GSL_CONST_CGS_KILOMETERS_PER_HOUR));
  rb_define_const(module, "KNOT", rb_float_new(GSL_CONST_CGS_KNOT));
  rb_define_const(module, "POUND_MASS", rb_float_new(GSL_CONST_CGS_POUND_MASS));
  rb_define_const(module, "POUND_OUNCE", rb_float_new(GSL_CONST_CGS_OUNCE_MASS));
  rb_define_const(module, "POUND_TON", rb_float_new(GSL_CONST_CGS_TON));
  rb_define_const(module, "POUND_METRIC_TON",
      rb_float_new(GSL_CONST_CGS_METRIC_TON));
  rb_define_const(module, "POUND_UK_TON", rb_float_new(GSL_CONST_CGS_UK_TON));
  rb_define_const(module, "POUND_TROY_OUNCE",
      rb_float_new(GSL_CONST_CGS_TROY_OUNCE));
  rb_define_const(module, "CARAT", rb_float_new(GSL_CONST_CGS_CARAT));
  rb_define_const(module, "UNIFIED_ATOMIC_MASS",
      rb_float_new(GSL_CONST_CGS_UNIFIED_ATOMIC_MASS));
  rb_define_const(module, "GRAM_FORCE", rb_float_new(GSL_CONST_CGS_GRAM_FORCE));
  rb_define_const(module, "POUND_FORCE", rb_float_new(GSL_CONST_CGS_POUND_FORCE));
  rb_define_const(module, "KILOPOUND_FORCE",
      rb_float_new(GSL_CONST_CGS_KILOPOUND_FORCE));
  rb_define_const(module, "POUNDAL", rb_float_new(GSL_CONST_CGS_POUNDAL));
  rb_define_const(module, "CALORIE", rb_float_new(GSL_CONST_CGS_CALORIE));
  rb_define_const(module, "THERM", rb_float_new(GSL_CONST_CGS_THERM));
  rb_define_const(module, "HORSEPOWER", rb_float_new(GSL_CONST_CGS_HORSEPOWER));
  rb_define_const(module, "BAR", rb_float_new(GSL_CONST_CGS_BAR));
  rb_define_const(module, "STD_ATMOSPHERE",
      rb_float_new(GSL_CONST_CGS_STD_ATMOSPHERE));
  rb_define_const(module, "TORR", rb_float_new(GSL_CONST_CGS_TORR));
  rb_define_const(module, "METER_OF_MERCURY",
      rb_float_new(GSL_CONST_CGS_METER_OF_MERCURY));
  rb_define_const(module, "INCH_OF_MERCURY",
      rb_float_new(GSL_CONST_CGS_INCH_OF_MERCURY));
  rb_define_const(module, "INCH_OF_WATER",
      rb_float_new(GSL_CONST_CGS_INCH_OF_WATER));
  rb_define_const(module, "PSI", rb_float_new(GSL_CONST_CGS_PSI));
  rb_define_const(module, "POISE", rb_float_new(GSL_CONST_CGS_POISE));
  rb_define_const(module, "STOKES", rb_float_new(GSL_CONST_CGS_STOKES));
  rb_define_const(module, "FARADAY", rb_float_new(GSL_CONST_CGS_FARADAY));
  rb_define_const(module, "ELECTRON_CHARGE",
      rb_float_new(GSL_CONST_CGS_ELECTRON_CHARGE));
  rb_define_const(module, "ELECTRON_CHARGE_ESU",
      rb_float_new(GSL_CONST_CGS_ELECTRON_CHARGE*GSL_CONST_CGS_SPEED_OF_LIGHT));
  rb_define_const(module, "GAUSS", rb_float_new(GSL_CONST_CGS_GAUSS));
  rb_define_const(module, "STILB", rb_float_new(GSL_CONST_CGS_STILB));
  rb_define_const(module, "LUMEN", rb_float_new(GSL_CONST_CGS_LUMEN));
  rb_define_const(module, "LUX", rb_float_new(GSL_CONST_CGS_LUX));
  rb_define_const(module, "PHOT", rb_float_new(GSL_CONST_CGS_PHOT));
  rb_define_const(module, "FOOTCANDLE", rb_float_new(GSL_CONST_CGS_FOOTCANDLE));
  rb_define_const(module, "LAMBERT", rb_float_new(GSL_CONST_CGS_LAMBERT));
  rb_define_const(module, "CURIE", rb_float_new(GSL_CONST_CGS_CURIE));
  rb_define_const(module, "ROENTGEN", rb_float_new(GSL_CONST_CGS_ROENTGEN));
  rb_define_const(module, "RAD", rb_float_new(GSL_CONST_CGS_RAD));

  rb_define_const(module, "BOLTZMANN", rb_float_new(GSL_CONST_CGS_BOLTZMANN));
  rb_define_const(module, "MOLAR_GAS", rb_float_new(GSL_CONST_CGS_MOLAR_GAS));

#ifdef GSL_1_1_LATER
  rb_define_const(module, "BOHR_RADIUS", rb_float_new(GSL_CONST_CGS_BOHR_RADIUS));
#endif
#ifdef GSL_1_2_LATER
  rb_define_const(module, "NEWTON", rb_float_new(GSL_CONST_CGS_NEWTON));
  rb_define_const(module, "DYNE", rb_float_new(GSL_CONST_CGS_DYNE));
  rb_define_const(module, "JOULE", rb_float_new(GSL_CONST_CGS_JOULE));
  rb_define_const(module, "ERG", rb_float_new(GSL_CONST_CGS_ERG));
#endif

}

static void rb_gsl_const_num(VALUE module)
{
  rb_define_const(module, "AVOGADRO", rb_float_new(GSL_CONST_NUM_AVOGADRO));
  rb_define_const(module, "FINE_STRUCTURE",
      rb_float_new(GSL_CONST_NUM_FINE_STRUCTURE));
#ifdef GSL_1_2_LATER
  rb_define_const(module, "YOTTA", rb_float_new(GSL_CONST_NUM_YOTTA));
  rb_define_const(module, "ZETTA", rb_float_new(GSL_CONST_NUM_ZETTA));
  rb_define_const(module, "EXA", rb_float_new(GSL_CONST_NUM_EXA));
  rb_define_const(module, "PETA", rb_float_new(GSL_CONST_NUM_PETA));
  rb_define_const(module, "TERA", rb_float_new(GSL_CONST_NUM_TERA));
  rb_define_const(module, "GIGA", rb_float_new(GSL_CONST_NUM_GIGA));
  rb_define_const(module, "MEGA", rb_float_new(GSL_CONST_NUM_MEGA));
  rb_define_const(module, "KILO", rb_float_new(GSL_CONST_NUM_KILO));
  rb_define_const(module, "MILLI", rb_float_new(GSL_CONST_NUM_MILLI));
  rb_define_const(module, "MICRO", rb_float_new(GSL_CONST_NUM_MICRO));
  rb_define_const(module, "NANO", rb_float_new(GSL_CONST_NUM_NANO));
  rb_define_const(module, "PICO", rb_float_new(GSL_CONST_NUM_PICO));
  rb_define_const(module, "FEMTO", rb_float_new(GSL_CONST_NUM_FEMTO));
  rb_define_const(module, "ATTO", rb_float_new(GSL_CONST_NUM_ATTO));
  rb_define_const(module, "ZEPTO", rb_float_new(GSL_CONST_NUM_ZEPTO));
  rb_define_const(module, "YOCTO", rb_float_new(GSL_CONST_NUM_YOCTO));
#endif
}

void Init_gsl_const(VALUE module)
{
  VALUE mgsl_const;
  VALUE mgsl_const_mks, mgsl_const_cgs, mgsl_const_num;

  mgsl_const = rb_define_module_under(module, "CONST");
  mgsl_const_mks = rb_define_module_under(mgsl_const, "MKSA");
  rb_gsl_const_mks(mgsl_const_mks);
  mgsl_const_cgs = rb_define_module_under(mgsl_const, "CGSM");
  rb_gsl_const_cgs(mgsl_const_cgs);
  mgsl_const_num = rb_define_module_under(mgsl_const, "NUM");
  rb_gsl_const_num(mgsl_const_num);
  Init_gsl_const_additional(mgsl_const_mks, mgsl_const_cgs, mgsl_const_num);
}

#else

static void rb_gsl_const_mks(VALUE module)
{
  rb_define_const(module, "SPEED_OF_LIGHT",
      rb_float_new(GSL_CONST_MKSA_SPEED_OF_LIGHT));
  rb_define_const(module, "GRAVITATIONAL_CONSTANT",
      rb_float_new(GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT));
  rb_define_const(module, "PLANCKS_CONSTANT_H",
      rb_float_new(GSL_CONST_MKSA_PLANCKS_CONSTANT_H));
  rb_define_const(module, "PLANCKS_CONSTANT_HBAR",
      rb_float_new(GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR));
  rb_define_const(module, "VACUUM_PERMEABILITY",
      rb_float_new(GSL_CONST_MKSA_VACUUM_PERMEABILITY));
  rb_define_const(module, "VACUUM_PERMITTIVITY",
      rb_float_new(GSL_CONST_MKSA_VACUUM_PERMITTIVITY));
  rb_define_const(module, "ASTRONOMICAL_UNIT",
      rb_float_new(GSL_CONST_MKSA_ASTRONOMICAL_UNIT));
  rb_define_const(module, "LIGHT_YEAR", rb_float_new(GSL_CONST_MKSA_LIGHT_YEAR));
  rb_define_const(module, "PARSEC", rb_float_new(GSL_CONST_MKSA_PARSEC));
  rb_define_const(module, "GRAV_ACCEL", rb_float_new(GSL_CONST_MKSA_GRAV_ACCEL));
  rb_define_const(module, "ELECTRON_VOLT",
      rb_float_new(GSL_CONST_MKSA_ELECTRON_VOLT));
  rb_define_const(module, "MASS_ELECTRON",
      rb_float_new(GSL_CONST_MKSA_MASS_ELECTRON));
  rb_define_const(module, "MASS_MUON", rb_float_new(GSL_CONST_MKSA_MASS_MUON));
  rb_define_const(module, "MASS_PROTON", rb_float_new(GSL_CONST_MKSA_MASS_PROTON));
  rb_define_const(module, "MASS_NEUTRON", rb_float_new(GSL_CONST_MKSA_MASS_NEUTRON));
  rb_define_const(module, "RYDBERG", rb_float_new(GSL_CONST_MKSA_RYDBERG));
  rb_define_const(module, "BOHR_MAGNETON",
      rb_float_new(GSL_CONST_MKSA_BOHR_MAGNETON));
  rb_define_const(module, "NUCLEAR_MAGNETON",
      rb_float_new(GSL_CONST_MKSA_NUCLEAR_MAGNETON));
  rb_define_const(module, "ELECTRON_MAGNETIC_MOMENT",
      rb_float_new(GSL_CONST_MKSA_ELECTRON_MAGNETIC_MOMENT));
  rb_define_const(module, "PROTON_MAGNETIC_MOMENT",
      rb_float_new(GSL_CONST_MKSA_PROTON_MAGNETIC_MOMENT));
  rb_define_const(module, "STANDARD_GAS_VOLUME",
      rb_float_new(GSL_CONST_MKSA_STANDARD_GAS_VOLUME));

  rb_define_const(module, "MINUTE", rb_float_new(GSL_CONST_MKSA_MINUTE));
  rb_define_const(module, "HOUR", rb_float_new(GSL_CONST_MKSA_HOUR));
  rb_define_const(module, "DAY", rb_float_new(GSL_CONST_MKSA_DAY));
  rb_define_const(module, "WEEK", rb_float_new(GSL_CONST_MKSA_WEEK));
  rb_define_const(module, "INCH", rb_float_new(GSL_CONST_MKSA_INCH));
  rb_define_const(module, "FOOT", rb_float_new(GSL_CONST_MKSA_FOOT));
  rb_define_const(module, "YARD", rb_float_new(GSL_CONST_MKSA_YARD));
  rb_define_const(module, "MILE", rb_float_new(GSL_CONST_MKSA_MILE));
  rb_define_const(module, "NAUTICAL_MILE",
      rb_float_new(GSL_CONST_MKSA_NAUTICAL_MILE));
  rb_define_const(module, "FATHOM", rb_float_new(GSL_CONST_MKSA_FATHOM));
  rb_define_const(module, "MIL", rb_float_new(GSL_CONST_MKSA_MIL));
  rb_define_const(module, "POINT", rb_float_new(GSL_CONST_MKSA_POINT));
  rb_define_const(module, "TEXPOINT", rb_float_new(GSL_CONST_MKSA_TEXPOINT));
  rb_define_const(module, "MICRON", rb_float_new(GSL_CONST_MKSA_MICRON));
  rb_define_const(module, "ANGSTROM", rb_float_new(GSL_CONST_MKSA_ANGSTROM));
  rb_define_const(module, "HECTARE", rb_float_new(GSL_CONST_MKSA_HECTARE));
  rb_define_const(module, "ACRE", rb_float_new(GSL_CONST_MKSA_ACRE));
  rb_define_const(module, "BARN", rb_float_new(GSL_CONST_MKSA_BARN));
  rb_define_const(module, "LITER", rb_float_new(GSL_CONST_MKSA_LITER));
  rb_define_const(module, "US_GALLON", rb_float_new(GSL_CONST_MKSA_US_GALLON));
  rb_define_const(module, "QUART", rb_float_new(GSL_CONST_MKSA_QUART));
  rb_define_const(module, "PINT", rb_float_new(GSL_CONST_MKSA_PINT));
  rb_define_const(module, "CUP", rb_float_new(GSL_CONST_MKSA_CUP));
  rb_define_const(module, "FLUID_OUNCE", rb_float_new(GSL_CONST_MKSA_FLUID_OUNCE));
  rb_define_const(module, "TABLESPOON", rb_float_new(GSL_CONST_MKSA_TABLESPOON));
  rb_define_const(module, "CANADIAN_GALLON",
      rb_float_new(GSL_CONST_MKSA_CANADIAN_GALLON));

  rb_define_const(module, "UK_GALLON", rb_float_new(GSL_CONST_MKSA_UK_GALLON));
  rb_define_const(module, "KILOMETERS_PER_HOUR",
      rb_float_new(GSL_CONST_MKSA_MILES_PER_HOUR));
  rb_define_const(module, "MILES_PER_HOUR",
      rb_float_new(GSL_CONST_MKSA_KILOMETERS_PER_HOUR));
  rb_define_const(module, "KNOT", rb_float_new(GSL_CONST_MKSA_KNOT));
  rb_define_const(module, "POUND_MASS", rb_float_new(GSL_CONST_MKSA_POUND_MASS));
  rb_define_const(module, "POUND_OUNCE", rb_float_new(GSL_CONST_MKSA_OUNCE_MASS));
  rb_define_const(module, "POUND_TON", rb_float_new(GSL_CONST_MKSA_TON));
  rb_define_const(module, "POUND_METRIC_TON",
      rb_float_new(GSL_CONST_MKSA_METRIC_TON));
  rb_define_const(module, "POUND_UK_TON", rb_float_new(GSL_CONST_MKSA_UK_TON));
  rb_define_const(module, "POUND_TROY_OUNCE",
      rb_float_new(GSL_CONST_MKSA_TROY_OUNCE));
  rb_define_const(module, "CARAT", rb_float_new(GSL_CONST_MKSA_CARAT));
  rb_define_const(module, "UNIFIED_ATOMIC_MASS",
      rb_float_new(GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS));
  rb_define_const(module, "GRAM_FORCE", rb_float_new(GSL_CONST_MKSA_GRAM_FORCE));
  rb_define_const(module, "POUND_FORCE", rb_float_new(GSL_CONST_MKSA_POUND_FORCE));
  rb_define_const(module, "KILOPOUND_FORCE",
      rb_float_new(GSL_CONST_MKSA_KILOPOUND_FORCE));
  rb_define_const(module, "POUNDAL", rb_float_new(GSL_CONST_MKSA_POUNDAL));
  rb_define_const(module, "CALORIE", rb_float_new(GSL_CONST_MKSA_CALORIE));
  rb_define_const(module, "BTU", rb_float_new(GSL_CONST_MKSA_BTU));
  rb_define_const(module, "THERM", rb_float_new(GSL_CONST_MKSA_THERM));
  rb_define_const(module, "HORSEPOWER", rb_float_new(GSL_CONST_MKSA_HORSEPOWER));
  rb_define_const(module, "BAR", rb_float_new(GSL_CONST_MKSA_BAR));
  rb_define_const(module, "STD_ATMOSPHERE",
      rb_float_new(GSL_CONST_MKSA_STD_ATMOSPHERE));
  rb_define_const(module, "TORR", rb_float_new(GSL_CONST_MKSA_TORR));
  rb_define_const(module, "METER_OF_MERCURY",
      rb_float_new(GSL_CONST_MKSA_METER_OF_MERCURY));
  rb_define_const(module, "INCH_OF_MERCURY",
      rb_float_new(GSL_CONST_MKSA_INCH_OF_MERCURY));
  rb_define_const(module, "INCH_OF_WATER",
      rb_float_new(GSL_CONST_MKSA_INCH_OF_WATER));
  rb_define_const(module, "PSI", rb_float_new(GSL_CONST_MKSA_PSI));
  rb_define_const(module, "POISE", rb_float_new(GSL_CONST_MKSA_POISE));
  rb_define_const(module, "STOKES", rb_float_new(GSL_CONST_MKSA_STOKES));
  rb_define_const(module, "FARADAY", rb_float_new(GSL_CONST_MKSA_FARADAY));
  rb_define_const(module, "ELECTRON_CHARGE",
      rb_float_new(GSL_CONST_MKSA_ELECTRON_CHARGE));
  rb_define_const(module, "GAUSS", rb_float_new(GSL_CONST_MKSA_GAUSS));
  rb_define_const(module, "STILB", rb_float_new(GSL_CONST_MKSA_STILB));
  rb_define_const(module, "LUMEN", rb_float_new(GSL_CONST_MKSA_LUMEN));
  rb_define_const(module, "LUX", rb_float_new(GSL_CONST_MKSA_LUX));
  rb_define_const(module, "PHOT", rb_float_new(GSL_CONST_MKSA_PHOT));
  rb_define_const(module, "FOOTCANDLE", rb_float_new(GSL_CONST_MKSA_FOOTCANDLE));
  rb_define_const(module, "LAMBERT", rb_float_new(GSL_CONST_MKSA_LAMBERT));
  rb_define_const(module, "CURIE", rb_float_new(GSL_CONST_MKSA_CURIE));
  rb_define_const(module, "ROENTGEN", rb_float_new(GSL_CONST_MKSA_ROENTGEN));
  rb_define_const(module, "RAD", rb_float_new(GSL_CONST_MKSA_RAD));
  rb_define_const(module, "SOLAR_MASS", rb_float_new(GSL_CONST_MKSA_SOLAR_MASS));

  rb_define_const(module, "BOLTZMANN", rb_float_new(GSL_CONST_MKSA_BOLTZMANN));
  rb_define_const(module, "MOLAR_GAS", rb_float_new(GSL_CONST_MKSA_MOLAR_GAS));

  rb_define_const(module, "BOHR_RADIUS", rb_float_new(GSL_CONST_MKSA_BOHR_RADIUS));
  rb_define_const(module, "NEWTON", rb_float_new(GSL_CONST_MKSA_NEWTON));
  rb_define_const(module, "DYNE", rb_float_new(GSL_CONST_MKSA_DYNE));
  rb_define_const(module, "JOULE", rb_float_new(GSL_CONST_MKSA_JOULE));
  rb_define_const(module, "ERG", rb_float_new(GSL_CONST_MKSA_ERG));

#ifdef GSL_1_4_9_LATER
  rb_define_const(module, "STEFAN_BOLTZMANN_CONSTANT",
      rb_float_new(GSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT));
  rb_define_const(module, "THOMSON_CROSS_SECTION",
      rb_float_new(GSL_CONST_MKSA_THOMSON_CROSS_SECTION));
#endif


#ifdef GSL_1_8_LATER
  rb_define_const(module, "DEBYE", rb_float_new(GSL_CONST_MKSA_DEBYE));
#endif

}


static void rb_gsl_const_cgs(VALUE module)
{
  rb_define_const(module, "SPEED_OF_LIGHT",
      rb_float_new(GSL_CONST_CGSM_SPEED_OF_LIGHT));
  rb_define_const(module, "GRAVITATIONAL_CONSTANT",
      rb_float_new(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT));
  rb_define_const(module, "PLANCKS_CONSTANT_H",
      rb_float_new(GSL_CONST_CGSM_PLANCKS_CONSTANT_H));
  rb_define_const(module, "PLANCKS_CONSTANT_HBAR",
      rb_float_new(GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR));
  rb_define_const(module, "ASTRONOMICAL_UNIT",
      rb_float_new(GSL_CONST_CGSM_ASTRONOMICAL_UNIT));
  rb_define_const(module, "LIGHT_YEAR",
      rb_float_new(GSL_CONST_CGSM_LIGHT_YEAR));
  rb_define_const(module, "PARSEC",
      rb_float_new(GSL_CONST_CGSM_PARSEC));
  rb_define_const(module, "GRAV_ACCEL",
      rb_float_new(GSL_CONST_CGSM_GRAV_ACCEL));
  rb_define_const(module, "ELECTRON_VOLT",
      rb_float_new(GSL_CONST_CGSM_ELECTRON_VOLT));
  rb_define_const(module, "MASS_ELECTRON",
      rb_float_new(GSL_CONST_CGSM_MASS_ELECTRON));
  rb_define_const(module, "MASS_MUON", rb_float_new(GSL_CONST_CGSM_MASS_MUON));
  rb_define_const(module, "MASS_PROTON", rb_float_new(GSL_CONST_CGSM_MASS_PROTON));
  rb_define_const(module, "MASS_NEUTRON", rb_float_new(GSL_CONST_CGSM_MASS_NEUTRON));
  rb_define_const(module, "RYDBERG", rb_float_new(GSL_CONST_CGSM_RYDBERG));

  rb_define_const(module, "BOHR_MAGNETON",
      rb_float_new(GSL_CONST_CGSM_BOHR_MAGNETON));

  rb_define_const(module, "NUCLEAR_MAGNETON",
      rb_float_new(GSL_CONST_CGSM_NUCLEAR_MAGNETON));
  rb_define_const(module, "ELECTRON_MAGNETIC_MOMENT",
      rb_float_new(GSL_CONST_CGSM_ELECTRON_MAGNETIC_MOMENT));
  rb_define_const(module, "PROTON_MAGNETIC_MOMENT",
      rb_float_new(GSL_CONST_CGSM_PROTON_MAGNETIC_MOMENT));
  rb_define_const(module, "STANDARD_GAS_VOLUME",
      rb_float_new(GSL_CONST_CGSM_STANDARD_GAS_VOLUME));

  rb_define_const(module, "MINUTE", rb_float_new(GSL_CONST_CGSM_MINUTE));
  rb_define_const(module, "HOUR", rb_float_new(GSL_CONST_CGSM_HOUR));
  rb_define_const(module, "DAY", rb_float_new(GSL_CONST_CGSM_DAY));
  rb_define_const(module, "WEEK", rb_float_new(GSL_CONST_CGSM_WEEK));
  rb_define_const(module, "INCH", rb_float_new(GSL_CONST_CGSM_INCH));
  rb_define_const(module, "FOOT", rb_float_new(GSL_CONST_CGSM_FOOT));
  rb_define_const(module, "YARD", rb_float_new(GSL_CONST_CGSM_YARD));
  rb_define_const(module, "MILE", rb_float_new(GSL_CONST_CGSM_MILE));
  rb_define_const(module, "NAUTICAL_MILE",
      rb_float_new(GSL_CONST_CGSM_NAUTICAL_MILE));
  rb_define_const(module, "FATHOM", rb_float_new(GSL_CONST_CGSM_FATHOM));
  rb_define_const(module, "MIL", rb_float_new(GSL_CONST_CGSM_MIL));
  rb_define_const(module, "POINT", rb_float_new(GSL_CONST_CGSM_POINT));
  rb_define_const(module, "TEXPOINT", rb_float_new(GSL_CONST_CGSM_TEXPOINT));
  rb_define_const(module, "MICRON", rb_float_new(GSL_CONST_CGSM_MICRON));
  rb_define_const(module, "ANGSTROM", rb_float_new(GSL_CONST_CGSM_ANGSTROM));
  rb_define_const(module, "HECTARE", rb_float_new(GSL_CONST_CGSM_HECTARE));
  rb_define_const(module, "ACRE", rb_float_new(GSL_CONST_CGSM_ACRE));
  rb_define_const(module, "BARN", rb_float_new(GSL_CONST_CGSM_BARN));
  rb_define_const(module, "LITER", rb_float_new(GSL_CONST_CGSM_LITER));
  rb_define_const(module, "US_GALLON", rb_float_new(GSL_CONST_CGSM_US_GALLON));
  rb_define_const(module, "QUART", rb_float_new(GSL_CONST_CGSM_QUART));
  rb_define_const(module, "PINT", rb_float_new(GSL_CONST_CGSM_PINT));
  rb_define_const(module, "CUP", rb_float_new(GSL_CONST_CGSM_CUP));
  rb_define_const(module, "FLUID_OUNCE", rb_float_new(GSL_CONST_CGSM_FLUID_OUNCE));
  rb_define_const(module, "TABLESPOON", rb_float_new(GSL_CONST_CGSM_TABLESPOON));
  rb_define_const(module, "CANADIAN_GALLON",
      rb_float_new(GSL_CONST_CGSM_CANADIAN_GALLON));

  rb_define_const(module, "UK_GALLON", rb_float_new(GSL_CONST_CGSM_UK_GALLON));
  rb_define_const(module, "KILOMETERS_PER_HOUR",
      rb_float_new(GSL_CONST_CGSM_MILES_PER_HOUR));
  rb_define_const(module, "MILES_PER_HOUR",
      rb_float_new(GSL_CONST_CGSM_KILOMETERS_PER_HOUR));
  rb_define_const(module, "KNOT", rb_float_new(GSL_CONST_CGSM_KNOT));
  rb_define_const(module, "POUND_MASS", rb_float_new(GSL_CONST_CGSM_POUND_MASS));
  rb_define_const(module, "POUND_OUNCE", rb_float_new(GSL_CONST_CGSM_OUNCE_MASS));
  rb_define_const(module, "POUND_TON", rb_float_new(GSL_CONST_CGSM_TON));
  rb_define_const(module, "POUND_METRIC_TON",
      rb_float_new(GSL_CONST_CGSM_METRIC_TON));
  rb_define_const(module, "POUND_UK_TON", rb_float_new(GSL_CONST_CGSM_UK_TON));
  rb_define_const(module, "POUND_TROY_OUNCE",
      rb_float_new(GSL_CONST_CGSM_TROY_OUNCE));
  rb_define_const(module, "CARAT", rb_float_new(GSL_CONST_CGSM_CARAT));
  rb_define_const(module, "UNIFIED_ATOMIC_MASS",
      rb_float_new(GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS));
  rb_define_const(module, "GRAM_FORCE", rb_float_new(GSL_CONST_CGSM_GRAM_FORCE));
  rb_define_const(module, "POUND_FORCE", rb_float_new(GSL_CONST_CGSM_POUND_FORCE));
  rb_define_const(module, "KILOPOUND_FORCE",
      rb_float_new(GSL_CONST_CGSM_KILOPOUND_FORCE));
  rb_define_const(module, "POUNDAL", rb_float_new(GSL_CONST_CGSM_POUNDAL));
  rb_define_const(module, "CALORIE", rb_float_new(GSL_CONST_CGSM_CALORIE));
  rb_define_const(module, "BTU", rb_float_new(GSL_CONST_CGSM_BTU));
  rb_define_const(module, "THERM", rb_float_new(GSL_CONST_CGSM_THERM));
  rb_define_const(module, "HORSEPOWER", rb_float_new(GSL_CONST_CGSM_HORSEPOWER));
  rb_define_const(module, "BAR", rb_float_new(GSL_CONST_CGSM_BAR));
  rb_define_const(module, "STD_ATMOSPHERE",
      rb_float_new(GSL_CONST_CGSM_STD_ATMOSPHERE));
  rb_define_const(module, "TORR", rb_float_new(GSL_CONST_CGSM_TORR));
  rb_define_const(module, "METER_OF_MERCURY",
      rb_float_new(GSL_CONST_CGSM_METER_OF_MERCURY));
  rb_define_const(module, "INCH_OF_MERCURY",
      rb_float_new(GSL_CONST_CGSM_INCH_OF_MERCURY));
  rb_define_const(module, "INCH_OF_WATER",
      rb_float_new(GSL_CONST_CGSM_INCH_OF_WATER));
  rb_define_const(module, "PSI", rb_float_new(GSL_CONST_CGSM_PSI));
  rb_define_const(module, "POISE", rb_float_new(GSL_CONST_CGSM_POISE));
  rb_define_const(module, "STOKES", rb_float_new(GSL_CONST_CGSM_STOKES));
  rb_define_const(module, "FARADAY", rb_float_new(GSL_CONST_CGSM_FARADAY));
  rb_define_const(module, "ELECTRON_CHARGE",
      rb_float_new(GSL_CONST_CGSM_ELECTRON_CHARGE));
  rb_define_const(module, "ELECTRON_CHARGE_ESU",
      rb_float_new(GSL_CONST_CGSM_ELECTRON_CHARGE*GSL_CONST_CGSM_SPEED_OF_LIGHT));
#ifndef GSL_1_13_LATER
  rb_define_const(module, "GAUSS", rb_float_new(GSL_CONST_CGSM_GAUSS));
#endif
  rb_define_const(module, "STILB", rb_float_new(GSL_CONST_CGSM_STILB));
  rb_define_const(module, "LUMEN", rb_float_new(GSL_CONST_CGSM_LUMEN));
  rb_define_const(module, "LUX", rb_float_new(GSL_CONST_CGSM_LUX));
  rb_define_const(module, "PHOT", rb_float_new(GSL_CONST_CGSM_PHOT));
  rb_define_const(module, "FOOTCANDLE", rb_float_new(GSL_CONST_CGSM_FOOTCANDLE));
  rb_define_const(module, "LAMBERT", rb_float_new(GSL_CONST_CGSM_LAMBERT));
  rb_define_const(module, "CURIE", rb_float_new(GSL_CONST_CGSM_CURIE));
  rb_define_const(module, "ROENTGEN", rb_float_new(GSL_CONST_CGSM_ROENTGEN));
  rb_define_const(module, "RAD", rb_float_new(GSL_CONST_CGSM_RAD));
  rb_define_const(module, "SOLAR_MASS", rb_float_new(GSL_CONST_CGSM_SOLAR_MASS));

  rb_define_const(module, "BOLTZMANN", rb_float_new(GSL_CONST_CGSM_BOLTZMANN));
  rb_define_const(module, "MOLAR_GAS", rb_float_new(GSL_CONST_CGSM_MOLAR_GAS));

  rb_define_const(module, "BOHR_RADIUS", rb_float_new(GSL_CONST_CGSM_BOHR_RADIUS));
  rb_define_const(module, "NEWTON", rb_float_new(GSL_CONST_CGSM_NEWTON));
  rb_define_const(module, "DYNE", rb_float_new(GSL_CONST_CGSM_DYNE));
  rb_define_const(module, "JOULE", rb_float_new(GSL_CONST_CGSM_JOULE));
  rb_define_const(module, "ERG", rb_float_new(GSL_CONST_CGSM_ERG));

#ifdef GSL_1_4_9_LATER
  rb_define_const(module, "STEFAN_BOLTZMANN_CONSTANT",
      rb_float_new(GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT));
  rb_define_const(module, "THOMSON_CROSS_SECTION",
      rb_float_new(GSL_CONST_CGSM_THOMSON_CROSS_SECTION));
#endif
}

static void rb_gsl_const_num(VALUE module)
{
  rb_define_const(module, "AVOGADRO", rb_float_new(GSL_CONST_NUM_AVOGADRO));
  rb_define_const(module, "FINE_STRUCTURE",
      rb_float_new(GSL_CONST_NUM_FINE_STRUCTURE));
  rb_define_const(module, "YOTTA", rb_float_new(GSL_CONST_NUM_YOTTA));
  rb_define_const(module, "ZETTA", rb_float_new(GSL_CONST_NUM_ZETTA));
  rb_define_const(module, "EXA", rb_float_new(GSL_CONST_NUM_EXA));
  rb_define_const(module, "PETA", rb_float_new(GSL_CONST_NUM_PETA));
  rb_define_const(module, "TERA", rb_float_new(GSL_CONST_NUM_TERA));
  rb_define_const(module, "GIGA", rb_float_new(GSL_CONST_NUM_GIGA));
  rb_define_const(module, "MEGA", rb_float_new(GSL_CONST_NUM_MEGA));
  rb_define_const(module, "KILO", rb_float_new(GSL_CONST_NUM_KILO));
  rb_define_const(module, "MILLI", rb_float_new(GSL_CONST_NUM_MILLI));
  rb_define_const(module, "MICRO", rb_float_new(GSL_CONST_NUM_MICRO));
  rb_define_const(module, "NANO", rb_float_new(GSL_CONST_NUM_NANO));
  rb_define_const(module, "PICO", rb_float_new(GSL_CONST_NUM_PICO));
  rb_define_const(module, "FEMTO", rb_float_new(GSL_CONST_NUM_FEMTO));
  rb_define_const(module, "ATTO", rb_float_new(GSL_CONST_NUM_ATTO));
  rb_define_const(module, "ZEPTO", rb_float_new(GSL_CONST_NUM_ZEPTO));
  rb_define_const(module, "YOCTO", rb_float_new(GSL_CONST_NUM_YOCTO));
}

void Init_gsl_const(VALUE module)
{
  VALUE mgsl_const;
  VALUE mgsl_const_mks, mgsl_const_cgs, mgsl_const_num;

  mgsl_const = rb_define_module_under(module, "CONST");
  mgsl_const_mks = rb_define_module_under(mgsl_const, "MKSA");
  rb_gsl_const_mks(mgsl_const_mks);
  mgsl_const_cgs = rb_define_module_under(mgsl_const, "CGSM");
  rb_gsl_const_cgs(mgsl_const_cgs);
  mgsl_const_num = rb_define_module_under(mgsl_const, "NUM");
  rb_gsl_const_num(mgsl_const_num);
  Init_gsl_const_additional(mgsl_const_mks, mgsl_const_cgs, mgsl_const_num);
}
#endif
