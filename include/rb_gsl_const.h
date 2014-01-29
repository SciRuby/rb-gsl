/*
  rb_gsl.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef ___RB_GSL_CONST_H___
#define ___RB_GSL_CONST_H___

#include "rb_gsl.h"
#ifdef GSL_CONST_OLD
#include <gsl/gsl_const_mks.h>
#include <gsl/gsl_const_cgs.h>
#else
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_cgsm.h>
#endif
#include <gsl/gsl_const_num.h>

EXTERN VALUE mgsl_const_mks, mgsl_const_cgs;

#endif
