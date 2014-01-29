/*
  rb_gsl_roots.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY
*/

#ifndef ___RB_GSL_ROOT_H___
#define ___RB_GSL_ROOT_H___

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "rb_gsl.h"

EXTERN VALUE cgsl_fsolver;
EXTERN VALUE cgsl_fdfsolver;

#endif
