/*
  tensor.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

/*
  The gsl_tensor package is developed by J.Burguet, and
  distributed separately as an add-on package.
 */

#ifdef HAVE_TENSOR_TENSOR_H

#include "include/rb_gsl_tensor.h"

#ifdef HAVE_NARRAY_H
#include "include/rb_gsl_with_narray.h"
#endif

#define BASE_DOUBLE
#include "include/templates_on.h"
#include "tensor_source.c"
#include "include/templates_off.h"
#undef  BASE_DOUBLE

#define BASE_INT
#include "include/templates_on.h"
#include "tensor_source.c"
#include "include/templates_off.h"
#undef  BASE_INT

#endif
