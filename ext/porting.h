/*
  porting.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or POLYNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef ___RB_GSL_PORTING_H___
#define ___RB_GSL_PORTING_H___

#include "ruby.h"

#if RUBY_API_VERSION_MAJOR > 2 || (RUBY_API_VERSION_MAJOR == 2 && RUBY_API_VERSION_MINOR >= 1)
// Ruby < v2.1

#define rb_obj_reveal(instance, klass_value) RBASIC(instance)->klass = klass_value

#endif  // Ruby < v2.1

#endif
