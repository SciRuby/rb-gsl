/////////////////////////////////////////////////////////////////////
// = NMatrix
//
// A linear algebra library for scientific computation in Ruby.
// NMatrix is part of SciRuby.
//
// NMatrix was originally inspired by and derived from NArray, by
// Masahiro Tanaka: http://narray.rubyforge.org
//
// == Copyright Information
//
// SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
// NMatrix is Copyright (c) 2012, Ruby Science Foundation
//
// Please see LICENSE.txt for additional copyright notices.
//
// == Contributing
//
// By contributing source code to SciRuby, you agree to be bound by
// our Contributor Agreement:
//
// * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
//
// == rational.c
//
// This file is largely based off of Ruby 1.9.3's rational.c. It
// contains functions for dealing with rational types which are
// not Ruby VALUE-based.

#ifndef RATIONAL_C
# define RATIONAL_C

#include "nmatrix.h"


inline int64_t nmrb_gcd(int64_t x, int64_t y) {
  int64_t t;

  if (x < 0) x = -x;
  if (y < 0) y = -y;

  if (x == 0) return y;
  if (y == 0) return x;

  while (x > 0) {
    t = x;
    x = y % x;
    y = t;
  }

  return y;
}

/*
static double f_gcd(double x, double y) {
  double z;

  if (x < 0.0) x = -x;
  if (y < 0.0) y = -y;
  if (x == 0.0) return y;
  if (y == 0.0) return x;

  for (;;) {
    z = x;
    x = y % x;
    y = z;
  }
  // NOTREACHED
}*/

/*
inline VALUE nmrb_rational_new_no_reduce1(VALUE klass, int64_t x) {
  return nurat_s_canonicalize_internal_no_reduce(klass, x, 1);
}

inline static VALUE nurat_s_canonicalize_internal_no_reduce(VALUE klass, int64_t num, int64_t den) {
  if (den < 0) {
    num = -num;
    den = -den;
  } else if (den == 0) {
    rb_raise_zerodiv();
  }
  return nurat_s_new_internal(klass, num, den);
}*/
/*
inline VALUE nmrb_rational_new(VALUE klass, int64_t num, int64_t den) {
  NEWOBJ(obj, struct RRational);
  OBJSETUP(obj, klass, T_RATIONAL);

  obj->num = INT2FIX(num);
  obj->den = INT2FIX(den);

  return (VALUE)obj;
}
*/


#endif
