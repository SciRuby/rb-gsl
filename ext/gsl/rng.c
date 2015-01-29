/*
  rng.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

/*
   Document-class: <i>GSL::Rng</i>
   Random number generator
*/

#include "include/rb_gsl_rng.h"
#ifdef HAVE_RNGEXTRA_RNGEXTRA_H
#include "rngextra/rngextra.h"
#endif

/* Global, since used in other source files */
VALUE cgsl_rng;

enum rb_gsl_rng_generator {
  GSL_RNG_DEFAULT,
  GSL_RNG_MT19937, GSL_RNG_MT19937_1999, GSL_RNG_MT19937_1998,
  GSL_RNG_RANLXS0, GSL_RNG_RANLXS1, GSL_RNG_RANLXS2,
  GSL_RNG_RANLXD1, GSL_RNG_RANLXD2,
  GSL_RNG_RANLUX, GSL_RNG_RANLUX389,
  GSL_RNG_CMRG, GSL_RNG_MRG,
  GSL_RNG_TAUS, GSL_RNG_TAUS2, GSL_RNG_TAUS113,  GSL_RNG_GFSR4,
  GSL_RNG_RAND,
  GSL_RNG_RANDOM_BSD, GSL_RNG_RANDOM_GLIBC2,
  GSL_RNG_RANDOM8_GLIBC2, GSL_RNG_RANDOM32_GLIBC2, GSL_RNG_RANDOM64_GLIBC2,
  GSL_RNG_RANDOM128_GLIBC2, GSL_RNG_RANDOM256_GLIBC2,
  GSL_RNG_RANDOM8_BSD, GSL_RNG_RANDOM32_BSD, GSL_RNG_RANDOM64_BSD,
  GSL_RNG_RANDOM128_BSD, GSL_RNG_RANDOM256_BSD,
  GSL_RNG_RANDOM_LIBC5, GSL_RNG_RANDOM8_LIBC5, GSL_RNG_RANDOM32_LIBC5,
  GSL_RNG_RANDOM64_LIBC5, GSL_RNG_RANDOM128_LIBC5, GSL_RNG_RANDOM256_LIBC5,
  GSL_RNG_RAND48,
  GSL_RNG_RAN0, GSL_RNG_RAN1, GSL_RNG_RAN2, GSL_RNG_RAN3,
  GSL_RNG_RANF, GSL_RNG_RANMAR, GSL_RNG_R250, GSL_RNG_TT800, GSL_RNG_VAX,
  GSL_RNG_TRANSPUTER, GSL_RNG_RANDU, GSL_RNG_MINSTD,
  GSL_RNG_UNI, GSL_RNG_UNI32, GSL_RNG_SLATEC, GSL_RNG_ZUF,
  /* from gsl-1.1 */
  GSL_RNG_BOROSH13, GSL_RNG_COVEYOU, GSL_RNG_FISHMAN18, GSL_RNG_FISHMAN20,
  GSL_RNG_FISHMAN2X, GSL_RNG_KNUTHRAN, GSL_RNG_KNUTHRAN2,
  GSL_RNG_LECUYER21, GSL_RNG_WATERMAN14,
  /* For the rngextra package */
  /* Tue Oct 19 22:55:56 JST 2004 by Y.Tsunesada */
  GSL_RNGEXTRA_RNG1, GSL_RNGEXTRA_RNG2,
  /* GSL-1.9 */
  GSL_RNG_KNUTHRAN2002,
};

static const gsl_rng_type* get_gsl_rng_type(VALUE t);
/*
  Document-method: <i>GSL::Rng.alloc</i>
    Constructor.
*/
static VALUE rb_gsl_rng_alloc(int argc, VALUE *argv, VALUE klass)
{
  gsl_rng *r = NULL;
  const gsl_rng_type *T;
  unsigned long seed;
  int itype;
  gsl_rng_env_setup();
  if (argc == 0) {
    T = gsl_rng_default;
    seed = gsl_rng_default_seed;
  } else {
    T = get_gsl_rng_type(argv[0]);
    if (argc == 1) {
      seed = gsl_rng_default_seed;
    } else if (argc == 2) {
      itype = TYPE(argv[1]);
      if (itype == T_FIXNUM || itype == T_BIGNUM) {
  seed = FIX2INT(argv[1]);
      } else {
  rb_raise(rb_eArgError,
     "bad argument 2, seed must be an integer.");
      }
    } else {
      rb_raise(rb_eArgError, "too many arguments (%d for 0 or 1)", argc);
    }
  }
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);
  return Data_Wrap_Struct(klass, 0, gsl_rng_free, r);
}

static const gsl_rng_type* get_gsl_rng_type_int(int itype);
static const gsl_rng_type* get_gsl_rng_type_name(char *name);
static const gsl_rng_type* get_gsl_rng_type(VALUE t)
{
  switch (TYPE(t)) {
  case T_STRING:
    return get_gsl_rng_type_name(STR2CSTR(t));
    break;
  case T_FIXNUM:
    return get_gsl_rng_type_int(FIX2INT(t));
    break;
  default:
    rb_raise(rb_eTypeError, "String or Fixnum expected");
  }
}

static const gsl_rng_type* get_gsl_rng_type_name(char *name)
{
  if (str_tail_grep(name, "default") == 0) return gsl_rng_default;
  else if (str_tail_grep(name, "mt19937") == 0) return gsl_rng_mt19937;
#ifdef GSL_1_1_LATER
  else if (str_tail_grep(name, "borosh13") == 0) return gsl_rng_borosh13;
  else if (str_tail_grep(name, "coveyou") == 0) return gsl_rng_coveyou;
  else if (str_tail_grep(name, "fishman18") == 0) return gsl_rng_fishman18;
  else if (str_tail_grep(name, "fishman20") == 0) return gsl_rng_fishman20;
  else if (str_tail_grep(name, "fishman2x") == 0) return gsl_rng_fishman2x;
  else if (str_tail_grep(name, "lecuyer21") == 0) return gsl_rng_lecuyer21;
  else if (str_tail_grep(name, "waterman14") == 0) return gsl_rng_waterman14;
  else if (str_tail_grep(name, "knuthran") == 0) return gsl_rng_knuthran;
  else if (str_tail_grep(name, "knuthran2") == 0) return gsl_rng_knuthran2;
#endif
#ifdef GSL_1_2_LATER
  else if (str_tail_grep(name, "mt19937_1999") == 0) return gsl_rng_mt19937_1999;
  else if (str_tail_grep(name, "mt19937-1999") == 0) return gsl_rng_mt19937_1999;
  else if (str_tail_grep(name, "mt19937_1998") == 0) return gsl_rng_mt19937_1998;
  else if (str_tail_grep(name, "mt19937-1998") == 0) return gsl_rng_mt19937_1998;
  else if (str_tail_grep(name, "taus113") == 0) return gsl_rng_taus113;
  else if (str_tail_grep(name, "taus2") == 0) return gsl_rng_taus2;
#endif
  else if (str_tail_grep(name, "mt19937") == 0) return gsl_rng_mt19937;
  else if (str_tail_grep(name, "ranlxs0") == 0) return gsl_rng_ranlxs0;
  else if (str_tail_grep(name, "ranlxs1") == 0) return gsl_rng_ranlxs1;
  else if (str_tail_grep(name, "ranlxs2") == 0) return gsl_rng_ranlxs2;
  else if (str_tail_grep(name, "ranlxd1") == 0) return gsl_rng_ranlxd1;
  else if (str_tail_grep(name, "ranlxd2") == 0) return gsl_rng_ranlxd2;
  else if (str_tail_grep(name, "ranlux") == 0) return gsl_rng_ranlux;
  else if (str_tail_grep(name, "ranlux389") == 0) return gsl_rng_ranlux389;
  else if (str_tail_grep(name, "cmrg") == 0) return gsl_rng_cmrg;
  else if (str_tail_grep(name, "mrg") == 0) return gsl_rng_mrg;
  else if (str_tail_grep(name, "taus") == 0) return gsl_rng_taus;
  else if (str_tail_grep(name, "gfsr4") == 0) return gsl_rng_gfsr4;
  else if (str_tail_grep(name, "rand") == 0) return gsl_rng_rand;
  else if (str_tail_grep(name, "random_libc5") == 0) return gsl_rng_random_libc5;
  else if (str_tail_grep(name, "random8_libc5") == 0) return gsl_rng_random8_libc5;
  else if (str_tail_grep(name, "random32_libc5") == 0) return gsl_rng_random32_libc5;
  else if (str_tail_grep(name, "random64_libc5") == 0) return gsl_rng_random64_libc5;
  else if (str_tail_grep(name, "random128_libc5") == 0) return gsl_rng_random128_libc5;
  else if (str_tail_grep(name, "random256_libc5") == 0) return gsl_rng_random256_libc5;
  else if (str_tail_grep(name, "random-libc5") == 0) return gsl_rng_random_libc5;
  else if (str_tail_grep(name, "random8-libc5") == 0) return gsl_rng_random8_libc5;
  else if (str_tail_grep(name, "random32-libc5") == 0) return gsl_rng_random32_libc5;
  else if (str_tail_grep(name, "random64-libc5") == 0) return gsl_rng_random64_libc5;
  else if (str_tail_grep(name, "random128-libc5") == 0) return gsl_rng_random128_libc5;
  else if (str_tail_grep(name, "random256-libc5") == 0) return gsl_rng_random256_libc5;
  else if (str_tail_grep(name, "random_glibc2") == 0) return gsl_rng_random_glibc2;
  else if (str_tail_grep(name, "random8_glibc2") == 0) return gsl_rng_random8_glibc2;
  else if (str_tail_grep(name, "random32_glibc2") == 0) return gsl_rng_random32_glibc2;
  else if (str_tail_grep(name, "random64_glibc2") == 0) return gsl_rng_random64_glibc2;
  else if (str_tail_grep(name, "random128_glibc2") == 0) return gsl_rng_random128_glibc2;
  else if (str_tail_grep(name, "random256_glibc2") == 0) return gsl_rng_random256_glibc2;
  else if (str_tail_grep(name, "random-glibc2") == 0) return gsl_rng_random_glibc2;
  else if (str_tail_grep(name, "random8-glibc2") == 0) return gsl_rng_random8_glibc2;
  else if (str_tail_grep(name, "random32-glibc2") == 0) return gsl_rng_random32_glibc2;
  else if (str_tail_grep(name, "random64-glibc2") == 0) return gsl_rng_random64_glibc2;
  else if (str_tail_grep(name, "random128-glibc2") == 0) return gsl_rng_random128_glibc2;
  else if (str_tail_grep(name, "random256-glibc2") == 0) return gsl_rng_random256_glibc2;
  else if (str_tail_grep(name, "random_bsd") == 0) return gsl_rng_random_bsd;
  else if (str_tail_grep(name, "random8_bsd") == 0) return gsl_rng_random8_bsd;
  else if (str_tail_grep(name, "random32_bsd") == 0) return gsl_rng_random32_bsd;
  else if (str_tail_grep(name, "random64_bsd") == 0) return gsl_rng_random64_bsd;
  else if (str_tail_grep(name, "random128_bsd") == 0) return gsl_rng_random128_bsd;
  else if (str_tail_grep(name, "random256_bsd") == 0) return gsl_rng_random256_bsd;
  else if (str_tail_grep(name, "random-bsd") == 0) return gsl_rng_random_bsd;
  else if (str_tail_grep(name, "random8-bsd") == 0) return gsl_rng_random8_bsd;
  else if (str_tail_grep(name, "random32-bsd") == 0) return gsl_rng_random32_bsd;
  else if (str_tail_grep(name, "random64-bsd") == 0) return gsl_rng_random64_bsd;
  else if (str_tail_grep(name, "random128-bsd") == 0) return gsl_rng_random128_bsd;
  else if (str_tail_grep(name, "random256-bsd") == 0) return gsl_rng_random256_bsd;
  else if (str_tail_grep(name, "rand48") == 0) return gsl_rng_rand48;
  else if (str_tail_grep(name, "ran0") == 0) return gsl_rng_ran0;
  else if (str_tail_grep(name, "ran1") == 0) return gsl_rng_ran1;
  else if (str_tail_grep(name, "ran2") == 0) return gsl_rng_ran2;
  else if (str_tail_grep(name, "ran3") == 0) return gsl_rng_ran3;
  else if (str_tail_grep(name, "ranf") == 0) return gsl_rng_ranf;
  else if (str_tail_grep(name, "ranmar") == 0) return gsl_rng_ranmar;
  else if (str_tail_grep(name, "r250") == 0) return gsl_rng_r250;
  else if (str_tail_grep(name, "tt800") == 0) return gsl_rng_tt800;
  else if (str_tail_grep(name, "vax") == 0) return gsl_rng_vax;
  else if (str_tail_grep(name, "transputer") == 0) return gsl_rng_transputer;
  else if (str_tail_grep(name, "randu") == 0) return gsl_rng_randu;
  else if (str_tail_grep(name, "minstd") == 0) return gsl_rng_minstd;
  else if (str_tail_grep(name, "uni") == 0) return gsl_rng_uni;
  else if (str_tail_grep(name, "uni32") == 0) return gsl_rng_uni32;
  else if (str_tail_grep(name, "slatec") == 0) return gsl_rng_slatec;
  else if (str_tail_grep(name, "zuf") == 0) return gsl_rng_zuf;
#ifdef HAVE_RNGEXTRA_RNGEXTRA_H
  else if (str_tail_grep(name, "rngextra_rng1") == 0) return rngextra_rng1;
  else if (str_tail_grep(name, "rngextra_rng2") == 0) return rngextra_rng2;
  else if (str_tail_grep(name, "rngextra-rng1") == 0) return rngextra_rng1;
  else if (str_tail_grep(name, "rngextra-rng2") == 0) return rngextra_rng2;
#else
  else if (str_tail_grep(name, "rngextra_rng1")*str_tail_grep(name, "rngextra_rng2") == 0)
    rb_raise(rb_eNotImpError, "Install the rngextra package found at <http://www.network-theory.co.uk/download/rngextra/>.");
  else if (str_tail_grep(name, "rngextra_rng2")*str_tail_grep(name, "rngextra_rng2") == 0)
    rb_raise(rb_eNotImpError, "Install the rngextra package found at <http://www.network-theory.co.uk/download/rngextra/>.");
  else if (str_tail_grep(name, "rngextra-rng1")*str_tail_grep(name, "rngextra_rng2") == 0)
    rb_raise(rb_eNotImpError, "Install the rngextra package found at <http://www.network-theory.co.uk/download/rngextra/>.");
  else if (str_tail_grep(name, "rngextra-rng2")*str_tail_grep(name, "rngextra_rng2") == 0)
    rb_raise(rb_eNotImpError, "Install the rngextra package found at <http://www.network-theory.co.uk/download/rngextra/>.");
#endif
#ifdef GSL_1_9_LATER
  else if (str_tail_grep(name, "knuthran2002") == 0) return gsl_rng_knuthran2002;
#endif
  else
    rb_raise(rb_eArgError, "unknown generator type \"%s\"", name);
}

static const gsl_rng_type* get_gsl_rng_type_int(int itype)
{
  const gsl_rng_type *T;

  switch (itype) {
  case GSL_RNG_DEFAULT: T = gsl_rng_default; break;
  case GSL_RNG_MT19937: T = gsl_rng_mt19937; break; /* default */
#ifdef GSL_1_2_LATER
  case GSL_RNG_MT19937_1999: T = gsl_rng_mt19937_1999; break;
  case GSL_RNG_MT19937_1998: T = gsl_rng_mt19937_1998; break;
  case GSL_RNG_TAUS113: T = gsl_rng_taus113; break;
  case GSL_RNG_TAUS2: T = gsl_rng_taus2; break;
#endif
  case GSL_RNG_RANLXS0: T = gsl_rng_ranlxs0; break;
  case GSL_RNG_RANLXS1: T = gsl_rng_ranlxs1; break;
  case GSL_RNG_RANLXS2: T = gsl_rng_ranlxs2; break;
  case GSL_RNG_RANLXD1: T = gsl_rng_ranlxd1; break;
  case GSL_RNG_RANLXD2: T = gsl_rng_ranlxd2; break;
  case GSL_RNG_RANLUX: T = gsl_rng_ranlux; break;
  case GSL_RNG_RANLUX389: T = gsl_rng_ranlux389; break;
  case GSL_RNG_CMRG: T = gsl_rng_cmrg; break;
  case GSL_RNG_MRG: T = gsl_rng_mrg; break;
  case GSL_RNG_TAUS: T = gsl_rng_taus; break;
  case GSL_RNG_GFSR4: T = gsl_rng_gfsr4; break;
  case GSL_RNG_RAND: T = gsl_rng_rand; break;
  case GSL_RNG_RANDOM_LIBC5: T = gsl_rng_random_libc5; break;
  case GSL_RNG_RANDOM8_LIBC5: T = gsl_rng_random8_libc5; break;
  case GSL_RNG_RANDOM32_LIBC5: T = gsl_rng_random32_libc5; break;
  case GSL_RNG_RANDOM64_LIBC5: T = gsl_rng_random64_libc5; break;
  case GSL_RNG_RANDOM128_LIBC5: T = gsl_rng_random128_libc5; break;
  case GSL_RNG_RANDOM256_LIBC5: T = gsl_rng_random256_libc5; break;
  case GSL_RNG_RANDOM_GLIBC2: T = gsl_rng_random_glibc2; break;
  case GSL_RNG_RANDOM8_GLIBC2: T = gsl_rng_random8_glibc2; break;
  case GSL_RNG_RANDOM32_GLIBC2: T = gsl_rng_random32_glibc2; break;
  case GSL_RNG_RANDOM64_GLIBC2: T = gsl_rng_random64_glibc2; break;
  case GSL_RNG_RANDOM128_GLIBC2: T = gsl_rng_random128_glibc2; break;
  case GSL_RNG_RANDOM256_GLIBC2: T = gsl_rng_random256_glibc2; break;
  case GSL_RNG_RANDOM_BSD: T = gsl_rng_random_bsd; break;
  case GSL_RNG_RANDOM8_BSD: T = gsl_rng_random8_bsd; break;
  case GSL_RNG_RANDOM32_BSD: T = gsl_rng_random32_bsd; break;
  case GSL_RNG_RANDOM64_BSD: T = gsl_rng_random64_bsd; break;
  case GSL_RNG_RANDOM128_BSD: T = gsl_rng_random128_bsd; break;
  case GSL_RNG_RANDOM256_BSD: T = gsl_rng_random256_bsd; break;
  case GSL_RNG_RAND48: T = gsl_rng_rand48; break;
  case GSL_RNG_RAN0: T = gsl_rng_ran0; break;
  case GSL_RNG_RAN1: T = gsl_rng_ran1; break;
  case GSL_RNG_RAN2: T = gsl_rng_ran2; break;
  case GSL_RNG_RAN3: T = gsl_rng_ran3; break;
  case GSL_RNG_RANF: T = gsl_rng_ranf; break;
  case GSL_RNG_RANMAR: T = gsl_rng_ranmar; break;
  case GSL_RNG_R250: T = gsl_rng_r250; break;
  case GSL_RNG_TT800: T = gsl_rng_tt800; break;
  case GSL_RNG_VAX: T = gsl_rng_vax; break;
  case GSL_RNG_TRANSPUTER: T = gsl_rng_transputer; break;
  case GSL_RNG_RANDU: T = gsl_rng_randu; break;
  case GSL_RNG_MINSTD: T = gsl_rng_minstd; break;
  case GSL_RNG_UNI: T = gsl_rng_uni; break;
  case GSL_RNG_UNI32: T = gsl_rng_uni32; break;
  case GSL_RNG_SLATEC: T = gsl_rng_slatec; break;
  case GSL_RNG_ZUF: T = gsl_rng_zuf; break;
#ifdef GSL_1_1_LATER
  case GSL_RNG_BOROSH13: T = gsl_rng_borosh13; break;
  case GSL_RNG_COVEYOU: T = gsl_rng_coveyou; break;
  case GSL_RNG_FISHMAN18: T = gsl_rng_fishman18; break;
  case GSL_RNG_FISHMAN20: T = gsl_rng_fishman20; break;
  case GSL_RNG_FISHMAN2X: T = gsl_rng_fishman2x; break;
  case GSL_RNG_KNUTHRAN: T = gsl_rng_knuthran; break;
  case GSL_RNG_KNUTHRAN2: T = gsl_rng_knuthran2; break;
  case GSL_RNG_LECUYER21: T = gsl_rng_lecuyer21; break;
  case GSL_RNG_WATERMAN14: T = gsl_rng_waterman14; break;
#endif
#ifdef HAVE_RNGEXTRA_RNGEXTRA_H
  case GSL_RNGEXTRA_RNG1: T = rngextra_rng1; break;
  case GSL_RNGEXTRA_RNG2: T = rngextra_rng2; break;
#else
  case GSL_RNGEXTRA_RNG1:
  case GSL_RNGEXTRA_RNG2:
    rb_raise(rb_eNotImpError, "Install the rngextra package found at <http://www.network-theory.co.uk/download/rngextra/>.");
    break;
#endif
#ifdef GSL_1_9_LATER
  case GSL_RNG_KNUTHRAN2002: T = gsl_rng_knuthran2002; break;
#endif
  default:
    rb_raise(rb_eTypeError, "wrong generator type");
  }

  return  T;
}

static void rb_gsl_rng_define_const_type(VALUE module)
{
  rb_define_const(cgsl_rng, "DEFAULT", INT2FIX(GSL_RNG_DEFAULT));
  rb_define_const(cgsl_rng, "MT19937", INT2FIX(GSL_RNG_MT19937));
  rb_define_const(cgsl_rng, "MT19937_1999", INT2FIX(GSL_RNG_MT19937_1999));
  rb_define_const(cgsl_rng, "MT19937_1998", INT2FIX(GSL_RNG_MT19937_1999));
  rb_define_const(cgsl_rng, "RANLXS0", INT2FIX(GSL_RNG_RANLXS0));
  rb_define_const(cgsl_rng, "RANLXS1", INT2FIX(GSL_RNG_RANLXS1));
  rb_define_const(cgsl_rng, "RANLXS2", INT2FIX(GSL_RNG_RANLXS2));
  rb_define_const(cgsl_rng, "RANLXD1", INT2FIX(GSL_RNG_RANLXD1));
  rb_define_const(cgsl_rng, "RANLXD2", INT2FIX(GSL_RNG_RANLXD2));
  rb_define_const(cgsl_rng, "RANLUX", INT2FIX(GSL_RNG_RANLUX));
  rb_define_const(cgsl_rng, "RANLUX389", INT2FIX(GSL_RNG_RANLUX389));
  rb_define_const(cgsl_rng, "CMRG", INT2FIX(GSL_RNG_CMRG));
  rb_define_const(cgsl_rng, "MRG", INT2FIX(GSL_RNG_MRG));
  rb_define_const(cgsl_rng, "TAUS", INT2FIX(GSL_RNG_TAUS));
  rb_define_const(cgsl_rng, "TAUS2", INT2FIX(GSL_RNG_TAUS2));
  rb_define_const(cgsl_rng, "TAUS113", INT2FIX(GSL_RNG_TAUS113));
  rb_define_const(cgsl_rng, "GFSR4", INT2FIX(GSL_RNG_GFSR4));
  rb_define_const(cgsl_rng, "RAND", INT2FIX(GSL_RNG_RAND));
  rb_define_const(cgsl_rng, "RANDOM_LIBC5", INT2FIX(GSL_RNG_RANDOM_LIBC5));
  rb_define_const(cgsl_rng, "RANDOM8_LIBC5", INT2FIX(GSL_RNG_RANDOM8_LIBC5));
  rb_define_const(cgsl_rng, "RANDOM32_LIBC5", INT2FIX(GSL_RNG_RANDOM32_LIBC5));
  rb_define_const(cgsl_rng, "RANDOM64_LIBC5", INT2FIX(GSL_RNG_RANDOM64_LIBC5));
  rb_define_const(cgsl_rng, "RANDOM128_LIBC5", INT2FIX(GSL_RNG_RANDOM128_LIBC5));
  rb_define_const(cgsl_rng, "RANDOM256_LIBC5", INT2FIX(GSL_RNG_RANDOM256_LIBC5));
  rb_define_const(cgsl_rng, "RANDOM_GLIBC2", INT2FIX(GSL_RNG_RANDOM_GLIBC2));
  rb_define_const(cgsl_rng, "RANDOM8_GLIBC2", INT2FIX(GSL_RNG_RANDOM8_GLIBC2));
  rb_define_const(cgsl_rng, "RANDOM32_GLIBC2", INT2FIX(GSL_RNG_RANDOM32_GLIBC2));
  rb_define_const(cgsl_rng, "RANDOM64_GLIBC2", INT2FIX(GSL_RNG_RANDOM64_GLIBC2));
  rb_define_const(cgsl_rng, "RANDOM128_GLIBC2", INT2FIX(GSL_RNG_RANDOM128_GLIBC2));
  rb_define_const(cgsl_rng, "RANDOM256_GLIBC2", INT2FIX(GSL_RNG_RANDOM256_GLIBC2));
  rb_define_const(cgsl_rng, "RANDOM_BSD", INT2FIX(GSL_RNG_RANDOM_BSD));
  rb_define_const(cgsl_rng, "RANDOM8_BSD", INT2FIX(GSL_RNG_RANDOM8_BSD));
  rb_define_const(cgsl_rng, "RANDOM32_BSD", INT2FIX(GSL_RNG_RANDOM32_BSD));
  rb_define_const(cgsl_rng, "RANDOM64_BSD", INT2FIX(GSL_RNG_RANDOM64_BSD));
  rb_define_const(cgsl_rng, "RANDOM128_BSD", INT2FIX(GSL_RNG_RANDOM128_BSD));
  rb_define_const(cgsl_rng, "RANDOM256_BSD", INT2FIX(GSL_RNG_RANDOM256_BSD));
  rb_define_const(cgsl_rng, "RAND48", INT2FIX(GSL_RNG_RAND48));
  rb_define_const(cgsl_rng, "RAN0", INT2FIX(GSL_RNG_RAN0));
  rb_define_const(cgsl_rng, "RAN1", INT2FIX(GSL_RNG_RAN1));
  rb_define_const(cgsl_rng, "RAN2", INT2FIX(GSL_RNG_RAN2));
  rb_define_const(cgsl_rng, "RAN3", INT2FIX(GSL_RNG_RAN3));
  rb_define_const(cgsl_rng, "RANF", INT2FIX(GSL_RNG_RANF));
  rb_define_const(cgsl_rng, "RANMAR", INT2FIX(GSL_RNG_RANMAR));
  rb_define_const(cgsl_rng, "R250", INT2FIX(GSL_RNG_R250));
  rb_define_const(cgsl_rng, "TT800", INT2FIX(GSL_RNG_TT800));
  rb_define_const(cgsl_rng, "VAX", INT2FIX(GSL_RNG_VAX));
  rb_define_const(cgsl_rng, "TRANSPUTER", INT2FIX(GSL_RNG_TRANSPUTER));
  rb_define_const(cgsl_rng, "RANDU", INT2FIX(GSL_RNG_RANDU));
  rb_define_const(cgsl_rng, "MINSTD", INT2FIX(GSL_RNG_MINSTD));
  rb_define_const(cgsl_rng, "UNI", INT2FIX(GSL_RNG_UNI));
  rb_define_const(cgsl_rng, "UNI32", INT2FIX(GSL_RNG_UNI32));
  rb_define_const(cgsl_rng, "SLATEC", INT2FIX(GSL_RNG_SLATEC));
  rb_define_const(cgsl_rng, "ZUF", INT2FIX(GSL_RNG_ZUF));
  /* from gsl-1.1 */
  rb_define_const(cgsl_rng, "BOROSH13", INT2FIX(GSL_RNG_BOROSH13));
  rb_define_const(cgsl_rng, "COVEYOU", INT2FIX(GSL_RNG_COVEYOU));
  rb_define_const(cgsl_rng, "FISHMAN18", INT2FIX(GSL_RNG_FISHMAN18));
  rb_define_const(cgsl_rng, "FISHMAN20", INT2FIX(GSL_RNG_FISHMAN20));
  rb_define_const(cgsl_rng, "FISHMAN2X", INT2FIX(GSL_RNG_FISHMAN2X));
  rb_define_const(cgsl_rng, "KNUTHRAN", INT2FIX(GSL_RNG_KNUTHRAN));
  rb_define_const(cgsl_rng, "KNUTHRAN2", INT2FIX(GSL_RNG_KNUTHRAN2));
  rb_define_const(cgsl_rng, "LECUYER21", INT2FIX(GSL_RNG_LECUYER21));
  rb_define_const(cgsl_rng, "WATERMAN14", INT2FIX(GSL_RNG_WATERMAN14));
  rb_define_const(cgsl_rng, "RNGEXTRA_RNG1", INT2FIX(GSL_RNGEXTRA_RNG1));
  rb_define_const(cgsl_rng, "RNGEXTRA_RNG2", INT2FIX(GSL_RNGEXTRA_RNG2));
  rb_define_const(module, "RNGEXTRA_RNG1", INT2FIX(GSL_RNGEXTRA_RNG1));
  rb_define_const(module, "RNGEXTRA_RNG2", INT2FIX(GSL_RNGEXTRA_RNG2));

}

/* singleton */
static VALUE rb_gsl_rng_default_seed(VALUE obj)
{
  return UINT2NUM(gsl_rng_default_seed);
}

/* singleton */
static VALUE rb_gsl_rng_set_default_seed(VALUE obj, VALUE seed)
{
  gsl_rng_default_seed = NUM2UINT(seed);
  return seed;
}

static VALUE rb_gsl_rng_set(VALUE obj, VALUE s)
{
  gsl_rng *r = NULL;
  unsigned long seed;
  seed = NUM2UINT(s);
  Data_Get_Struct(obj, gsl_rng, r);
  gsl_rng_set(r, seed);
  return obj;
}

/*
  Document-method: <i>GSL::Rng#get</i>
    Returns a random integer from the generator.
    The minimum and maximum values depend on the algorithm used,
    but all integers in the range [min,max] are equally likely.
    The values of min and max can determined using the auxiliary
    methodss GSL::Rng#max and GSL::Rng#min.
*/
static VALUE rb_gsl_rng_get(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;
  gsl_vector_int *v;
  size_t n, i;
  Data_Get_Struct(obj, gsl_rng, r);
  switch (argc) {
  case 0:
    return UINT2NUM(gsl_rng_get(r));
    break;
  case 1:
    n = NUM2INT(argv[0]);
    v = gsl_vector_int_alloc(n);
    for (i = 0; i < n; i++) gsl_vector_int_set(v, i, (int) gsl_rng_get(r));
    return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free, v);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    break;
  }
}

static VALUE rb_gsl_rng_uniform(int argc, VALUE *argv, VALUE obj)
{
  gsl_rng *r = NULL;
  gsl_vector *v;
  size_t n, i;
  Data_Get_Struct(obj, gsl_rng, r);
  switch (argc) {
  case 0:
    return rb_float_new(gsl_rng_uniform(r));
    break;
  case 1:
    n = NUM2INT(argv[0]);
    v = gsl_vector_alloc(n);
    for (i = 0; i < n; i++) gsl_vector_set(v, i, gsl_rng_uniform(r));
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0 or 1)", argc);
    break;
  }
}

static VALUE rb_gsl_rng_uniform_pos(VALUE obj)
{
  gsl_rng *r = NULL;
  Data_Get_Struct(obj, gsl_rng, r);
  return rb_float_new(gsl_rng_uniform_pos(r));
}

static VALUE rb_gsl_rng_uniform_int(VALUE obj, VALUE n)
{
  gsl_rng *r = NULL;
  unsigned long int nn;
  nn = NUM2UINT(n);
  Data_Get_Struct(obj, gsl_rng, r);
  return UINT2NUM(gsl_rng_uniform_int(r, nn));
}

static VALUE rb_gsl_rng_name(VALUE obj)
{
  gsl_rng *r = NULL;
  Data_Get_Struct(obj, gsl_rng, r);
  return rb_str_new2(gsl_rng_name(r));
}

static VALUE rb_gsl_rng_max(VALUE obj, VALUE s)
{
  gsl_rng *r = NULL;
  Data_Get_Struct(obj, gsl_rng, r);
  return UINT2NUM(gsl_rng_max(r));
}

static VALUE rb_gsl_rng_min(VALUE obj, VALUE s)
{
  gsl_rng *r = NULL;
  Data_Get_Struct(obj, gsl_rng, r);
  return UINT2NUM(gsl_rng_min(r));
}

static VALUE rb_gsl_rng_size(VALUE obj, VALUE s)
{
  gsl_rng *r = NULL;
  Data_Get_Struct(obj, gsl_rng, r);
  return UINT2NUM(gsl_rng_size(r));
}

/* get a list of the available generator types */
/* module function */
static VALUE rb_gsl_rng_types_setup(VALUE obj)
{
  const gsl_rng_type **t, **t0;
  VALUE ary;
  t0 = gsl_rng_types_setup();
  ary = rb_ary_new();
  for (t = t0; *t != 0; t++) rb_ary_push(ary, rb_str_new2((*t)->name));
  return ary;
}

/* module function */
static VALUE rb_gsl_rng_env_setup(VALUE obj)
{
  gsl_rng_env_setup();
  return obj;
}

static VALUE rb_gsl_rng_clone(VALUE obj)
{
  gsl_rng *r = NULL, *rnew = NULL;
  Data_Get_Struct(obj, gsl_rng, r);
  rnew = gsl_rng_clone(r);
  return Data_Wrap_Struct(CLASS_OF(obj), 0, gsl_rng_free, rnew);
}

static VALUE rb_gsl_rng_print_state(VALUE obj)
{
  gsl_rng *r = NULL;
  Data_Get_Struct(obj, gsl_rng, r);
  gsl_rng_print_state(r);
  return obj;
}

#ifdef GSL_1_4_LATER
static VALUE rb_gsl_rng_fwrite(VALUE obj, VALUE io)
{
  gsl_rng *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_rng, h);
  f = rb_gsl_open_writefile(io, &flag);
  status = gsl_rng_fwrite(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE rb_gsl_rng_fread(VALUE obj, VALUE io)
{
  gsl_rng *h = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, gsl_rng, h);
  f = rb_gsl_open_readfile(io, &flag);
  status = gsl_rng_fread(f, h);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}
#endif

static VALUE rb_gsl_rng_memcpy(VALUE obj, VALUE dst, VALUE org)
{
  gsl_rng *dest, *src;
  CHECK_RNG(dst); CHECK_RNG(org);
  Data_Get_Struct(dst, gsl_rng, dest);
  Data_Get_Struct(org, gsl_rng, src);
  gsl_rng_memcpy(dest, src);
  return dst;
}

void Init_gsl_rng(VALUE module)
{
  cgsl_rng = rb_define_class_under(module, "Rng", cGSL_Object);

  rb_gsl_rng_define_const_type(module);

  rb_define_singleton_method(cgsl_rng, "alloc", rb_gsl_rng_alloc, -1);

  rb_define_singleton_method(cgsl_rng, "default_seed", rb_gsl_rng_default_seed, 0);
  rb_define_singleton_method(cgsl_rng, "set_default_seed", rb_gsl_rng_set_default_seed, 1);
  rb_define_singleton_method(cgsl_rng, "default_seed=", rb_gsl_rng_set_default_seed, 1);

  rb_define_method(cgsl_rng, "set", rb_gsl_rng_set, 1);
  rb_define_alias(cgsl_rng, "set_seed", "set");
  rb_define_alias(cgsl_rng, "seed=", "set");
  rb_define_method(cgsl_rng, "get", rb_gsl_rng_get, -1);
  rb_define_alias(cgsl_rng, "gen", "get");
  rb_define_method(cgsl_rng, "uniform", rb_gsl_rng_uniform, -1);
  rb_define_method(cgsl_rng, "uniform_pos", rb_gsl_rng_uniform_pos, 0);
  rb_define_method(cgsl_rng, "uniform_int", rb_gsl_rng_uniform_int, 1);

  rb_define_method(cgsl_rng, "name", rb_gsl_rng_name, 0);
  rb_define_method(cgsl_rng, "max", rb_gsl_rng_max, 0);
  rb_define_method(cgsl_rng, "min", rb_gsl_rng_min, 0);
  rb_define_method(cgsl_rng, "size", rb_gsl_rng_size, 0);

  rb_define_singleton_method(cgsl_rng, "types_setup", rb_gsl_rng_types_setup, 0);
  rb_define_singleton_method(cgsl_rng, "types", rb_gsl_rng_types_setup, 0);
  rb_define_singleton_method(cgsl_rng, "env_setup", rb_gsl_rng_env_setup, 0);

  rb_define_method(cgsl_rng, "clone", rb_gsl_rng_clone, 0);
  rb_define_alias(cgsl_rng, "duplicate", "clone");
  rb_define_method(cgsl_rng, "print_state", rb_gsl_rng_print_state, 0);

#ifdef GSL_1_4_LATER
  rb_define_method(cgsl_rng, "fwrite", rb_gsl_rng_fwrite, 1);
  rb_define_method(cgsl_rng, "fread", rb_gsl_rng_fread, 1);
#endif
  rb_define_singleton_method(cgsl_rng, "memcpy", rb_gsl_rng_memcpy, 2);
}
