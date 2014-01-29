/*
  rb_gsl_fft.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FFTNESS FOR A PARTICULAR PURPOSE.
*/

#ifndef ___RB_GSL_FFT_H___
#define ___RB_GSL_FFT_H___

#include <gsl/gsl_fft.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include "rb_gsl.h"
#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif

typedef struct
{
  size_t n;
  size_t nf;
  size_t factor[64];
  gsl_complex *twiddle[64];
  gsl_complex *trig;
} GSL_FFT_Wavetable;

typedef struct
{
  size_t n;
  double *scratch;
} GSL_FFT_Workspace;

enum {
  RB_GSL_FFT_INPLACE,
  RB_GSL_FFT_COPY,
};

EXTERN VALUE mgsl_fft;
EXTERN VALUE cgsl_fft_wavetable;
EXTERN VALUE cgsl_fft_wavetable_factor;
EXTERN VALUE cgsl_fft_complex_wavetable, cgsl_fft_complex_workspace;
EXTERN VALUE cgsl_fft_real_wavetable, cgsl_fft_halfcomplex_wavetable;
EXTERN VALUE cgsl_fft_real_workspace;

#endif
