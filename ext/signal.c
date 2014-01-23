/*
  signal.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_config.h"
#include "rb_gsl_fft.h"

enum FFT_CONV_CORR {
  RB_GSL_FFT_CONVOLVE = 0,
  RB_GSL_FFT_CORRELATE = 1,
  RB_GSL_FFT_REAL = 2,
  RB_GSL_FFT_HALFCOMPLEX = 3,
  RB_GSL_FFT_DECONVOLVE = 4,  
};

#ifndef WAVETABLE_P
#define WAVETABLE_P(x) (rb_obj_is_kind_of(x,cgsl_fft_halfcomplex_wavetable)?1:0)
#endif

#ifndef CHECK_WAVETABLE
#define CHECK_WAVETABLE(x) if(!rb_obj_is_kind_of(x,cgsl_fft_halfcomplex_wavetable))\
    rb_raise(rb_eTypeError, "wrong argument type (FFT::HalfComplex::Wavetable expected)");
#endif

#ifndef WORKSPACE_P
#define WORKSPACE_P(x) (rb_obj_is_kind_of(x,cgsl_fft_real_workspace)?1:0)
#endif

#ifndef CHECK_WORKSPACE
#define CHECK_WORKSPACE(x) if(!rb_obj_is_kind_of(x,cgsl_fft_real_workspace))\
    rb_raise(rb_eTypeError, "wrong argument type (FFT::Real::Workspace expected)");
#endif

static void complex_mul(double re1, double im1, double re2, double im2,
      double *re, double *im)
{
  *re = re1*re2 - im1*im2;
  *im = re1*im2 + im1*re2;
}

static void complex_conj_mul(double re1, double im1, double re2, double im2,
      double *re, double *im)
{
  *re = re1*re2 + im1*im2;
  *im = -re1*im2 + im1*re2;
}

static void complex_div(double re1, double im1, double re2, double im2,
      double *re, double *im)
{
  double factor = re2*re2 + im2*im2;
  complex_conj_mul(re1, im1, re2, im2, re, im);
  *re /= factor;
  *im /= factor;
}

/* data1, data2: FFTed data */
static void rbgsl_calc_conv_corr_c(const double *data1, const double *data2, 
     double *data3, size_t size,
     gsl_fft_halfcomplex_wavetable *table,
     gsl_fft_real_workspace *space, enum FFT_CONV_CORR calcflag)
{
  size_t i;
  double re1, re2, im1, im2;
  void (*complex_cal)(double, double, double, double, double*, double*);
  
  switch (calcflag) {
  case RB_GSL_FFT_CONVOLVE:
    complex_cal = complex_mul;
    data3[0] = data1[0]*data2[0];
    data3[size-1] = data1[size-1]*data2[size-1];    
    break;
  case RB_GSL_FFT_CORRELATE:  
    data3[0] = data1[0]*data2[0];  
    data3[size-1] = data1[size-1]*data2[size-1];        
    complex_cal = complex_conj_mul;
    break;
  case RB_GSL_FFT_DECONVOLVE:
    complex_cal = complex_div;
    data3[0] = data1[0]/data2[0];
    data3[size-1] = data1[size-1]/data2[size-1];       
    break;
  default:
    rb_raise(rb_eArgError, "Wrong flag.");
    break;
  }
  
  for (i = 1; i < size-1; i+=2) {
    re1 = data1[i];    im1 = data1[i+1];
    re2 = data2[i];    im2 = data2[i+1];
    (*complex_cal)(re1, im1, re2, im2, &data3[i], &data3[i+1]);
  }

}

static VALUE rb_gsl_fft_conv_corr(int argc, VALUE *argv, VALUE obj,
        enum FFT_CONV_CORR flag1,
        enum FFT_CONV_CORR flag2)
{
  double *data1, *data2, *data3 = NULL;
  size_t stride1, stride2, size1, size2;
#ifdef HAVE_NARRAY_H
  int naflag1, naflag2, shape;
#else
  int naflag1, naflag2;
#endif
  gsl_vector *v = NULL;
  gsl_fft_halfcomplex_wavetable *table = NULL;
  gsl_fft_real_wavetable *rtable = NULL;
  gsl_fft_real_workspace *space = NULL, *space2 = NULL;
  int flagt = 0, flagw = 0;
  //  size_t i;
  gsl_vector *vtmp1 = NULL, *vtmp2 = NULL;
  VALUE ary = NULL;
  switch (argc) {
  case 3:
    data1 = get_ptr_double3(obj, &size1, &stride1, &naflag1);
    data2 = get_ptr_double3(argv[0], &size2, &stride2, &naflag2);
    CHECK_WAVETABLE(argv[1]);
    Data_Get_Struct(argv[1], gsl_fft_halfcomplex_wavetable, table);
    CHECK_WORKSPACE(argv[2]);
    Data_Get_Struct(argv[2], gsl_fft_real_workspace, space);
    break;
  case 2:
    data1 = get_ptr_double3(obj, &size1, &stride1, &naflag1);
    data2 = get_ptr_double3(argv[0], &size2, &stride2, &naflag2);
    if (WAVETABLE_P(argv[1])) {
      Data_Get_Struct(argv[1], gsl_fft_halfcomplex_wavetable, table);
      space = gsl_fft_real_workspace_alloc(size1);
      flagw = 1;
    } else if (WORKSPACE_P(argv[1])) {
      Data_Get_Struct(argv[1], gsl_fft_real_workspace, space);
      table = gsl_fft_halfcomplex_wavetable_alloc(size1);
      flagt = 1;
    } else {
      rb_raise(rb_eTypeError, 
         "wrong argument type %s "
         "(FFT::HalfComplex::Wavetable or FFT::Real::Workspace expected)", 
         rb_class2name(CLASS_OF(argv[2])));
    }
    break;
  case 1:
    data1 = get_ptr_double3(obj, &size1, &stride1, &naflag1);
    data2 = get_ptr_double3(argv[0], &size2, &stride2, &naflag2);
    table = gsl_fft_halfcomplex_wavetable_alloc(size1);
    space = gsl_fft_real_workspace_alloc(size1);
    flagt = 1;
    flagw = 1;
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1-3)", argc);
  }

  switch (naflag1*naflag2) {
  case 0:
    v = gsl_vector_alloc(size1);
    switch (flag1) {
    case RB_GSL_FFT_REAL:
      ary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    default:      
      ary = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
      break;
    }
    data3 = v->data;
    break;
  case 1:
#ifdef HAVE_NARRAY_H
    shape = (int) size1;
    ary = na_make_object(NA_DFLOAT, 1, &shape, cNArray);
    data3 = NA_PTR_TYPE(ary, double*);
#endif
    break;
  default:
    break;
  }
  switch (flag1) {
  case RB_GSL_FFT_REAL:  // do FFT
    vtmp1 = gsl_vector_alloc(size1);
    vtmp2 = gsl_vector_alloc(size2);
    memcpy(vtmp1->data, data1, sizeof(double)*size1);
    memcpy(vtmp2->data, data2, sizeof(double)*size2);
    data1 = vtmp1->data;
    data2 = vtmp2->data;
    rtable = gsl_fft_real_wavetable_alloc(size1);
    if (size1 == space->n) {
      gsl_fft_real_transform(data1, stride1, size1, rtable, space); 
    } else {
      space2 = gsl_fft_real_workspace_alloc(size1);
      gsl_fft_real_transform(data1, stride1, size1, rtable, space2); 
      /* no freeing space2 here */
    }
    if (size1 != size2) {
      if (rtable) gsl_fft_real_wavetable_free(rtable);
      rtable = gsl_fft_real_wavetable_alloc(size2);
    }
    if (size2 == space->n) {
      gsl_fft_real_transform(data2, stride2, size2, rtable, space); 
    } else if (size2 == size1) {
      gsl_fft_real_transform(data2, stride2, size2, rtable, space2); 
      gsl_fft_real_workspace_free(space2);
    } else {
      if (space2) gsl_fft_real_workspace_free(space2);
      space2 = gsl_fft_real_workspace_alloc(size2);
      gsl_fft_real_transform(data2, stride2, size2, rtable, space2); 
      gsl_fft_real_workspace_free(space2);
    }
    gsl_fft_real_wavetable_free(rtable);
    space2 = NULL;
    rtable = NULL;
    break;
  case RB_GSL_FFT_HALFCOMPLEX:
    /* do nothing */
    break;
  default:
    /* not occur */
    break;
  }
  
  rbgsl_calc_conv_corr_c(data1, data2, data3, size1, table, space, flag2);

  if (flag1 == RB_GSL_FFT_REAL) {
    gsl_fft_halfcomplex_inverse(data3, 1, size1, table, space);
//    for (i = 0; i < size1; i++) data3[i] /= size1;
  }
  
  if (flagt == 1) gsl_fft_halfcomplex_wavetable_free(table);
  if (flagw == 1) gsl_fft_real_workspace_free(space);
  if (vtmp1) gsl_vector_free(vtmp1);
  if (vtmp2) gsl_vector_free(vtmp2);
  return ary;
}

/* GSL::Vector#convolve */
static VALUE rb_gsl_fft_real_convolve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_fft_conv_corr(argc, argv, obj, 
          RB_GSL_FFT_REAL,
          RB_GSL_FFT_CONVOLVE);
}
/* GSL::Vector#deconvolve */
static VALUE rb_gsl_fft_real_deconvolve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_fft_conv_corr(argc, argv, obj, 
            RB_GSL_FFT_REAL,
            RB_GSL_FFT_DECONVOLVE);
}

/* GSL::Vector#correlate */
static VALUE rb_gsl_fft_real_correlate(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_fft_conv_corr(argc, argv, obj, 
            RB_GSL_FFT_REAL,
            RB_GSL_FFT_CORRELATE);
}

/* GSL::Vector#halfcomplex_convolve */
static VALUE rb_gsl_fft_halfcomplex_convolve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_fft_conv_corr(argc, argv, obj, 
          RB_GSL_FFT_HALFCOMPLEX,
          RB_GSL_FFT_CONVOLVE);
}
/* GSL::Vector#halfcomplex_deconvolve */
static VALUE rb_gsl_fft_halfcomplex_deconvolve(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_fft_conv_corr(argc, argv, obj, 
            RB_GSL_FFT_HALFCOMPLEX,
            RB_GSL_FFT_DECONVOLVE);
}
/* GSL::Vector#halfcomplex_correlate */
static VALUE rb_gsl_fft_halfcomplex_correlate(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_fft_conv_corr(argc, argv, obj, 
            RB_GSL_FFT_HALFCOMPLEX,
            RB_GSL_FFT_CORRELATE);
}

void Init_gsl_signal(VALUE module)
{
  rb_define_method(cgsl_vector, "real_convolve", rb_gsl_fft_real_convolve, -1);
  rb_define_method(cgsl_vector, "real_deconvolve", rb_gsl_fft_real_deconvolve, -1);			     
  rb_define_method(cgsl_vector, "real_correlate", rb_gsl_fft_real_correlate, -1);	

  rb_define_alias(cgsl_vector, "convolve", "real_convolve");
  rb_define_alias(cgsl_vector, "deconvolve", "real_deconvolve");
  rb_define_alias(cgsl_vector, "correlate", "real_correlate");

  rb_define_method(cgsl_vector, "halfcomplex_convolve", rb_gsl_fft_halfcomplex_convolve, -1);
  rb_define_method(cgsl_vector, "halfcomplex_deconvolve", rb_gsl_fft_halfcomplex_deconvolve, -1);			     
  rb_define_method(cgsl_vector, "halfcomplex_correlate", rb_gsl_fft_halfcomplex_correlate, -1);			     

  rb_define_alias(cgsl_vector, "hc_convolve", "halfcomplex_convolve");
  rb_define_alias(cgsl_vector, "hc_deconvolve", "halfcomplex_deconvolve");
  rb_define_alias(cgsl_vector, "hc_correlate", "halfcomplex_correlate");
}

#undef WAVETABLE_P
#undef CHECK_WAVETABLE
#undef WORKSPACE_P
#undef CHECK_WORKSPACE
