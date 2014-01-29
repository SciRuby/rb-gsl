/*
  stats.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_array.h"
#include "rb_gsl_statistics.h"
#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif

static double* get_vector_stats2(int argc, VALUE *argv, VALUE obj, 
				 size_t *stride, size_t *size)
{
  double *v = NULL;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments");
    v = get_vector_ptr(argv[0], stride, size);
    break;
  default:
    v = get_vector_ptr(obj, stride, size);
    break;
  }
  return v;
}

static VALUE rb_gsl_stats_XXX(int argc, VALUE *argv, VALUE obj,
			      double (*f)(const double*, size_t, size_t))
{
  size_t stride, size;
  double *data = NULL;
  data = get_vector_stats2(argc, argv, obj, &stride, &size);
  return rb_float_new((*f)(data, stride, size));
}

static VALUE rb_gsl_stats_XXX1(int argc, VALUE *argv, VALUE obj, 
			      double (*f)(const double*, size_t, size_t, double))
{
  size_t stride, size;
  double *data = NULL;
  if (argc < 1) rb_raise(rb_eArgError, "too few arguments");
  data = get_vector_stats2(argc, argv, obj, &stride, &size);
  return rb_float_new((*f)(data, stride, size, NUM2DBL(argv[argc-1])));
}

static VALUE rb_gsl_stats_XXX2(int argc, VALUE *argv, VALUE obj, 
			       double (*f)(const double*, size_t, size_t),
			       double (*fm)(const double*, size_t, size_t, double, double))
{
  double x, a, b, *data = NULL;
  size_t stride, size;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    switch (argc) {
    case 2:
      data = get_vector_ptr(argv[0], &stride, &size);
      a = NUM2DBL(argv[1]); b = NUM2DBL(argv[2]);
      x = (*fm)(data, stride, size, a, b);
      break;
    case 1:
      data = get_vector_ptr(argv[0], &stride, &size);
      x = (*f)(data, stride, size);
      break;
    default:
      rb_raise(rb_eArgError, 
	       "wrong number of arguments (%d for 1 or 2)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 0:
      data = get_vector_ptr(obj, &stride, &size);
      x = (*f)(data, stride, size);
      break;
    case 1:
      data = get_vector_ptr(obj, &stride, &size);
      a = NUM2DBL(argv[0]); b = NUM2DBL(argv[1]);
      x = (*fm)(data, stride, size, a, b);
      break;
    default:
      rb_raise(rb_eArgError, 
	       "wrong number of arguments (%d for 0 or 1)", argc);
      break;
    }
    break;
  }
  return rb_float_new(x);
}

static VALUE rb_gsl_stats_mean(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_XXX(argc, argv, obj, gsl_stats_mean);
}

static VALUE rb_gsl_stats_XXX_m(int argc, VALUE *argv, VALUE obj,
				double (*f)(const double*, size_t, size_t),
				double (*fm)(const double*, size_t, size_t, double))
{
  double x, mean, *data = NULL;
  size_t stride, size;
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    switch (argc) {
    case 2:
      data = get_vector_ptr(argv[0], &stride, &size);
      mean = NUM2DBL(argv[1]);
      x = (*fm)(data, stride, size, mean);
      break;
    case 1:
      data = get_vector_ptr(argv[0], &stride, &size);
      x = (*f)(data, stride, size);
      break;
    default:
      rb_raise(rb_eArgError, 
	       "wrong number of arguments (%d for 1 or 2)", argc);
      break;
    }
    break;
  default:
    switch (argc) {
    case 0:
      data = get_vector_ptr(obj, &stride, &size);
      x = (*f)(data, stride, size);
      break;
    case 1:
      data = get_vector_ptr(obj, &stride, &size);
      mean = NUM2DBL(argv[0]);
      x = (*fm)(data, stride, size, mean);
      break;
    default:
      rb_raise(rb_eArgError, 
	       "wrong number of arguments (%d for 0 or 1)", argc);
      break;
    }
    break;
  }
  return rb_float_new(x);
}

static VALUE rb_gsl_stats_variance_m(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_XXX_m(argc, argv, obj,
			    gsl_stats_variance, gsl_stats_variance_m);
}

static VALUE rb_gsl_stats_sd_m(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_XXX_m(argc, argv, obj,
			    gsl_stats_sd, gsl_stats_sd_m);
}

#ifdef GSL_1_11_LATER
static VALUE rb_gsl_stats_tss_m(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_XXX_m(argc, argv, obj,
			    gsl_stats_tss, gsl_stats_tss_m);
}
#endif

static VALUE rb_gsl_stats_variance_with_fixed_mean(int argc, VALUE *argv, 
						   VALUE obj)
{
  return rb_gsl_stats_XXX1(argc, argv, obj,
			   gsl_stats_variance_with_fixed_mean);
}

static VALUE rb_gsl_stats_sd_with_fixed_mean(int argc, VALUE *argv, 
						   VALUE obj)
{
  return rb_gsl_stats_XXX1(argc, argv, obj,
			   gsl_stats_sd_with_fixed_mean);
}

static VALUE rb_gsl_stats_absdev_m(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_XXX_m(argc, argv, obj,
			    gsl_stats_absdev, gsl_stats_absdev_m);
}

static VALUE rb_gsl_stats_skew(int argc, VALUE *argv, 
			       VALUE obj)
{
  return rb_gsl_stats_XXX2(argc, argv, obj,
			   gsl_stats_skew,
			   gsl_stats_skew_m_sd);
}

static VALUE rb_gsl_stats_kurtosis(int argc, VALUE *argv, 
						   VALUE obj)
{
  return rb_gsl_stats_XXX2(argc, argv, obj,
			   gsl_stats_kurtosis,
			   gsl_stats_kurtosis_m_sd);
}

static VALUE rb_gsl_stats_lag1_autocorrelation(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_XXX_m(argc, argv, obj,
			    gsl_stats_lag1_autocorrelation, gsl_stats_lag1_autocorrelation_m);
}

/****************************/

static void get_vector_stats3(int argc, VALUE *argv, VALUE obj, 
			      double **w, size_t *stridew, size_t *sizew,
			      double **x, size_t *stridex, size_t *sizex)
{
  switch (TYPE(obj)) {
  case T_MODULE:  case T_CLASS:  case T_OBJECT:
    if (argc < 2) rb_raise(rb_eArgError, "too few arguments");
    *w = get_vector_ptr(argv[0], stridew, sizew);
    *x = get_vector_ptr(argv[1], stridex, sizex);
    break;
  default:
    if (argc < 1) rb_raise(rb_eArgError, "too few arguments");
    *x = get_vector_ptr(obj, stridex, sizex);
    *w = get_vector_ptr(argv[0], stridew, sizew);
    break;
  }
}

static VALUE rb_gsl_stats_wXXX(int argc, VALUE *argv, VALUE obj,
			       double (*f)(const double*, size_t, const double*, 
					   size_t, size_t))
{
  double *w, *x;
  double mean;
  size_t sizew, stridew, sizex, stridex;
  get_vector_stats3(argc, argv, obj, &w, &stridew, &sizew, &x, &stridex, &sizex);
  mean = (*f)(w, stridew, x, stridex, sizex);
  return rb_float_new(mean);
}

static VALUE rb_gsl_stats_wXXX_m(int argc, VALUE *argv, VALUE obj,
			       double (*f)(const double*, size_t, const double*, 
					   size_t, size_t, double))
{
  double *w, *x;
  double mean;
  size_t sizew, stridew, sizex, stridex;
  get_vector_stats3(argc, argv, obj, &w, &stridew, &sizew, &x, &stridex, &sizex);
  mean = (*f)(w, stridew, x, stridex, sizex, NUM2DBL(argv[argc-1]));
  return rb_float_new(mean);
}

static VALUE rb_gsl_stats_wmean(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_wXXX(argc, argv, obj, gsl_stats_wmean);
}

static VALUE rb_gsl_stats_wvariance(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_wXXX(argc, argv, obj, gsl_stats_wvariance);
}

static VALUE rb_gsl_stats_wsd(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_wXXX(argc, argv, obj, gsl_stats_wsd);
}

static VALUE rb_gsl_stats_wvariance_m(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_wXXX_m(argc, argv, obj, gsl_stats_wvariance_m);
}

static VALUE rb_gsl_stats_wsd_m(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_wXXX_m(argc, argv, obj, gsl_stats_wsd_m);
}

static VALUE rb_gsl_stats_wvariance_with_fixed_mean(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_wXXX_m(argc, argv, obj, gsl_stats_wvariance_with_fixed_mean);
}

static VALUE rb_gsl_stats_wsd_with_fixed_mean(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_wXXX_m(argc, argv, obj, gsl_stats_wsd_with_fixed_mean);
}

static VALUE rb_gsl_stats_wabsdev(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_wXXX(argc, argv, obj, gsl_stats_wabsdev);
}

static VALUE rb_gsl_stats_wabsdev_m(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_wXXX_m(argc, argv, obj, gsl_stats_wabsdev_m);
}

static VALUE rb_gsl_stats_wskew(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_wXXX(argc, argv, obj, gsl_stats_wskew);
}

static VALUE rb_gsl_stats_wkurtosis(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_wXXX(argc, argv, obj, gsl_stats_wkurtosis);
}

static VALUE rb_gsl_stats_wskew_m_sd(VALUE obj, VALUE ww, VALUE wm, VALUE wsd)
{
  double *w, *x;
  size_t stridew, stridex, sizew, sizex;
  double skew_m;
  x = get_vector_ptr(obj, &stridex, &sizex);
  w = get_vector_ptr(ww, &stridew, &sizew);
  skew_m = gsl_stats_wskew_m_sd(w, stridew, x, stridex, sizex, NUM2DBL(wm),
				NUM2DBL(wsd));
  return rb_float_new(skew_m);
}

static VALUE rb_gsl_stats_wkurtosis_m_sd(VALUE obj, VALUE ww, VALUE wm, VALUE wsd)
{
  double *w, *x;
  size_t stridew, stridex, sizew, sizex;
  double kurtosis_m;
  x = get_vector_ptr(obj, &stridex, &sizex);
  w = get_vector_ptr(ww, &stridew, &sizew);
  kurtosis_m = gsl_stats_wkurtosis_m_sd(w, stridew, x, stridex, sizex, NUM2DBL(wm),
					NUM2DBL(wsd));
  return rb_float_new(kurtosis_m);
}

static VALUE rb_gsl_stats_max(int argc, VALUE *argv, VALUE obj)
{
  size_t stride, size;
  double max, *data = NULL;
  data = get_vector_stats2(argc, argv, obj, &stride, &size);
  max = gsl_stats_max(data, stride, size);
  return rb_float_new(max);
}

static VALUE rb_gsl_stats_min(int argc, VALUE *argv, VALUE obj)
{
  double min, *data = NULL;
  size_t stride, size;
  data = get_vector_stats2(argc, argv, obj, &stride, &size);
  min = gsl_stats_min(data, stride, size);
  return rb_float_new(min);
}

static VALUE rb_gsl_stats_minmax(int argc, VALUE *argv, VALUE obj)
{
  size_t stride, size;
  double min, max, *data = NULL;
  data = get_vector_stats2(argc, argv, obj, &stride, &size);
  gsl_stats_minmax(&min, &max, data, stride, size);
  return rb_ary_new3(2, rb_float_new(min), rb_float_new(max));
}

static VALUE rb_gsl_stats_max_index(int argc, VALUE *argv, VALUE obj)
{
  double *data = NULL;
  size_t index, stride, size;
  data = get_vector_stats2(argc, argv, obj, &stride, &size);
  index = gsl_stats_max_index(data, stride, size);
  return INT2FIX(index);
}

static VALUE rb_gsl_stats_min_index(int argc, VALUE *argv, VALUE obj)
{
  double *data = NULL;
  size_t index, stride, size;
  data = get_vector_stats2(argc, argv, obj, &stride, &size);
  index = gsl_stats_min_index(data, stride, size);
  return INT2FIX(index);
}

static VALUE rb_gsl_stats_minmax_index(int argc, VALUE *argv, VALUE obj)
{
  double *data = NULL;
  size_t imin, imax, stride, size;
  data = get_vector_stats2(argc, argv, obj, &stride, &size);
  gsl_stats_minmax_index(&imin, &imax, data, stride, size);
  return rb_ary_new3(2, INT2FIX(imin), INT2FIX(imax));
}

static VALUE rb_gsl_stats_median_from_sorted_data(int argc, VALUE *argv, VALUE obj)
{
  size_t stride, size;
  double median, *data = NULL;
  data = get_vector_stats2(argc, argv, obj, &stride, &size);
  median = gsl_stats_median_from_sorted_data(data, stride, size);
  return rb_float_new(median);
}

static VALUE rb_gsl_stats_median(int argc, VALUE *argv, VALUE obj)
{
  size_t stride, size;
  double median, *data = NULL, *data2 = NULL;
  data = get_vector_stats2(argc, argv, obj, &stride, &size);
  data2 = (double *) malloc(sizeof(double)*size*stride);
  memcpy(data2, data, sizeof(double)*size*stride);
  gsl_sort(data2, stride, size);
  median = gsl_stats_median_from_sorted_data(data2, stride, size);
  free(data2);
  return rb_float_new(median);
}

static VALUE rb_gsl_stats_quantile_from_sorted_data(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_stats_XXX1(argc, argv, obj,
			   gsl_stats_quantile_from_sorted_data);
  /*  size_t stride, size;
  double quantile, *data = NULL;
  data = get_vector_ptr(obj, &stride, &size);
  quantile = gsl_stats_quantile_from_sorted_data(data, stride, size, NUM2DBL(f));
  return rb_float_new(quantile);*/
}
/*
static VALUE rb_gsl_stats_quantile(int argc, VALUE *argv, VALUE obj)
{
  size_t stride, size;
  double quantile, *data = NULL, *data2 = NULL;
  data = get_vector_stats2(argc, argv, obj, &stride, &size);
  data2 = (double *) malloc(sizeof(double)*size*stride);
  memcpy(data2, data, sizeof(double)*size*stride);
  gsl_sort(data2, stride, size);
  quantile = gsl_stats_quantile_from_sorted_data(data2, stride, size);
  free(data2);
  return rb_float_new(quantile);
}
*/
static VALUE rb_gsl_stats_covariance2(VALUE obj, VALUE vv1, VALUE vv2)
{
  double *data1, *data2;
  size_t stride1, stride2, size;
  data1 = get_vector_ptr(vv1, &stride1, &size);
  data2 = get_vector_ptr(vv2, &stride2, &size);
  return rb_float_new(gsl_stats_covariance(data1, stride1, data2,
					   stride2, size));
}

static VALUE rb_gsl_stats_covariance_m2(VALUE obj, VALUE vv1, VALUE vv2,
					VALUE m1, VALUE m2)
{
  double *data1, *data2;
  size_t stride1, stride2, size;
  data1 = get_vector_ptr(vv1, &stride1, &size);
  data2 = get_vector_ptr(vv2, &stride2, &size);
  return rb_float_new(gsl_stats_covariance_m(data1, stride1, data2,
					   stride2, size, NUM2DBL(m1), NUM2DBL(m2)));
}

#ifdef GSL_1_10_LATER
static VALUE rb_gsl_stats_correlation(VALUE obj, VALUE vv1, VALUE vv2)
{
  double *data1, *data2;
  size_t stride1, stride2, size;
  data1 = get_vector_ptr(vv1, &stride1, &size);
  data2 = get_vector_ptr(vv2, &stride2, &size);
  return rb_float_new(gsl_stats_correlation(data1, stride1, data2,
					   stride2, size));
}
static VALUE rb_gsl_stats_pvariance(VALUE obj, VALUE vv1, VALUE vv2)
{
  double *data1, *data2;
  size_t stride1, stride2, size1, size2;
  data1 = get_vector_ptr(vv1, &stride1, &size1);
  data2 = get_vector_ptr(vv2, &stride2, &size2);
  return rb_float_new(gsl_stats_pvariance(data1, stride1, size1, data2,
					   stride2, size2));
}
static VALUE rb_gsl_stats_ttest(VALUE obj, VALUE vv1, VALUE vv2)
{
  double *data1, *data2;
  size_t stride1, stride2, size1, size2;
  data1 = get_vector_ptr(vv1, &stride1, &size1);
  data2 = get_vector_ptr(vv2, &stride2, &size2);
  return rb_float_new(gsl_stats_ttest(data1, stride1, size1, data2,
					   stride2, size2));
}
#endif

static VALUE rb_gsl_stats_wmean2(VALUE obj, VALUE ww, VALUE dd)
{
  double wmean, *dataw = NULL, *data = NULL;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  wmean = gsl_stats_wmean(dataw, stridew, data, strided, size);
  return rb_float_new(wmean);
}

static VALUE rb_gsl_stats_wvariance2(VALUE obj, VALUE ww, VALUE dd)
{
  double wvariance, *dataw = NULL, *data = NULL;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  wvariance = gsl_stats_wvariance(dataw, stridew, data, strided, size);
  return rb_float_new(wvariance);
}

static VALUE rb_gsl_stats_wvariance_m2(VALUE obj, VALUE ww, VALUE dd, VALUE mm)
{
  double *dataw = NULL, *data = NULL;
  double wvariance, m;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  m = NUM2DBL(mm);
  wvariance = gsl_stats_wvariance_m(dataw, stridew, data, strided, size, m);
  return rb_float_new(wvariance);
}

static VALUE rb_gsl_stats_wsd2(VALUE obj, VALUE ww, VALUE dd)
{
  double wsd, *dataw = NULL, *data = NULL;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  wsd = gsl_stats_wsd(dataw, stridew, data, strided, size);
  return rb_float_new(wsd);
}

static VALUE rb_gsl_stats_wsd_m2(VALUE obj, VALUE ww, VALUE dd, VALUE mm)
{
  double *dataw = NULL, *data = NULL;
  double wsd, m;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  m = NUM2DBL(mm);
  wsd = gsl_stats_wsd_m(dataw, stridew, data, strided, size, m);
  return rb_float_new(wsd);
}

static VALUE rb_gsl_stats_wvariance_with_fixed_mean2(VALUE obj, VALUE ww, VALUE dd, 
						     VALUE mm)
{
  double wvariance, m;
  double *dataw = NULL, *data = NULL;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  m = NUM2DBL(mm);
  wvariance = gsl_stats_wvariance_with_fixed_mean(dataw, stridew, 
						  data, strided, size, m);

  return rb_float_new(wvariance);
}

static VALUE rb_gsl_stats_wsd_with_fixed_mean2(VALUE obj, VALUE ww, VALUE dd, 
						     VALUE mm)
{
  double wsd, m;
  double *dataw = NULL, *data = NULL;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  m = NUM2DBL(mm);
  wsd = gsl_stats_wsd_with_fixed_mean(dataw, stridew, data, strided, size, m);
  return rb_float_new(wsd);
}


static VALUE rb_gsl_stats_wabsdev2(VALUE obj, VALUE ww, VALUE dd)
{
  double wabsdev, *dataw = NULL, *data = NULL;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  wabsdev = gsl_stats_wabsdev(dataw, stridew, data, strided, size);
  return rb_float_new(wabsdev);
}

static VALUE rb_gsl_stats_wabsdev_m2(VALUE obj, VALUE ww, VALUE dd, VALUE mm)
{
  double *dataw = NULL, *data = NULL;
  double wabsdev, m;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  m = NUM2DBL(mm);
  wabsdev = gsl_stats_wabsdev_m(dataw, stridew, data, strided, size, m);
  return rb_float_new(wabsdev);
}
static VALUE rb_gsl_stats_wskew2(VALUE obj, VALUE ww, VALUE dd)
{
  double wskew, *dataw = NULL, *data = NULL;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  wskew = gsl_stats_wskew(dataw, stridew, data, strided, size);
  return rb_float_new(wskew);
}

static VALUE rb_gsl_stats_wskew_m2(VALUE obj, VALUE ww, VALUE dd, VALUE mm, VALUE ss)
{
  double *dataw = NULL, *data = NULL;
  double wskew, m, sd;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  m = NUM2DBL(mm);
  sd = NUM2DBL(ss);
  wskew = gsl_stats_wskew_m_sd(dataw, stridew, data, strided, size, m, sd);
  return rb_float_new(wskew);
}

static VALUE rb_gsl_stats_wkurtosis2(VALUE obj, VALUE ww, VALUE dd)
{
  double wkurtosis, *dataw = NULL, *data = NULL;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  wkurtosis = gsl_stats_wkurtosis(dataw, stridew, data, strided, size);
  return rb_float_new(wkurtosis);
}

static VALUE rb_gsl_stats_wkurtosis_m2(VALUE obj, VALUE ww, VALUE dd, VALUE mm, VALUE ss)
{
  double *dataw = NULL, *data = NULL;
  double wkurtosis, m, sd;
  size_t stridew, strided, size;
  dataw = get_vector_ptr(ww, &stridew, &size);
  data = get_vector_ptr(dd, &strided, &size);
  m = NUM2DBL(mm);
  sd = NUM2DBL(ss);
  wkurtosis = gsl_stats_wkurtosis_m_sd(dataw, stridew, data, strided, size, m, sd);
  return rb_float_new(wkurtosis);
}

void Init_gsl_stats(VALUE module)
{
  VALUE mgsl_stats;

  mgsl_stats = rb_define_module_under(module, "Stats");

  rb_define_singleton_method(mgsl_stats, "mean", rb_gsl_stats_mean, -1);
  rb_define_method(cgsl_vector, "stats_mean", rb_gsl_stats_mean, -1);
  rb_define_alias(cgsl_vector, "mean", "stats_mean");
  rb_define_alias(cgsl_vector, "average", "stats_mean");

  rb_define_singleton_method(mgsl_stats, "variance", rb_gsl_stats_variance_m, -1);
  rb_define_singleton_method(mgsl_stats, "variance_m", rb_gsl_stats_variance_m, -1);
  rb_define_method(cgsl_vector, "stats_variance_m", rb_gsl_stats_variance_m, -1);
  rb_define_alias(cgsl_vector, "variance_m", "stats_variance_m");
  rb_define_alias(cgsl_vector, "variance", "stats_variance_m");
  rb_define_alias(cgsl_vector, "var", "stats_variance_m");

  rb_define_singleton_method(mgsl_stats, "sd", rb_gsl_stats_sd_m, -1);
  rb_define_singleton_method(mgsl_stats, "sd_m", rb_gsl_stats_sd_m, -1);
#ifdef GSL_1_11_LATER
  rb_define_singleton_method(mgsl_stats, "tss", rb_gsl_stats_tss_m, -1);
  rb_define_singleton_method(mgsl_stats, "tss_m", rb_gsl_stats_tss_m, -1);
#endif
  rb_define_singleton_method(mgsl_stats, "sdev", rb_gsl_stats_sd_m, -1);
  rb_define_singleton_method(mgsl_stats, "sigma", rb_gsl_stats_sd_m, -1);
  rb_define_method(cgsl_vector, "stats_sd_m", rb_gsl_stats_sd_m, -1);
  rb_define_alias(cgsl_vector, "sd_m", "stats_sd_m");
  rb_define_alias(cgsl_vector, "stats_sd", "stats_sd_m");
  rb_define_alias(cgsl_vector, "sd", "stats_sd_m");
  rb_define_alias(cgsl_vector, "sigma", "stats_sd_m");
  rb_define_alias(cgsl_vector, "sdev", "stats_sd_m");
#ifdef GSL_1_11_LATER
  rb_define_method(cgsl_vector, "stats_tss_m", rb_gsl_stats_tss_m, -1);
  rb_define_alias(cgsl_vector, "stats_tss", "stats_tss_m");
  rb_define_alias(cgsl_vector, "tss_m", "stats_tss_m");
  rb_define_alias(cgsl_vector, "tss", "stats_tss_m");
#endif

  rb_define_singleton_method(mgsl_stats, "variance_with_fixed_mean", 
			     rb_gsl_stats_variance_with_fixed_mean, -1);
  rb_define_method(cgsl_vector, "stats_variance_with_fixed_mean", 
		   rb_gsl_stats_variance_with_fixed_mean, -1);
  rb_define_alias(cgsl_vector, "variance_with_fixed_mean", 
		  "stats_variance_with_fixed_mean");

  rb_define_singleton_method(mgsl_stats, "sd_with_fixed_mean", 
			     rb_gsl_stats_sd_with_fixed_mean, -1);
  rb_define_method(cgsl_vector, "stats_sd_with_fixed_mean", 
		   rb_gsl_stats_sd_with_fixed_mean, -1);
  rb_define_alias(cgsl_vector, "sd_with_fixed_mean", 
		  "stats_sd_with_fixed_mean");

  rb_define_singleton_method(mgsl_stats, "absdev", rb_gsl_stats_absdev_m, -1);
  rb_define_singleton_method(mgsl_stats, "absdev_m", rb_gsl_stats_absdev_m, -1);
  rb_define_method(cgsl_vector, "stats_absdev_m", rb_gsl_stats_absdev_m, -1);
  rb_define_alias(cgsl_vector, "absdev_m", "stats_absdev_m");
  rb_define_alias(cgsl_vector, "absdev", "stats_absdev_m");

  rb_define_singleton_method(mgsl_stats, "skew", rb_gsl_stats_skew, -1);
  rb_define_singleton_method(mgsl_stats, "skew_m", rb_gsl_stats_skew, -1);
  rb_define_method(cgsl_vector, "stats_skew_m", rb_gsl_stats_skew, -1);
  rb_define_alias(cgsl_vector, "skew_m", "stats_skew_m");
  rb_define_alias(cgsl_vector, "skew", "stats_skew_m");

  rb_define_singleton_method(mgsl_stats, "kurtosis", rb_gsl_stats_kurtosis, -1);
  rb_define_singleton_method(mgsl_stats, "kurtosis_m", rb_gsl_stats_kurtosis, -1);
  rb_define_method(cgsl_vector, "stats_kurtosis_m", rb_gsl_stats_kurtosis, -1);
  rb_define_alias(cgsl_vector, "kurtosis_m", "stats_kurtosis_m");
  rb_define_alias(cgsl_vector, "kurtosis", "stats_kurtosis_m");

  rb_define_singleton_method(mgsl_stats, "lag1_autocorrelation", rb_gsl_stats_lag1_autocorrelation, -1);
  rb_define_singleton_method(mgsl_stats, "lag1_autocorrelation_m", rb_gsl_stats_lag1_autocorrelation, -1);
  rb_define_method(cgsl_vector, "stats_lag1_autocorrelation_m", rb_gsl_stats_lag1_autocorrelation, -1);
  rb_define_alias(cgsl_vector, "lag1_autocorrelation_m", "stats_lag1_autocorrelation_m");
  rb_define_alias(cgsl_vector, "lag1_autocorrelation", "stats_lag1_autocorrelation_m");

  rb_define_singleton_method(mgsl_stats, "covariance", rb_gsl_stats_covariance2, 2);
  rb_define_singleton_method(mgsl_stats, "covariance_m", rb_gsl_stats_covariance_m2, 4);
  
#ifdef GSL_1_10_LATER
  rb_define_singleton_method(mgsl_stats, "correlation", rb_gsl_stats_correlation, 2);  
  rb_define_singleton_method(mgsl_stats, "pvariance", rb_gsl_stats_pvariance, 2);    
  rb_define_singleton_method(mgsl_stats, "ttest", rb_gsl_stats_ttest, 2);      
#endif

  /*****/
  
  rb_define_singleton_method(mgsl_stats, "wmean", rb_gsl_stats_wmean2, -1);
  rb_define_singleton_method(mgsl_stats, "wvariance", rb_gsl_stats_wvariance2, -1);
  rb_define_singleton_method(mgsl_stats, "wvariance_m", rb_gsl_stats_wvariance_m2, -1);
  rb_define_singleton_method(mgsl_stats, "wsd", rb_gsl_stats_wsd2, -1);
  rb_define_singleton_method(mgsl_stats, "wsd_m", rb_gsl_stats_wsd_m2, -1);
  rb_define_singleton_method(mgsl_stats, "wvariance_with_fixed_mean", 
			     rb_gsl_stats_wvariance_with_fixed_mean2, -1);
  rb_define_singleton_method(mgsl_stats, "wsd_with_fixed_mean", 
			     rb_gsl_stats_wsd_with_fixed_mean2, -1);
  rb_define_singleton_method(mgsl_stats, "wabsdev", rb_gsl_stats_wabsdev2, -1);
  rb_define_singleton_method(mgsl_stats, "wabsdev_m", rb_gsl_stats_wabsdev_m2, -1);
  rb_define_singleton_method(mgsl_stats, "wskew", rb_gsl_stats_wskew2, -1);
  rb_define_singleton_method(mgsl_stats, "wskew_m_sd", rb_gsl_stats_wskew_m2, 4);
  rb_define_singleton_method(mgsl_stats, "wkurtosis", rb_gsl_stats_wkurtosis2, -1);
  rb_define_singleton_method(mgsl_stats, "wkurtosis_m_sd", rb_gsl_stats_wkurtosis_m2, 4);

  /*****/

  rb_define_method(cgsl_vector, "stats_wmean", rb_gsl_stats_wmean, -1);
  rb_define_alias(cgsl_vector, "wmean", "stats_wmean");
  rb_define_method(cgsl_vector, "stats_wvariance", rb_gsl_stats_wvariance, -1);
  rb_define_alias(cgsl_vector, "wvariance", "stats_wvariance");
  rb_define_method(cgsl_vector, "stats_wvariance_m", rb_gsl_stats_wvariance_m, -1);
  rb_define_alias(cgsl_vector, "wvariance_m", "stats_wvariance_m");
  rb_define_method(cgsl_vector, "stats_wsd", rb_gsl_stats_wsd, -1);
  rb_define_alias(cgsl_vector, "wsd", "stats_wsd");
  rb_define_method(cgsl_vector, "stats_wsd_m", rb_gsl_stats_wsd_m, -1);
  rb_define_alias(cgsl_vector, "wsd_m", "stats_wsd_m");
  rb_define_method(cgsl_vector, "stats_wvariance_with_fixed_mean", 
		   rb_gsl_stats_wvariance_with_fixed_mean, -1);
  rb_define_alias(cgsl_vector, "wvariance_with_fixed_mean", 
		  "stats_wvariance_with_fixed_mean");
  rb_define_method(cgsl_vector, "stats_wsd_with_fixed_mean", 
		   rb_gsl_stats_wsd_with_fixed_mean, -1);
  rb_define_alias(cgsl_vector, "wsd_with_fixed_mean", "stats_wsd_with_fixed_mean");
  rb_define_method(cgsl_vector, "stats_wabsdev", rb_gsl_stats_wabsdev, -1);
  rb_define_alias(cgsl_vector, "wabsdev", "stats_wabsdev");
  rb_define_method(cgsl_vector, "stats_wabsdev_m", rb_gsl_stats_wabsdev_m, -1);
  rb_define_alias(cgsl_vector, "wabsdev_m", "stats_wabsdev_m");
  rb_define_method(cgsl_vector, "stats_wskew", rb_gsl_stats_wskew, -1);
  rb_define_alias(cgsl_vector, "wskew", "stats_wskew");
  rb_define_method(cgsl_vector, "stats_wskew_m_sd", rb_gsl_stats_wskew_m_sd, 2);
  rb_define_alias(cgsl_vector, "wskew_m_sd", "stats_wskew_m_sd");
  rb_define_method(cgsl_vector, "stats_wkurtosis", rb_gsl_stats_wkurtosis, -1);
  rb_define_alias(cgsl_vector, "wkurtosis", "stats_wkurtosis");
  rb_define_method(cgsl_vector, "stats_wkurtosis_m_sd", 
		   rb_gsl_stats_wkurtosis_m_sd, 2);
  rb_define_alias(cgsl_vector, "wkurtosis_m_sd", "stats_wkurtosis_m_sd");

  /*****/
  rb_define_singleton_method(mgsl_stats, "max", rb_gsl_stats_max, -1);
  rb_define_singleton_method(mgsl_stats, "min", rb_gsl_stats_min, -1);
  rb_define_singleton_method(mgsl_stats, "minmax", rb_gsl_stats_minmax, -1);
  rb_define_singleton_method(mgsl_stats, "max_index", rb_gsl_stats_max_index, -1);
  rb_define_singleton_method(mgsl_stats, "min_index", rb_gsl_stats_min_index, -1);
  rb_define_singleton_method(mgsl_stats, "minmax_index", 
			     rb_gsl_stats_minmax_index, -1);

  rb_define_method(cgsl_vector, "stats_max", rb_gsl_stats_max, -1);
  rb_define_method(cgsl_vector, "stats_min", rb_gsl_stats_min, -1);
  rb_define_method(cgsl_vector, "stats_minmax", rb_gsl_stats_minmax, -1);
  rb_define_method(cgsl_vector, "stats_max_index", rb_gsl_stats_max_index, -1);
  rb_define_method(cgsl_vector, "stats_min_index", rb_gsl_stats_min_index, -1);
  rb_define_method(cgsl_vector, "stats_minmax_index", rb_gsl_stats_minmax_index, -1);

  rb_define_singleton_method(mgsl_stats, "median_from_sorted_data", 
			     rb_gsl_stats_median_from_sorted_data, -1);
  rb_define_method(cgsl_vector, "stats_median_from_sorted_data", 
		   rb_gsl_stats_median_from_sorted_data, -1);
  rb_define_alias(cgsl_vector, "median_from_sorted_data", 
		  "stats_median_from_sorted_data");
  rb_define_method(cgsl_vector, "median", rb_gsl_stats_median, -1);

  rb_define_method(cgsl_vector, "stats_quantile_from_sorted_data", 
		   rb_gsl_stats_quantile_from_sorted_data, -1);
  rb_define_alias(cgsl_vector, "quantile_from_sorted_data", 
		  "stats_quantile_from_sorted_data");

}
