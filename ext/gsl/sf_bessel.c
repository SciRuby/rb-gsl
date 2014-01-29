/*
  sf_bessel.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "include/rb_gsl_sf.h"
EXTERN VALUE cgsl_vector;

/* Cylindrical Bessel Functions */
static VALUE rb_gsl_sf_bessel_J0(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_J0, x);
}

static VALUE rb_gsl_sf_bessel_J0_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_J0_e, x);
}

static VALUE rb_gsl_sf_bessel_J1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_J1, x);
}

static VALUE rb_gsl_sf_bessel_J1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_J1_e, x);
}

static VALUE rb_gsl_sf_bessel_Jn(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_bessel_Jn, n, x);
}

static VALUE rb_gsl_sf_bessel_Jn_e(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_bessel_Jn_e, n, x);
}

static VALUE rb_gsl_sf_bessel_Xn_array(VALUE obj, VALUE n0, VALUE n1, VALUE x,
				       int (*f)(int, int, double, double[]))
{
  int nmin, nmax, n;
  gsl_vector *v = NULL;
  CHECK_FIXNUM(n0); CHECK_FIXNUM(n1);
  Need_Float(x);
  nmin = FIX2INT(n0);
  nmax = FIX2INT(n1);
  n = nmax - nmin + 1;
  v = gsl_vector_alloc(n);
  (*f)(nmin, nmax, NUM2DBL(x), v->data);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

static VALUE rb_gsl_sf_bessel_Jn_array(VALUE obj, VALUE n0, VALUE n1, VALUE x)
{
  return rb_gsl_sf_bessel_Xn_array(obj, n0, n1, x, gsl_sf_bessel_Jn_array);
}

/* Irregular Cylindrical Bessel Functions */
static VALUE rb_gsl_sf_bessel_Y0(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_Y0, x);
}

static VALUE rb_gsl_sf_bessel_Y0_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_Y0_e, x);
}

static VALUE rb_gsl_sf_bessel_Y1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_Y1, x);
}

static VALUE rb_gsl_sf_bessel_Y1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_Y1_e, x);
}

static VALUE rb_gsl_sf_bessel_Yn(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_bessel_Yn, n, x);
}

static VALUE rb_gsl_sf_bessel_Yn_e(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_bessel_Yn_e, n, x);
}

static VALUE rb_gsl_sf_bessel_Yn_array(VALUE obj, VALUE n0, VALUE n1, VALUE x)
{
  return rb_gsl_sf_bessel_Xn_array(obj, n0, n1, x, gsl_sf_bessel_Yn_array);
}

/* Regular Modified Cylindrical Bessel Functions */
static VALUE rb_gsl_sf_bessel_I0(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_I0, x);
}

static VALUE rb_gsl_sf_bessel_I0_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_I0_e, x);
}

static VALUE rb_gsl_sf_bessel_I1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_I1, x);
}

static VALUE rb_gsl_sf_bessel_I1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_I1_e, x);
}

static VALUE rb_gsl_sf_bessel_In(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_bessel_In, n, x);
}

static VALUE rb_gsl_sf_bessel_In_e(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_bessel_In_e, n, x);
}

static VALUE rb_gsl_sf_bessel_In_array(VALUE obj, VALUE n0, VALUE n1, VALUE x)
{
  return rb_gsl_sf_bessel_Xn_array(obj, n0, n1, x, gsl_sf_bessel_In_array);
}

static VALUE rb_gsl_sf_bessel_I0_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_I0_scaled, x);
}

static VALUE rb_gsl_sf_bessel_I0_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_I0_scaled_e, x);
}

static VALUE rb_gsl_sf_bessel_I1_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_I1_scaled, x);
}

static VALUE rb_gsl_sf_bessel_I1_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_I1_scaled_e, x);
}

static VALUE rb_gsl_sf_bessel_In_scaled(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_bessel_In_scaled, n, x);
}

static VALUE rb_gsl_sf_bessel_In_scaled_e(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_bessel_In_scaled_e, n, x);
}

static VALUE rb_gsl_sf_bessel_In_scaled_array(VALUE obj, VALUE n0, VALUE n1,
					      VALUE x)
{
  return rb_gsl_sf_bessel_Xn_array(obj, n0, n1, x, gsl_sf_bessel_In_scaled_array);
}

/* Irregular Modified Cylindrical Bessel Functions */
static VALUE rb_gsl_sf_bessel_K0(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_K0, x);
}

static VALUE rb_gsl_sf_bessel_K0_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_K0_e, x);
}

static VALUE rb_gsl_sf_bessel_K1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_K1, x);
}

static VALUE rb_gsl_sf_bessel_K1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_K1_e, x);
}

static VALUE rb_gsl_sf_bessel_Kn(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_bessel_Kn, n, x);
}

static VALUE rb_gsl_sf_bessel_Kn_e(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_bessel_Kn_e, n, x);
}

static VALUE rb_gsl_sf_bessel_Kn_array(VALUE obj, VALUE n0, VALUE n1, VALUE x)
{
  return rb_gsl_sf_bessel_Xn_array(obj, n0, n1, x, gsl_sf_bessel_Kn_array);
}

static VALUE rb_gsl_sf_bessel_K0_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_K0_scaled, x);
}

static VALUE rb_gsl_sf_bessel_K0_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_K0_scaled_e, x);
}

static VALUE rb_gsl_sf_bessel_K1_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_K1_scaled, x);
}

static VALUE rb_gsl_sf_bessel_K1_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_K1_scaled_e, x);
}

static VALUE rb_gsl_sf_bessel_Kn_scaled(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_bessel_Kn_scaled, n, x);
}

static VALUE rb_gsl_sf_bessel_Kn_scaled_e(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_bessel_Kn_scaled_e, n, x);
}

static VALUE rb_gsl_sf_bessel_Kn_scaled_array(VALUE obj, VALUE n0, VALUE n1,
					      VALUE x)
{
  return rb_gsl_sf_bessel_Xn_array(obj, n0, n1, x, gsl_sf_bessel_Kn_scaled_array);
}

/* Spherical Bessel Functions */
static VALUE rb_gsl_sf_bessel_j0(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_j0, x);
}

static VALUE rb_gsl_sf_bessel_j0_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_j0_e, x);
}

static VALUE rb_gsl_sf_bessel_j1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_j1, x);
}

static VALUE rb_gsl_sf_bessel_j1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_j1_e, x);
}

static VALUE rb_gsl_sf_bessel_j2(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_j2, x);
}

static VALUE rb_gsl_sf_bessel_j2_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_j2_e, x);
}

static VALUE rb_gsl_sf_bessel_jl(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_bessel_jl, n, x);
}

static VALUE rb_gsl_sf_bessel_jl_e(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_bessel_jl_e, n, x);
}

static VALUE rb_gsl_sf_bessel_xl_array(VALUE obj, VALUE n1, VALUE x,
				       int (*f)(int, double, double[]))
{
  int nmax, n;
  // local variable "status" declared and set, but never used
  //int status;
  gsl_vector *v = NULL;
  CHECK_FIXNUM(n1);
  Need_Float(x);
  nmax = FIX2INT(n1);
  n = nmax  + 1;
  v = gsl_vector_alloc(n);
  /*status =*/ (*f)(nmax, NUM2DBL(x), v->data);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

static VALUE rb_gsl_sf_bessel_jl_array(VALUE obj, VALUE n1, VALUE x)
{
  return rb_gsl_sf_bessel_xl_array(obj, n1, x, gsl_sf_bessel_jl_array);
}

static VALUE rb_gsl_sf_bessel_jl_steed_array(VALUE obj, VALUE n1, VALUE x)
{
  return rb_gsl_sf_bessel_xl_array(obj, n1, x, gsl_sf_bessel_jl_steed_array);
}

/* Irregular Cylindrical Bessel Functions */
static VALUE rb_gsl_sf_bessel_y0(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_y0, x);
}

static VALUE rb_gsl_sf_bessel_y0_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_y0_e, x);
}

static VALUE rb_gsl_sf_bessel_y1(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_y1, x);
}

static VALUE rb_gsl_sf_bessel_y1_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_y1_e, x);
}

static VALUE rb_gsl_sf_bessel_y2(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_y2, x);
}

static VALUE rb_gsl_sf_bessel_y2_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_y2_e, x);
}

static VALUE rb_gsl_sf_bessel_yl(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_bessel_yl, n, x);
}

static VALUE rb_gsl_sf_bessel_yl_e(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_bessel_yl_e, n, x);
}

static VALUE rb_gsl_sf_bessel_yl_array(VALUE obj, VALUE n1, VALUE x)
{
  return rb_gsl_sf_bessel_xl_array(obj, n1, x, gsl_sf_bessel_yl_array);
}

/* Regular Modified Cylindrical Bessel Functions */
static VALUE rb_gsl_sf_bessel_i0_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_i0_scaled, x);
}

static VALUE rb_gsl_sf_bessel_i0_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_i0_scaled_e, x);
}

static VALUE rb_gsl_sf_bessel_i1_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_i1_scaled, x);
}

static VALUE rb_gsl_sf_bessel_i1_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_i1_scaled_e, x);
}

static VALUE rb_gsl_sf_bessel_i2_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_i2_scaled, x);
}

static VALUE rb_gsl_sf_bessel_i2_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_i2_scaled_e, x);
}

static VALUE rb_gsl_sf_bessel_il_scaled(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_bessel_il_scaled, n, x);
}

static VALUE rb_gsl_sf_bessel_il_scaled_e(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_bessel_il_scaled_e, n, x);
}

static VALUE rb_gsl_sf_bessel_il_scaled_array(VALUE obj, VALUE n1, VALUE x)
{
  return rb_gsl_sf_bessel_xl_array(obj, n1, x, gsl_sf_bessel_il_scaled_array);
}

/* Irregular Modified Cylindrical Bessel Functions */

static VALUE rb_gsl_sf_bessel_k0_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_k0_scaled, x);
}

static VALUE rb_gsl_sf_bessel_k0_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_k0_scaled_e, x);
}

static VALUE rb_gsl_sf_bessel_k1_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_k1_scaled, x);
}

static VALUE rb_gsl_sf_bessel_k1_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_k1_scaled_e, x);
}

static VALUE rb_gsl_sf_bessel_k2_scaled(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval1(gsl_sf_bessel_k2_scaled, x);
}

static VALUE rb_gsl_sf_bessel_k2_scaled_e(VALUE obj, VALUE x)
{
  return rb_gsl_sf_eval_e(gsl_sf_bessel_k2_scaled_e, x);
}

static VALUE rb_gsl_sf_bessel_kl_scaled(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_int_double(gsl_sf_bessel_kl_scaled, n, x);
}

static VALUE rb_gsl_sf_bessel_kl_scaled_e(VALUE obj, VALUE n, VALUE x)
{
  return rb_gsl_sf_eval_e_int_double(gsl_sf_bessel_kl_scaled_e, n, x);
}

static VALUE rb_gsl_sf_bessel_kl_scaled_array(VALUE obj, VALUE n1,
					      VALUE x)
{
  return rb_gsl_sf_bessel_xl_array(obj, n1, x, gsl_sf_bessel_kl_scaled_array);
}

/* Regular Bessel Function - Fractional Order */

static VALUE rb_gsl_sf_bessel_Jnu(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_bessel_Jnu, nu, x);
}

static VALUE rb_gsl_sf_bessel_Jnu_e(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_bessel_Jnu_e, nu, x);
}

/* The new array will be returned, the original array "ary" is not modified */
static VALUE rb_gsl_sf_bessel_sequence_Jnu_e(int argc, VALUE *argv, VALUE obj)
{
  size_t i, size;
  // local variable "status" declared and set, but never used
  //int status;
  int flag = 0;
  gsl_vector *v = NULL;
  gsl_mode_t mode;
  char c;
  VALUE nu, m, ary;

  nu = argv[0];
  switch (argc) {
  case 2:
    ary = argv[1];
    mode = GSL_PREC_DOUBLE;
    break;
  case 3:
    m = argv[1];
    ary = argv[2];
    switch (TYPE(m)) {
    case T_STRING:
      c = tolower(NUM2CHR(m));
      if (c == 'd') mode = GSL_PREC_DOUBLE;
      else if (c == 's') mode = GSL_PREC_SINGLE;
      else if (c == 'a') mode = GSL_PREC_APPROX;
      else mode = GSL_PREC_DOUBLE;
      break;
    case T_FIXNUM:
      mode = FIX2INT(m);
      break;
    default:
      rb_raise(rb_eTypeError, "wrong argument type %s (String or Fixnum expected)",
	       rb_class2name(CLASS_OF(m)));
      break;
    }
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }

  switch (TYPE(ary)) {
  case T_ARRAY:
    size = RARRAY_LEN(ary);
    v = gsl_vector_alloc(size);
    for (i = 0; i < size; i++) gsl_vector_set(v, i, NUM2DBL(rb_ary_entry(ary, i)));
    flag = 1;
    break;
  default:
    CHECK_VECTOR(ary);
    Data_Get_Struct(ary, gsl_vector, v);
    size = v->size;
    flag = 0;
    break;
  }
  /*status =*/ gsl_sf_bessel_sequence_Jnu_e(NUM2DBL(nu), mode, size, v->data);
  if (flag == 1) return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
  else return ary;
}

/* Irregular Bessel Function - Fractional Order */
static VALUE rb_gsl_sf_bessel_Ynu(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_bessel_Ynu, nu, x);
}

static VALUE rb_gsl_sf_bessel_Ynu_e(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_bessel_Ynu_e, nu, x);
}

/* Regular Modified Bessel Function - Fractional Order */
static VALUE rb_gsl_sf_bessel_Inu(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_bessel_Inu, nu, x);
}

static VALUE rb_gsl_sf_bessel_Inu_e(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_bessel_Inu_e, nu, x);
}

static VALUE rb_gsl_sf_bessel_Inu_scaled(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_bessel_Inu_scaled, nu, x);
}

static VALUE rb_gsl_sf_bessel_Inu_scaled_e(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_bessel_Inu_scaled_e, nu, x);
}

/* Irregular Modified Bessel Function - Fractional Order */
static VALUE rb_gsl_sf_bessel_Knu(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_bessel_Knu, nu, x);
}

static VALUE rb_gsl_sf_bessel_Knu_e(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_bessel_Knu_e, nu, x);
}

static VALUE rb_gsl_sf_bessel_lnKnu(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_bessel_lnKnu, nu, x);
}

static VALUE rb_gsl_sf_bessel_lnKnu_e(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_bessel_lnKnu_e, nu, x);
}

static VALUE rb_gsl_sf_bessel_Knu_scaled(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_double_double(gsl_sf_bessel_Knu_scaled, nu, x);
}

static VALUE rb_gsl_sf_bessel_Knu_scaled_e(VALUE obj, VALUE nu, VALUE x)
{
  return rb_gsl_sf_eval_e_double2(gsl_sf_bessel_Knu_scaled_e, nu, x);
}

static VALUE rb_gsl_sf_eval_double_uint(double (*func)(double, unsigned int), VALUE ff, VALUE argv);

static VALUE rb_gsl_sf_eval_double_uint(double (*func)(double, unsigned int), VALUE ff, VALUE argv)
{
  gsl_vector *v, *vnew;
  VALUE ary;
  size_t i, n;
  double val, f;

  f = NUM2DBL(ff);
  switch (TYPE(argv)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new((*func)(f, NUM2UINT(argv)));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      val = (*func)(f, NUM2UINT(rb_ary_entry(argv, i)));
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
    CHECK_VECTOR(argv);
    Data_Get_Struct(argv, gsl_vector, v);
    n = v->size;
    vnew = gsl_vector_alloc(n);
    for (i = 0; i < n; i++) {
      val = (*func)(f, (unsigned int) gsl_vector_get(v, i));
      gsl_vector_set(vnew, i, val);
    }
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    break;
  }
}

/* find s-th 0 point of J0 */
static VALUE rb_gsl_sf_bessel_zero_J0(VALUE obj, VALUE s)
{
  return rb_gsl_sf_eval1_uint(gsl_sf_bessel_zero_J0, s);
}

static VALUE rb_gsl_sf_bessel_zero_J0_e(VALUE obj, VALUE s)
{
  return rb_gsl_sf_eval_e_uint(gsl_sf_bessel_zero_J0_e, s);
}

static VALUE rb_gsl_sf_bessel_zero_J1(VALUE obj, VALUE s)
{
  return rb_gsl_sf_eval1_uint(gsl_sf_bessel_zero_J1, s);
}

static VALUE rb_gsl_sf_bessel_zero_J1_e(VALUE obj, VALUE s)
{
  return rb_gsl_sf_eval_e_uint(gsl_sf_bessel_zero_J1_e, s);
}

static VALUE rb_gsl_sf_bessel_zero_Jnu(VALUE obj, VALUE n, VALUE s)
{
  return rb_gsl_sf_eval_double_uint(gsl_sf_bessel_zero_Jnu, n, s);
}

static VALUE rb_gsl_sf_bessel_zero_Jnu_e(VALUE obj, VALUE n, VALUE s)
{
  return rb_gsl_sf_eval_e_double_uint(gsl_sf_bessel_zero_Jnu_e, n, s);
}

void Init_gsl_sf_bessel(VALUE module)
{
  VALUE mgsl_sf_bessel;

  rb_define_module_function(module, "bessel_J0",  rb_gsl_sf_bessel_J0, 1);
  rb_define_module_function(module, "bessel_J0_e",  rb_gsl_sf_bessel_J0_e, 1);
  rb_define_module_function(module, "bessel_J1",  rb_gsl_sf_bessel_J1, 1);
  rb_define_module_function(module, "bessel_J1_e",  rb_gsl_sf_bessel_J1_e, 1);
  rb_define_module_function(module, "bessel_Jn",  rb_gsl_sf_bessel_Jn, 2);
  rb_define_module_function(module, "bessel_Jn_e",  rb_gsl_sf_bessel_Jn_e, 2);
  rb_define_module_function(module, "bessel_Jn_array",  rb_gsl_sf_bessel_Jn_array, 3);
  rb_define_module_function(module, "bessel_Y0",  rb_gsl_sf_bessel_Y0, 1);
  rb_define_module_function(module, "bessel_Y0_e",  rb_gsl_sf_bessel_Y0_e, 1);
  rb_define_module_function(module, "bessel_Y1",  rb_gsl_sf_bessel_Y1, 1);
  rb_define_module_function(module, "bessel_Y1_e",  rb_gsl_sf_bessel_Y1_e, 1);
  rb_define_module_function(module, "bessel_Yn",  rb_gsl_sf_bessel_Yn, 2);
  rb_define_module_function(module, "bessel_Yn_e",  rb_gsl_sf_bessel_Yn_e, 2);
  rb_define_module_function(module, "bessel_Yn_array",  rb_gsl_sf_bessel_Yn_array, 3);
  rb_define_module_function(module, "bessel_I0",  rb_gsl_sf_bessel_I0, 1);
  rb_define_module_function(module, "bessel_I0_e",  rb_gsl_sf_bessel_I0_e, 1);
  rb_define_module_function(module, "bessel_I1",  rb_gsl_sf_bessel_I1, 1);
  rb_define_module_function(module, "bessel_I1_e",  rb_gsl_sf_bessel_I1_e, 1);
  rb_define_module_function(module, "bessel_In",  rb_gsl_sf_bessel_In, 2);
  rb_define_module_function(module, "bessel_In_e",  rb_gsl_sf_bessel_In_e, 2);
  rb_define_module_function(module, "bessel_In_array",  rb_gsl_sf_bessel_In_array, 3);
  rb_define_module_function(module, "bessel_I0_scaled",  rb_gsl_sf_bessel_I0_scaled, 1);
  rb_define_module_function(module, "bessel_I0_scaled_e",  rb_gsl_sf_bessel_I0_scaled_e, 1);
  rb_define_module_function(module, "bessel_I1_scaled",  rb_gsl_sf_bessel_I1_scaled, 1);
  rb_define_module_function(module, "bessel_I1_scaled_e",  rb_gsl_sf_bessel_I1_scaled_e, 1);
  rb_define_module_function(module, "bessel_In_scaled",  rb_gsl_sf_bessel_In_scaled, 2);
  rb_define_module_function(module, "bessel_In_scaled_e",  rb_gsl_sf_bessel_In_scaled_e, 2);
  rb_define_module_function(module, "bessel_In_scaled_array",  rb_gsl_sf_bessel_In_scaled_array, 3);
  rb_define_module_function(module, "bessel_K0",  rb_gsl_sf_bessel_K0, 1);
  rb_define_module_function(module, "bessel_K0_e",  rb_gsl_sf_bessel_K0_e, 1);
  rb_define_module_function(module, "bessel_K1",  rb_gsl_sf_bessel_K1, 1);
  rb_define_module_function(module, "bessel_K1_e",  rb_gsl_sf_bessel_K1_e, 1);
  rb_define_module_function(module, "bessel_Kn",  rb_gsl_sf_bessel_Kn, 2);
  rb_define_module_function(module, "bessel_Kn_e",  rb_gsl_sf_bessel_Kn_e, 2);
  rb_define_module_function(module, "bessel_Kn_array",  rb_gsl_sf_bessel_Kn_array, 3);
  rb_define_module_function(module, "bessel_K0_scaled",  rb_gsl_sf_bessel_K0_scaled, 1);
  rb_define_module_function(module, "bessel_K0_scaled_e",  rb_gsl_sf_bessel_K0_scaled_e, 1);
  rb_define_module_function(module, "bessel_K1_scaled",  rb_gsl_sf_bessel_K1_scaled, 1);
  rb_define_module_function(module, "bessel_K1_scaled_e",  rb_gsl_sf_bessel_K1_scaled_e, 1);
  rb_define_module_function(module, "bessel_Kn_scaled",  rb_gsl_sf_bessel_Kn_scaled, 2);
  rb_define_module_function(module, "bessel_Kn_scaled_e",  rb_gsl_sf_bessel_Kn_scaled_e, 2);
  rb_define_module_function(module, "bessel_Kn_scaled_array",  rb_gsl_sf_bessel_Kn_scaled_array, 3);
  rb_define_module_function(module, "bessel_j0",  rb_gsl_sf_bessel_j0, 1);
  rb_define_module_function(module, "bessel_j0_e",  rb_gsl_sf_bessel_j0_e, 1);
  rb_define_module_function(module, "bessel_j1",  rb_gsl_sf_bessel_j1, 1);
  rb_define_module_function(module, "bessel_j1_e",  rb_gsl_sf_bessel_j1_e, 1);
  rb_define_module_function(module, "bessel_j2",  rb_gsl_sf_bessel_j2, 1);
  rb_define_module_function(module, "bessel_j2_e",  rb_gsl_sf_bessel_j2_e, 1);
  rb_define_module_function(module, "bessel_jl",  rb_gsl_sf_bessel_jl, 2);
  rb_define_module_function(module, "bessel_jl_e",  rb_gsl_sf_bessel_jl_e, 2);
  rb_define_module_function(module, "bessel_jl_array",  rb_gsl_sf_bessel_jl_array, 2);
  rb_define_module_function(module, "bessel_jl_steed_array",  rb_gsl_sf_bessel_jl_steed_array, 2);
  rb_define_module_function(module, "bessel_y0",  rb_gsl_sf_bessel_y0, 1);
  rb_define_module_function(module, "bessel_y0_e",  rb_gsl_sf_bessel_y0_e, 1);
  rb_define_module_function(module, "bessel_y1",  rb_gsl_sf_bessel_y1, 1);
  rb_define_module_function(module, "bessel_y1_e",  rb_gsl_sf_bessel_y1_e, 1);
  rb_define_module_function(module, "bessel_y2",  rb_gsl_sf_bessel_y2, 1);
  rb_define_module_function(module, "bessel_y2_e",  rb_gsl_sf_bessel_y2_e, 1);
  rb_define_module_function(module, "bessel_yl",  rb_gsl_sf_bessel_yl, 2);
  rb_define_module_function(module, "bessel_yl_e",  rb_gsl_sf_bessel_yl_e, 2);
  rb_define_module_function(module, "bessel_yl_array",  rb_gsl_sf_bessel_yl_array, 2);
  rb_define_module_function(module, "bessel_i0_scaled",  rb_gsl_sf_bessel_i0_scaled, 1);
  rb_define_module_function(module, "bessel_i0_scaled_e",  rb_gsl_sf_bessel_i0_scaled_e, 1);
  rb_define_module_function(module, "bessel_i1_scaled",  rb_gsl_sf_bessel_i1_scaled, 1);
  rb_define_module_function(module, "bessel_i1_scaled_e",  rb_gsl_sf_bessel_i1_scaled_e, 1);
  rb_define_module_function(module, "bessel_i2_scaled",  rb_gsl_sf_bessel_i2_scaled, 1);
  rb_define_module_function(module, "bessel_i2_scaled_e",  rb_gsl_sf_bessel_i2_scaled_e, 1);
  rb_define_module_function(module, "bessel_il_scaled",  rb_gsl_sf_bessel_il_scaled, 2);
  rb_define_module_function(module, "bessel_il_scaled_e",  rb_gsl_sf_bessel_il_scaled_e, 2);
  rb_define_module_function(module, "bessel_il_scaled_array",  rb_gsl_sf_bessel_il_scaled_array, 2);
  rb_define_module_function(module, "bessel_k0_scaled",  rb_gsl_sf_bessel_k0_scaled, 1);
  rb_define_module_function(module, "bessel_k0_scaled_e",  rb_gsl_sf_bessel_k0_scaled_e, 1);
  rb_define_module_function(module, "bessel_k1_scaled",  rb_gsl_sf_bessel_k1_scaled, 1);
  rb_define_module_function(module, "bessel_k1_scaled_e",  rb_gsl_sf_bessel_k1_scaled_e, 1);
  rb_define_module_function(module, "bessel_k2_scaled",  rb_gsl_sf_bessel_k2_scaled, 1);
  rb_define_module_function(module, "bessel_k2_scaled_e",  rb_gsl_sf_bessel_k2_scaled_e, 1);
  rb_define_module_function(module, "bessel_kl_scaled",  rb_gsl_sf_bessel_kl_scaled, 2);
  rb_define_module_function(module, "bessel_kl_scaled_e",  rb_gsl_sf_bessel_kl_scaled_e, 2);
  rb_define_module_function(module, "bessel_kl_scaled_array",  rb_gsl_sf_bessel_kl_scaled_array, 2);
  rb_define_module_function(module, "bessel_Jnu",  rb_gsl_sf_bessel_Jnu, 2);
  rb_define_module_function(module, "bessel_Jnu_e",  rb_gsl_sf_bessel_Jnu_e, 2);
  rb_define_module_function(module, "bessel_sequence_Jnu_e",  rb_gsl_sf_bessel_sequence_Jnu_e, -1);

  rb_define_module_function(module, "bessel_Ynu",  rb_gsl_sf_bessel_Ynu, 2);
  rb_define_module_function(module, "bessel_Ynu_e",  rb_gsl_sf_bessel_Ynu_e, 2);
  rb_define_module_function(module, "bessel_Inu",  rb_gsl_sf_bessel_Inu, 2);
  rb_define_module_function(module, "bessel_Inu_e",  rb_gsl_sf_bessel_Inu_e, 2);
  rb_define_module_function(module, "bessel_Inu_scaled",  rb_gsl_sf_bessel_Inu_scaled, 2);
  rb_define_module_function(module, "bessel_Inu_scaled_e",  rb_gsl_sf_bessel_Inu_scaled_e, 2);
  rb_define_module_function(module, "bessel_Knu",  rb_gsl_sf_bessel_Knu, 2);
  rb_define_module_function(module, "bessel_Knu_e",  rb_gsl_sf_bessel_Knu_e, 2);
  rb_define_module_function(module, "bessel_lnKnu",  rb_gsl_sf_bessel_lnKnu, 2);
  rb_define_module_function(module, "bessel_lnKnu_e",  rb_gsl_sf_bessel_lnKnu_e, 2);
  rb_define_module_function(module, "bessel_Knu_scaled",  rb_gsl_sf_bessel_Knu_scaled, 2);
  rb_define_module_function(module, "bessel_Knu_scaled_e",  rb_gsl_sf_bessel_Knu_scaled_e, 2);
  rb_define_module_function(module, "bessel_zero_J0",  rb_gsl_sf_bessel_zero_J0, 1);
  rb_define_module_function(module, "bessel_zero_J0_e",  rb_gsl_sf_bessel_zero_J0_e, 1);
  rb_define_module_function(module, "bessel_zero_J1",  rb_gsl_sf_bessel_zero_J1, 1);
  rb_define_module_function(module, "bessel_zero_J1_e",  rb_gsl_sf_bessel_zero_J1_e, 1);
  rb_define_module_function(module, "bessel_zero_Jnu",  rb_gsl_sf_bessel_zero_Jnu, 2);
  rb_define_module_function(module, "bessel_zero_Jnu_e",  rb_gsl_sf_bessel_zero_Jnu_e, 2);

  /*******************************/
  mgsl_sf_bessel = rb_define_module_under(module, "Bessel");

  rb_define_module_function(mgsl_sf_bessel, "J0",  rb_gsl_sf_bessel_J0, 1);
  rb_define_module_function(mgsl_sf_bessel, "J0_e",  rb_gsl_sf_bessel_J0_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "J1",  rb_gsl_sf_bessel_J1, 1);
  rb_define_module_function(mgsl_sf_bessel, "J1_e",  rb_gsl_sf_bessel_J1_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "Jn",  rb_gsl_sf_bessel_Jn, 2);
  rb_define_module_function(mgsl_sf_bessel, "Jn_e",  rb_gsl_sf_bessel_Jn_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "Jn_array",  rb_gsl_sf_bessel_Jn_array, 3);
  rb_define_module_function(mgsl_sf_bessel, "Y0",  rb_gsl_sf_bessel_Y0, 1);
  rb_define_module_function(mgsl_sf_bessel, "Y0_e",  rb_gsl_sf_bessel_Y0_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "Y1",  rb_gsl_sf_bessel_Y1, 1);
  rb_define_module_function(mgsl_sf_bessel, "Y1_e",  rb_gsl_sf_bessel_Y1_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "Yn",  rb_gsl_sf_bessel_Yn, 2);
  rb_define_module_function(mgsl_sf_bessel, "Yn_e",  rb_gsl_sf_bessel_Yn_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "Yn_array",  rb_gsl_sf_bessel_Yn_array, 3);
  rb_define_module_function(mgsl_sf_bessel, "I0",  rb_gsl_sf_bessel_I0, 1);
  rb_define_module_function(mgsl_sf_bessel, "I0_e",  rb_gsl_sf_bessel_I0_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "I1",  rb_gsl_sf_bessel_I1, 1);
  rb_define_module_function(mgsl_sf_bessel, "I1_e",  rb_gsl_sf_bessel_I1_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "In",  rb_gsl_sf_bessel_In, 2);
  rb_define_module_function(mgsl_sf_bessel, "In_e",  rb_gsl_sf_bessel_In_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "In_array",  rb_gsl_sf_bessel_In_array, 3);
  rb_define_module_function(mgsl_sf_bessel, "I0_scaled",  rb_gsl_sf_bessel_I0_scaled, 1);
  rb_define_module_function(mgsl_sf_bessel, "I0_scaled_e",  rb_gsl_sf_bessel_I0_scaled_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "I1_scaled",  rb_gsl_sf_bessel_I1_scaled, 1);
  rb_define_module_function(mgsl_sf_bessel, "I1_scaled_e",  rb_gsl_sf_bessel_I1_scaled_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "In_scaled",  rb_gsl_sf_bessel_In_scaled, 2);
  rb_define_module_function(mgsl_sf_bessel, "In_scaled_e",  rb_gsl_sf_bessel_In_scaled_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "In_scaled_array",  rb_gsl_sf_bessel_In_scaled_array, 3);
  rb_define_module_function(mgsl_sf_bessel, "K0",  rb_gsl_sf_bessel_K0, 1);
  rb_define_module_function(mgsl_sf_bessel, "K0_e",  rb_gsl_sf_bessel_K0_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "K1",  rb_gsl_sf_bessel_K1, 1);
  rb_define_module_function(mgsl_sf_bessel, "K1_e",  rb_gsl_sf_bessel_K1_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "Kn",  rb_gsl_sf_bessel_Kn, 2);
  rb_define_module_function(mgsl_sf_bessel, "Kn_e",  rb_gsl_sf_bessel_Kn_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "Kn_array",  rb_gsl_sf_bessel_Kn_array, 3);
  rb_define_module_function(mgsl_sf_bessel, "K0_scaled",  rb_gsl_sf_bessel_K0_scaled, 1);
  rb_define_module_function(mgsl_sf_bessel, "K0_scaled_e",  rb_gsl_sf_bessel_K0_scaled_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "K1_scaled",  rb_gsl_sf_bessel_K1_scaled, 1);
  rb_define_module_function(mgsl_sf_bessel, "K1_scaled_e",  rb_gsl_sf_bessel_K1_scaled_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "Kn_scaled",  rb_gsl_sf_bessel_Kn_scaled, 2);
  rb_define_module_function(mgsl_sf_bessel, "Kn_scaled_e",  rb_gsl_sf_bessel_Kn_scaled_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "Kn_scaled_array",  rb_gsl_sf_bessel_Kn_scaled_array, 3);
  rb_define_module_function(mgsl_sf_bessel, "j0",  rb_gsl_sf_bessel_j0, 1);
  rb_define_module_function(mgsl_sf_bessel, "j0_e",  rb_gsl_sf_bessel_j0_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "j1",  rb_gsl_sf_bessel_j1, 1);
  rb_define_module_function(mgsl_sf_bessel, "j1_e",  rb_gsl_sf_bessel_j1_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "j2",  rb_gsl_sf_bessel_j2, 1);
  rb_define_module_function(mgsl_sf_bessel, "j2_e",  rb_gsl_sf_bessel_j2_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "jl",  rb_gsl_sf_bessel_jl, 2);
  rb_define_module_function(mgsl_sf_bessel, "jl_e",  rb_gsl_sf_bessel_jl_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "jl_array",  rb_gsl_sf_bessel_jl_array, 2);
  rb_define_module_function(mgsl_sf_bessel, "jl_steed_array",  rb_gsl_sf_bessel_jl_steed_array, 2);
  rb_define_module_function(mgsl_sf_bessel, "y0",  rb_gsl_sf_bessel_y0, 1);
  rb_define_module_function(mgsl_sf_bessel, "y0_e",  rb_gsl_sf_bessel_y0_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "y1",  rb_gsl_sf_bessel_y1, 1);
  rb_define_module_function(mgsl_sf_bessel, "y1_e",  rb_gsl_sf_bessel_y1_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "y2",  rb_gsl_sf_bessel_y2, 1);
  rb_define_module_function(mgsl_sf_bessel, "y2_e",  rb_gsl_sf_bessel_y2_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "yl",  rb_gsl_sf_bessel_yl, 2);
  rb_define_module_function(mgsl_sf_bessel, "yl_e",  rb_gsl_sf_bessel_yl_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "yl_array",  rb_gsl_sf_bessel_yl_array, 2);
  rb_define_module_function(mgsl_sf_bessel, "i0_scaled",  rb_gsl_sf_bessel_i0_scaled, 1);
  rb_define_module_function(mgsl_sf_bessel, "i0_scaled_e",  rb_gsl_sf_bessel_i0_scaled_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "i1_scaled",  rb_gsl_sf_bessel_i1_scaled, 1);
  rb_define_module_function(mgsl_sf_bessel, "i1_scaled_e",  rb_gsl_sf_bessel_i1_scaled_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "i2_scaled",  rb_gsl_sf_bessel_i2_scaled, 1);
  rb_define_module_function(mgsl_sf_bessel, "i2_scaled_e",  rb_gsl_sf_bessel_i2_scaled_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "il_scaled",  rb_gsl_sf_bessel_il_scaled, 2);
  rb_define_module_function(mgsl_sf_bessel, "il_scaled_e",  rb_gsl_sf_bessel_il_scaled_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "il_scaled_array",  rb_gsl_sf_bessel_il_scaled_array, 2);
  rb_define_module_function(mgsl_sf_bessel, "k0_scaled",  rb_gsl_sf_bessel_k0_scaled, 1);
  rb_define_module_function(mgsl_sf_bessel, "k0_scaled_e",  rb_gsl_sf_bessel_k0_scaled_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "k1_scaled",  rb_gsl_sf_bessel_k1_scaled, 1);
  rb_define_module_function(mgsl_sf_bessel, "k1_scaled_e",  rb_gsl_sf_bessel_k1_scaled_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "k2_scaled",  rb_gsl_sf_bessel_k2_scaled, 1);
  rb_define_module_function(mgsl_sf_bessel, "k2_scaled_e",  rb_gsl_sf_bessel_k2_scaled_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "kl_scaled",  rb_gsl_sf_bessel_kl_scaled, 2);
  rb_define_module_function(mgsl_sf_bessel, "kl_scaled_e",  rb_gsl_sf_bessel_kl_scaled_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "kl_scaled_array",  rb_gsl_sf_bessel_kl_scaled_array, 2);
  rb_define_module_function(mgsl_sf_bessel, "Jnu",  rb_gsl_sf_bessel_Jnu, 2);
  rb_define_module_function(mgsl_sf_bessel, "Jnu_e",  rb_gsl_sf_bessel_Jnu_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "sequence_Jnu_e",  rb_gsl_sf_bessel_sequence_Jnu_e, 3);

  rb_define_module_function(mgsl_sf_bessel, "Ynu",  rb_gsl_sf_bessel_Ynu, 2);
  rb_define_module_function(mgsl_sf_bessel, "Ynu_e",  rb_gsl_sf_bessel_Ynu_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "Inu",  rb_gsl_sf_bessel_Inu, 2);
  rb_define_module_function(mgsl_sf_bessel, "Inu_e",  rb_gsl_sf_bessel_Inu_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "Inu_scaled",  rb_gsl_sf_bessel_Inu_scaled, 2);
  rb_define_module_function(mgsl_sf_bessel, "Inu_scaled_e",  rb_gsl_sf_bessel_Inu_scaled_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "Knu",  rb_gsl_sf_bessel_Knu, 2);
  rb_define_module_function(mgsl_sf_bessel, "Knu_e",  rb_gsl_sf_bessel_Knu_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "lnKnu",  rb_gsl_sf_bessel_lnKnu, 2);
  rb_define_module_function(mgsl_sf_bessel, "lnKnu_e",  rb_gsl_sf_bessel_lnKnu_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "Knu_scaled",  rb_gsl_sf_bessel_Knu_scaled, 2);
  rb_define_module_function(mgsl_sf_bessel, "Knu_scaled_e",  rb_gsl_sf_bessel_Knu_scaled_e, 2);
  rb_define_module_function(mgsl_sf_bessel, "zero_J0",  rb_gsl_sf_bessel_zero_J0, 1);
  rb_define_module_function(mgsl_sf_bessel, "zero_J0_e",  rb_gsl_sf_bessel_zero_J0_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "zero_J1",  rb_gsl_sf_bessel_zero_J1, 1);
  rb_define_module_function(mgsl_sf_bessel, "zero_J1_e",  rb_gsl_sf_bessel_zero_J1_e, 1);
  rb_define_module_function(mgsl_sf_bessel, "zero_Jnu",  rb_gsl_sf_bessel_zero_Jnu, 2);
  rb_define_module_function(mgsl_sf_bessel, "zero_Jnu_e",  rb_gsl_sf_bessel_zero_Jnu_e, 2);

}
