/*
  sf.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "rb_gsl_array.h"
#include "rb_gsl_sf.h"
#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif

VALUE cgsl_sf_result, cgsl_sf_result_e10;

VALUE rb_gsl_sf_result_new(VALUE klass)
{
  gsl_sf_result *rslt = NULL;
  return Data_Make_Struct(klass, gsl_sf_result, 0, free, rslt);
}

static VALUE rb_gsl_sf_result_print(VALUE obj)
{
  gsl_sf_result *rslt = NULL;
  Data_Get_Struct(obj, gsl_sf_result, rslt);
  printf("%10.9e %10.9e\n", rslt->val, rslt->err);
  return obj;
}

static VALUE rb_gsl_sf_result_to_s(VALUE obj);

static VALUE rb_gsl_sf_result_inspect(VALUE obj)
{
  char buf[64];
  VALUE str;
  sprintf(buf, "%s\n", rb_class2name(CLASS_OF(obj)));
  str = rb_str_new2(buf);
  return rb_str_concat(str, rb_gsl_sf_result_to_s(obj));
}

static VALUE rb_gsl_sf_result_val(VALUE obj)
{
  gsl_sf_result *rslt = NULL;
  Data_Get_Struct(obj, gsl_sf_result, rslt);
  return rb_float_new(rslt->val);
}

static VALUE rb_gsl_sf_result_err(VALUE obj)
{
  gsl_sf_result *rslt = NULL;
  Data_Get_Struct(obj, gsl_sf_result, rslt);
  return rb_float_new(rslt->err);
}

static VALUE rb_gsl_sf_result_to_a(VALUE obj)
{
  gsl_sf_result *rslt = NULL;
  Data_Get_Struct(obj, gsl_sf_result, rslt);
  return rb_ary_new3(2, rb_float_new(rslt->val), rb_float_new(rslt->err));
}

static VALUE rb_gsl_sf_result_to_s(VALUE obj)
{
  gsl_sf_result *rslt = NULL;
  char str[32];
  Data_Get_Struct(obj, gsl_sf_result, rslt);
  sprintf(str, "%10.9e %10.9e", rslt->val, rslt->err);
  return rb_str_new2(str);
}

static VALUE rb_gsl_sf_result_e10_new(VALUE klass)
{
  gsl_sf_result_e10 *rslt = NULL;
  return Data_Make_Struct(cgsl_sf_result, gsl_sf_result_e10, 0, free, rslt);
}

static VALUE rb_gsl_sf_result_e10_val(VALUE obj)
{
  gsl_sf_result_e10 *rslt = NULL;
  Data_Get_Struct(obj, gsl_sf_result_e10, rslt);
  return rb_float_new(rslt->val);
}

static VALUE rb_gsl_sf_result_e10_err(VALUE obj)
{
  gsl_sf_result_e10 *rslt = NULL;
  Data_Get_Struct(obj, gsl_sf_result_e10, rslt);
  return rb_float_new(rslt->err);
}

static VALUE rb_gsl_sf_result_e10_e10(VALUE obj)
{
  gsl_sf_result_e10 *rslt = NULL;
  Data_Get_Struct(obj, gsl_sf_result_e10, rslt);
  return INT2FIX(rslt->e10);
}

static VALUE rb_gsl_sf_result_e10_to_a(VALUE obj)
{
  gsl_sf_result_e10 *rslt = NULL;
  Data_Get_Struct(obj, gsl_sf_result_e10, rslt);
  return rb_ary_new3(2, rb_float_new(rslt->val), rb_float_new(rslt->err));
}

static VALUE rb_gsl_sf_result_e10_to_s(VALUE obj)
{
  gsl_sf_result_e10 *rslt = NULL;
  char str[32];
  Data_Get_Struct(obj, gsl_sf_result_e10, rslt);
  sprintf(str, "%10.9e %10.9e\n", rslt->val, rslt->err);
  return rb_str_new2(str);
}

VALUE rb_gsl_sf_eval1(double (*func)(double), VALUE argv)
{
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(NUM2DBL(argv)));
    break;
  case T_ARRAY:
    return rb_gsl_ary_eval1(argv, func);
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      return rb_gsl_nary_eval1(argv, func);
    }
#endif
    if (MATRIX_P(argv)) {
      return matrix_eval_create(argv, func);
    } else if (VECTOR_P(argv)) {
      return vector_eval_create(argv, func);
    } else if (COMPLEX_P(argv) || VECTOR_COMPLEX_P(argv) || MATRIX_COMPLEX_P(argv)) {
      return rb_gsl_sf_eval_complex(func, argv);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s", rb_class2name(CLASS_OF(argv)));
    }
    break;
  }
}

VALUE rb_gsl_sf_eval_int_double(double (*func)(int, double), VALUE jj, VALUE argv)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  VALUE ary, xx;
  size_t i, j, k, n;
  double val;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  CHECK_FIXNUM(jj);
  j = FIX2INT(jj);
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(j, NUM2DBL(argv)));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      xx = rb_ary_entry(argv, i);
      Need_Float(xx);
      val = (*func)(j, NUM2DBL(xx));
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      argv = na_change_type(argv, NA_DFLOAT);
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)(j, ptr1[i]);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (k = 0; k < m->size2; k++) {
	  val = (*func)(j, gsl_matrix_get(m, i, k));
	  gsl_matrix_set(mnew, i, k, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
	val = (*func)(j, gsl_vector_get(v, i));
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval_double_int(double (*func)(double, int), VALUE argv, VALUE jj)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  VALUE ary, xx;
  size_t i, j, k, n;
  double val;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  CHECK_FIXNUM(jj);
  j = FIX2INT(jj);
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(NUM2DBL(argv), j));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      xx = rb_ary_entry(argv, i);
      Need_Float(xx);
      val = (*func)(NUM2DBL(xx), j);
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      argv = na_change_type(argv, NA_DFLOAT);
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)(ptr1[i], j);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (k = 0; k < m->size2; k++) {
	  val = (*func)(gsl_matrix_get(m, i, k), j);
	  gsl_matrix_set(mnew, i, k, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
	val = (*func)(gsl_vector_get(v, i), j);
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
  return Qnil;
}

VALUE rb_gsl_sf_eval_int_int_double(double (*func)(int, int, double), VALUE jj,
				    VALUE jj2, VALUE argv)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  VALUE ary, xx;
  size_t i, j, k, j2, n;
  double val;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  CHECK_FIXNUM(jj); CHECK_FIXNUM(jj2);
  j = FIX2INT(jj);
  j2 = FIX2INT(jj2);
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(j, j2, NUM2DBL(argv)));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      xx = rb_ary_entry(argv, i);
      Need_Float(xx);
      val = (*func)(j, j2, NUM2DBL(xx));
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      argv = na_change_type(argv, NA_DFLOAT);
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)(j, j2, ptr1[i]);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
        for (k = 0; k < m->size2; k++) {
          val = (*func)(j, j2, gsl_matrix_get(m, i, k));
          gsl_matrix_set(mnew, i, k, val);
        }
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
        val = (*func)(j, j2, gsl_vector_get(v, i));
        gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval_int_double_double(double (*func)(int, double, double), VALUE jj, 
				       VALUE ff, VALUE argv)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  VALUE ary, xx;
  size_t i, j, k, n;
  double f, val;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  CHECK_FIXNUM(jj);
  Need_Float(ff);
  j = FIX2INT(jj);
  f = NUM2DBL(ff);
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(j, f, NUM2DBL(argv)));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      xx = rb_ary_entry(argv, i);
      Need_Float(xx);
      val = (*func)(j, f, NUM2DBL(xx));
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      argv = na_change_type(argv, NA_DFLOAT);
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)(j, f, ptr1[i]);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (k = 0; k < m->size2; k++) {
	  val = (*func)(j, f, gsl_matrix_get(m, i, k));
	  gsl_matrix_set(mnew, i, k, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
	val = (*func)(j, f, gsl_vector_get(v, i));
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval_double_double(double (*func)(double, double), VALUE ff, VALUE argv)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  VALUE ary, xx;
  size_t i, k, n;
  double val, f;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  Need_Float(ff);
  f = NUM2DBL(ff);
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(f, NUM2DBL(argv)));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      xx = rb_ary_entry(argv, i);
      Need_Float(xx);
      val = (*func)(f, NUM2DBL(xx));
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      argv = na_change_type(argv, NA_DFLOAT);
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)(f, ptr1[i]);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (k = 0; k < m->size2; k++) {
	  gsl_matrix_set(mnew, i, k, (*func)(f, gsl_matrix_get(m, i, k)));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	gsl_vector_set(vnew, i, (*func)(f, gsl_vector_get(v, i)));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval_double3(double (*func)(double, double, double), 
			     VALUE ff, VALUE ff2, VALUE argv)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  VALUE ary, xx;
  size_t i, k, n;
  double val, f, f2;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  Need_Float(ff); Need_Float(ff2);
  f = NUM2DBL(ff);
  f2 = NUM2DBL(ff2);
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(f, f2, NUM2DBL(argv)));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      xx = rb_ary_entry(argv, i);
      Need_Float(xx);
      val = (*func)(f, f2, NUM2DBL(xx));
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      argv = na_change_type(argv, NA_DFLOAT);
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)(f, f2, ptr1[i]);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (k = 0; k < m->size2; k++) {
	  val = (*func)(f, f2, gsl_matrix_get(m, i, k));
	  gsl_matrix_set(mnew, i, k, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
	val = (*func)(f, f2, gsl_vector_get(v, i));
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval_double4(double (*func)(double, double, double, double), 
			     VALUE ff, VALUE ff2, VALUE ff3, VALUE argv)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  VALUE ary, xx;
  size_t i, k, n;
  double val, f, f2, f3;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  Need_Float(ff); Need_Float(ff2); Need_Float(ff3);
  f = NUM2DBL(ff);
  f2 = NUM2DBL(ff2);
  f3 = NUM2DBL(ff3);
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(f, f2, f3, NUM2DBL(argv)));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      xx = rb_ary_entry(argv, i);
      Need_Float(xx);
      val = (*func)(f, f2, f3, NUM2DBL(xx));
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      argv = na_change_type(argv, NA_DFLOAT);
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)(f, f2, f3, ptr1[i]);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (k = 0; k < m->size2; k++) {
	  val = (*func)(f, f2, f3, gsl_matrix_get(m, i, k));
	  gsl_matrix_set(mnew, i, k, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
	val = (*func)(f, f2, f3, gsl_vector_get(v, i));
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval1_int(double (*func)(int), VALUE argv)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  VALUE ary;
  size_t i, k, n;
  double val;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new((*func)(NUM2INT(argv)));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      val = (*func)(NUM2INT(rb_ary_entry(argv, i)));
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      argv = na_change_type(argv, NA_DFLOAT);
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)((int) ptr1[i]);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (k = 0; k < m->size2; k++) {
	  val = (*func)((int) gsl_matrix_get(m, i, k));
	  gsl_matrix_set(mnew, i, k, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
	val = (*func)((int) gsl_vector_get(v, i));
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval1_uint(double (*func)(unsigned int), VALUE argv)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
  VALUE ary;
  size_t i, k, n;
  double val;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new((*func)(NUM2UINT(argv)));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      val = (*func)(NUM2UINT(rb_ary_entry(argv, i)));
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      argv = na_change_type(argv, NA_DFLOAT);
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)((unsigned int) ptr1[i]);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (k = 0; k < m->size2; k++) {
	  val = (*func)((unsigned int) gsl_matrix_get(m, i, k));
	  gsl_matrix_set(mnew, i, k, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
	val = (*func)((unsigned int) gsl_vector_get(v, i));
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval_double_m(double (*func)(double, gsl_mode_t), VALUE argv, VALUE m)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *mm = NULL, *mnew = NULL;
  VALUE ary, xx;
  size_t i, k, n;
  double val;
  /*gsl_mode_t mode;
  char c;*/
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  switch (TYPE(m)) {
  case T_STRING:
    /*c = tolower(NUM2CHR(m));
    if (c == 'd') mode = GSL_PREC_DOUBLE;
    else if (c == 's') mode = GSL_PREC_SINGLE;
    else if (c == 'a') mode = GSL_PREC_APPROX;
    else mode = GSL_PREC_DOUBLE;*/
    break;
  case T_FIXNUM:
    /*mode = FIX2INT(m);*/
    break;
  default:
    rb_raise(rb_eArgError, "wrong type argument %s (String or Fixnum expected)",
	     rb_class2name(CLASS_OF(m)));
  }
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(NUM2DBL(argv), m));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      xx = rb_ary_entry(argv, i);
      Need_Float(xx);
      val = (*func)(NUM2DBL(xx), m);
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)(ptr1[i], m);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, mm);
      mnew = gsl_matrix_alloc(mm->size1, mm->size2);
      for (i = 0; i < mm->size1; i++) {
	for (k = 0; k < mm->size2; k++) {
	  val = (*func)(gsl_matrix_get(mm, i, k), m);
	  gsl_matrix_set(mnew, i, k, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
	val = (*func)(gsl_vector_get(v, i), m);
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval_double2_m(double (*func)(double, double, gsl_mode_t), 
			       VALUE argv, VALUE x2, VALUE m)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *mm = NULL, *mnew = NULL;
  VALUE ary, xx;
  size_t i, k, n;
  double val, xx2;
  /*gsl_mode_t mode;
  char c;*/
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  Need_Float(x2);
  xx2 = NUM2DBL(x2);
  /*c = tolower(NUM2CHR(m));
  if (c == 'd') mode = GSL_PREC_DOUBLE;
  else if (c == 's') mode = GSL_PREC_SINGLE;
  else if (c == 'a') mode = GSL_PREC_APPROX;
  else mode = GSL_PREC_DOUBLE;*/
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(NUM2DBL(argv), xx2, m));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      xx = rb_ary_entry(argv, i);
      Need_Float(xx);
      val = (*func)(NUM2DBL(xx), xx2, m);
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)(ptr1[i], xx2, m);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, mm);
      mnew = gsl_matrix_alloc(mm->size1, mm->size2);
      for (i = 0; i < mm->size1; i++) {
	for (k = 0; k < mm->size2; k++) {
	  val = (*func)(gsl_matrix_get(mm, i, k), xx2, m);
	  gsl_matrix_set(mnew, i, k, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
	val = (*func)(gsl_vector_get(v, i), xx2, m);
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval_double3_m(double (*func)(double, double, double, gsl_mode_t), 
			       VALUE argv, VALUE x2, VALUE x3, VALUE m)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *mm = NULL, *mnew = NULL;
  VALUE ary, xx;
  size_t i, k, n;
  double val, xx2, xx3;
  /*gsl_mode_t mode;
  char c;*/
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  Need_Float(x2); Need_Float(x3);
  xx2 = NUM2DBL(x2);
  xx3 = NUM2DBL(x3);
  /*c = tolower(NUM2CHR(m));
  if (c == 'd') mode = GSL_PREC_DOUBLE;
  else if (c == 's') mode = GSL_PREC_SINGLE;
  else if (c == 'a') mode = GSL_PREC_APPROX;
  else mode = GSL_PREC_DOUBLE;*/
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(NUM2DBL(argv), xx2, xx3, m));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      xx = rb_ary_entry(argv, i);
      Need_Float(xx);
      val = (*func)(NUM2DBL(xx), xx2, xx3, m);
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)(ptr1[i], xx2, xx3, m);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, mm);
      mnew = gsl_matrix_alloc(mm->size1, mm->size2);
      for (i = 0; i < mm->size1; i++) {
	for (k = 0; k < mm->size2; k++) {
	  val = (*func)(gsl_matrix_get(mm, i, k), xx2, xx3, m);
	  gsl_matrix_set(mnew, i, k, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
	val = (*func)(gsl_vector_get(v, i), xx2, xx3, m);
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval_double4_m(double (*func)(double, double, double, double,
					      gsl_mode_t), 
			       VALUE argv, VALUE x2, VALUE x3, VALUE x4, VALUE m)
{
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *mm = NULL, *mnew = NULL;
  VALUE ary, xx;
  size_t i, k, n;
  double val, xx2, xx3, xx4;
  /*gsl_mode_t mode;
  char c;*/
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  Need_Float(x2); Need_Float(x3); Need_Float(x4);
  xx2 = NUM2DBL(x2);  xx3 = NUM2DBL(x3); xx4 = NUM2DBL(x4);
  /*c = tolower(NUM2CHR(m));
  if (c == 'd') mode = GSL_PREC_DOUBLE;
  else if (c == 's') mode = GSL_PREC_SINGLE;
  else if (c == 'a') mode = GSL_PREC_APPROX;
  else mode = GSL_PREC_DOUBLE;*/
  if (CLASS_OF(argv) == rb_cRange) argv = rb_gsl_range2ary(argv);
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(NUM2DBL(argv), NUM2DBL(x2), NUM2DBL(x3),
				NUM2DBL(x4), m));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      xx = rb_ary_entry(argv, i);
      Need_Float(xx);
      val = (*func)(NUM2DBL(xx), xx2, xx3, xx4, m);
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)(ptr1[i], xx2, xx3, xx4, m);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, mm);
      mnew = gsl_matrix_alloc(mm->size1, mm->size2);
      for (i = 0; i < mm->size1; i++) {
	for (k = 0; k < mm->size2; k++) {
	  val = (*func)(gsl_matrix_get(mm, i, k), xx2, xx3, xx4, m);
	  gsl_matrix_set(mnew, i, k, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
	val = (*func)(gsl_vector_get(v, i),  xx2, xx3, xx4, m);
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval_e(int (*func)(double, gsl_sf_result*), VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(NUM2DBL(x), rslt);
  return v;
}

VALUE rb_gsl_sf_eval_e_int(int (*func)(int, gsl_sf_result*), VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  CHECK_FIXNUM(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(NUM2INT(x), rslt);
  return v;
}

VALUE rb_gsl_sf_eval_e_uint(int (*func)(unsigned int, gsl_sf_result*), VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(NUM2UINT(x), rslt);
  return v;
}

VALUE rb_gsl_sf_eval_e_int_uint(int (*func)(int, unsigned int, gsl_sf_result*), 
				VALUE n, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  CHECK_FIXNUM(n);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(FIX2INT(n), NUM2UINT(x), rslt);
  return v;
}

VALUE rb_gsl_sf_eval_e_double_uint(int (*func)(double, unsigned int, gsl_sf_result*), 
				VALUE y, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  Need_Float(y);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(NUM2DBL(y), NUM2UINT(x), rslt);
  return v;
}

VALUE rb_gsl_sf_eval_e_int_double(int (*func)(int, double, gsl_sf_result*), 
				  VALUE n, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  CHECK_FIXNUM(n);
  Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(FIX2INT(n), NUM2DBL(x), rslt);
  return v;
}

VALUE rb_gsl_sf_eval_e_int_double2(int (*func)(int, double, double, gsl_sf_result*), 
				  VALUE n, VALUE x1, VALUE x2)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  CHECK_FIXNUM(n);
  Need_Float(x1); Need_Float(x2);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(FIX2INT(n), NUM2DBL(x1), NUM2DBL(x2), rslt);
  return v;
}


VALUE rb_gsl_sf_eval_e_int_int_double(int (*func)(int, int, double, gsl_sf_result*), 
				      VALUE n1, VALUE n2, VALUE x)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  CHECK_FIXNUM(n1); CHECK_FIXNUM(n2);
  Need_Float(x);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(FIX2INT(n1), FIX2INT(n2), NUM2DBL(x), rslt);
  return v;
}

VALUE rb_gsl_sf_eval_e_double2(int (*func)(double, double, gsl_sf_result*), 
				  VALUE x1, VALUE x2)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  Need_Float(x1); Need_Float(x2); 
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(NUM2DBL(x1), NUM2DBL(x2), rslt);
  return v;
}


VALUE rb_gsl_sf_eval_e_double3(int (*func)(double, double, double, gsl_sf_result*), 
				  VALUE x1, VALUE x2, VALUE x3)
{
  gsl_sf_result *rslt = NULL;
  VALUE v;
  Need_Float(x1); Need_Float(x2); Need_Float(x3);
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(NUM2DBL(x1), NUM2DBL(x2),NUM2DBL(x3), rslt);
  return v;
}

VALUE rb_gsl_sf_eval_e_m(int (*func)(double, gsl_mode_t, gsl_sf_result*), 
			 VALUE x, VALUE m)
{
  gsl_mode_t mode;
  char c;
  gsl_sf_result *rslt = NULL;
  VALUE v;
  Need_Float(x);
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
    rb_raise(rb_eArgError, "wrong type argument %s (String or Fixnum expected)",
	     rb_class2name(CLASS_OF(m)));
    break;
  }
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(NUM2DBL(x), mode, rslt);
  return v;
}


VALUE rb_gsl_sf_eval_e_double2_m(int (*func)(double, double, gsl_mode_t, gsl_sf_result*), 
			 VALUE x1, VALUE x2, VALUE m)
{
  gsl_mode_t mode;
  char c;
  gsl_sf_result *rslt = NULL;
  VALUE v;
  Need_Float(x1);  Need_Float(x2);
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
    rb_raise(rb_eArgError, "wrong type argument %s (String or Fixnum expected)",
	     rb_class2name(CLASS_OF(m)));
    break;
  }
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(NUM2DBL(x1), NUM2DBL(x2), mode, rslt);
  return v;
}

VALUE rb_gsl_sf_eval_e_double3_m(int (*func)(double, double, double, gsl_mode_t, gsl_sf_result*), 
			 VALUE x1, VALUE x2, VALUE x3, VALUE m)
{
  gsl_mode_t mode;
  char c;
  gsl_sf_result *rslt = NULL;
  VALUE v;
  Need_Float(x1); Need_Float(x2); Need_Float(x3);
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
    rb_raise(rb_eArgError, "wrong type argument %s (String or Fixnum expected)",
	     rb_class2name(CLASS_OF(m)));
    break;
  }
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(NUM2DBL(x1), NUM2DBL(x2),NUM2DBL(x3), mode, rslt);
  return v;
}


VALUE rb_gsl_sf_eval_e_double4_m(int (*func)(double, double, double, double, gsl_mode_t, gsl_sf_result*), 
			 VALUE x1, VALUE x2, VALUE x3, VALUE x4, VALUE m)
{
  gsl_mode_t mode;
  char c;
  gsl_sf_result *rslt = NULL;
  VALUE v;
  Need_Float(x1); Need_Float(x2); Need_Float(x3); Need_Float(x4);
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
    rb_raise(rb_eArgError, "wrong type argument %s (String or Fixnum expected)",
	     rb_class2name(CLASS_OF(m)));
    break;
  }    
  v = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, rslt);
  (*func)(NUM2DBL(x1), NUM2DBL(x2),NUM2DBL(x3), NUM2DBL(x4), mode, rslt);
  return v;
}


VALUE eval_sf(double (*func)(double, gsl_mode_t), VALUE argv)
{
  VALUE ary, xx;
  size_t i, k, n;
  double val;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *mm = NULL, *mnew = NULL;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  switch (TYPE(argv)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    return rb_float_new((*func)(NUM2DBL(argv), GSL_PREC_DOUBLE));
    break;
  case T_ARRAY:
    n = RARRAY_LEN(argv);
    ary = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      xx = rb_ary_entry(argv, i);
      Need_Float(xx);
      val = (*func)(NUM2DBL(xx), GSL_PREC_DOUBLE);
      rb_ary_store(ary, i, rb_float_new(val));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(argv)) {
      ptr1 = NA_PTR_TYPE(argv, double*);
      GetNArray(argv, na);
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < n; i++) ptr2[i] = (*func)(ptr1[i], GSL_PREC_DOUBLE);
      return ary;
    }
#endif
    if (rb_obj_is_kind_of(argv, cgsl_matrix)) {
      Data_Get_Struct(argv, gsl_matrix, mm);
      mnew = gsl_matrix_alloc(mm->size1, mm->size2);
      for (i = 0; i < mm->size1; i++) {
	for (k = 0; k < mm->size2; k++) {
	  val = (*func)(gsl_matrix_get(mm, i, k), GSL_PREC_DOUBLE);
	  gsl_matrix_set(mnew, i, k, val);
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      CHECK_VECTOR(argv);
      Data_Get_Struct(argv, gsl_vector, v);
      n = v->size;
      vnew = gsl_vector_alloc(n);
      for (i = 0; i < n; i++) {
	val = (*func)(gsl_vector_get(v, i), GSL_PREC_DOUBLE);
	gsl_vector_set(vnew, i, val);
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    }
    break;
  }
}

VALUE rb_gsl_sf_eval_complex(double (*f)(double), VALUE obj)
{
  gsl_complex *z, *znew, c;
  gsl_vector_complex *v, *vnew;
  gsl_matrix_complex *m, *mnew;
  size_t i, j;
  if (COMPLEX_P(obj)) {
    Data_Get_Struct(obj, gsl_complex, z);
    znew = xmalloc(sizeof(gsl_complex));
    GSL_SET_REAL(znew, (*f)(GSL_REAL(*z)));
    GSL_SET_IMAG(znew, (*f)(GSL_IMAG(*z)));
    return Data_Wrap_Struct(cgsl_complex, 0, free, znew);
  } else if (VECTOR_COMPLEX_P(obj)) {
    Data_Get_Struct(obj, gsl_vector_complex, v);
    vnew = gsl_vector_complex_alloc(v->size);
    for (i = 0; i < v->size; i++) {
      z = GSL_COMPLEX_AT(v, i);
      GSL_SET_REAL(&c, (*f)(GSL_REAL(*z)));
      GSL_SET_IMAG(&c, (*f)(GSL_IMAG(*z)));
      gsl_vector_complex_set(vnew, i, c);
    }
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
  } else if (MATRIX_COMPLEX_P(obj)) {
    Data_Get_Struct(obj, gsl_matrix_complex, m);
    mnew = gsl_matrix_complex_alloc(m->size1, m->size2);
    for (i = 0; i < m->size1; i++) {
      for (j = 0; j < m->size2; j++) {
	c = gsl_matrix_complex_get(m, i, j);
	GSL_SET_REAL(&c, (*f)(GSL_REAL(c)));
	GSL_SET_IMAG(&c, (*f)(GSL_IMAG(c)));
	gsl_matrix_complex_set(mnew, i, j, c);
      }
    }
    return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, mnew);
  } else {
    rb_raise(rb_eTypeError, 
	     "wrong argument type %s "
	     " (GSL::Complex or GSL::Vector::Complex expected)", 
	     rb_class2name(CLASS_OF(obj)));
  }
}

void Init_gsl_sf(VALUE module)
{
  VALUE mgsl_sf;
  mgsl_sf = rb_define_module_under(module, "Sf");

  cgsl_sf_result = rb_define_class_under(mgsl_sf, "Result", 
					 cGSL_Object);
  rb_define_singleton_method(cgsl_sf_result, "new", rb_gsl_sf_result_new,
			     0);
  rb_define_method(cgsl_sf_result, "print", rb_gsl_sf_result_print, 0);
  rb_define_method(cgsl_sf_result, "inspect", rb_gsl_sf_result_inspect, 0);
  rb_define_method(cgsl_sf_result, "val", rb_gsl_sf_result_val, 0);
  rb_define_method(cgsl_sf_result, "err", rb_gsl_sf_result_err, 0);
  rb_define_method(cgsl_sf_result, "to_a", rb_gsl_sf_result_to_a, 0);
  rb_define_method(cgsl_sf_result, "to_s", rb_gsl_sf_result_to_s, 0);

  cgsl_sf_result_e10 = rb_define_class_under(mgsl_sf, "Result_e10", 
					     cGSL_Object);
  rb_define_singleton_method(cgsl_sf_result_e10, "new", 
			     rb_gsl_sf_result_e10_new, 0);
  rb_define_method(cgsl_sf_result_e10, "val", rb_gsl_sf_result_e10_val, 0);
  rb_define_method(cgsl_sf_result_e10, "err", rb_gsl_sf_result_e10_err, 0);
  rb_define_method(cgsl_sf_result_e10, "e10", rb_gsl_sf_result_e10_e10, 0);
  rb_define_method(cgsl_sf_result_e10, "to_a", rb_gsl_sf_result_e10_to_a, 0);
  rb_define_method(cgsl_sf_result_e10, "to_s", rb_gsl_sf_result_e10_to_s, 0);

  Init_gsl_sf_airy(mgsl_sf);
  Init_gsl_sf_bessel(mgsl_sf);
  Init_gsl_sf_clausen(mgsl_sf);
  Init_gsl_sf_coulomb(mgsl_sf);
  Init_gsl_sf_coupling(mgsl_sf);
  Init_gsl_sf_dawson(mgsl_sf);
  Init_gsl_sf_debye(mgsl_sf);
  Init_gsl_sf_dilog(mgsl_sf);
  Init_gsl_sf_elementary(mgsl_sf);
  Init_gsl_sf_ellint(mgsl_sf);
  Init_gsl_sf_elljac(mgsl_sf);
  Init_gsl_sf_erfc(mgsl_sf);
  Init_gsl_sf_exp(mgsl_sf);
  Init_gsl_sf_expint(mgsl_sf);
  Init_gsl_sf_fermi_dirac(mgsl_sf);
  Init_gsl_sf_gamma(mgsl_sf);
  Init_gsl_sf_gegenbauer(mgsl_sf);
  Init_gsl_sf_hyperg(mgsl_sf);
  Init_gsl_sf_laguerre(mgsl_sf);
  Init_gsl_sf_lambert(mgsl_sf);
  Init_gsl_sf_legendre(mgsl_sf);
  Init_gsl_sf_log(mgsl_sf);
  Init_gsl_sf_power(mgsl_sf);
  Init_gsl_sf_psi(mgsl_sf);
  Init_gsl_sf_synchrotron(mgsl_sf);
  Init_gsl_sf_transport(mgsl_sf);
  Init_gsl_sf_trigonometric(mgsl_sf);
  Init_gsl_sf_zeta(mgsl_sf);
  
#ifdef GSL_1_9_LATER  
	Init_sf_mathieu(mgsl_sf);
#endif
}
