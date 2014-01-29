/*
  poly_source.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#ifdef BASE_DOUBLE
#define NUMCONV(x) NUM2DBL(x)
#define NUMCONV2(x) NUM2DBL(x)
#define MAT_ROW_COL MATRIX_ROW_COL
#define MAT_P MATRIX_P
#define MAT_ROW_P MATRIX_ROW_P
#define MAT_COL_P MATRIX_COL_P
#define C_TO_VALUE rb_float_new
#define C_TO_VALUE2 rb_float_new
#define CHECK_MAT CHECK_MATRIX
#define MAT_VIEW_P MATRIX_VIEW_P
#define PRINTF_FORMAT "%4.3e "
#define VEC_ROW_COL VECTOR_ROW_COL
#define VEC_P VECTOR_P
#define VEC_ROW_P VECTOR_ROW_P
#define VEC_COL_P VECTOR_COL_P
#define CHECK_VEC CHECK_VECTOR
#define VEC_VIEW_P VECTOR_VIEW_P
#elif defined(BASE_INT)
#define NUMCONV(x) NUM2DBL(x)
#define NUMCONV2(x) NUM2INT(x)
#define PRINTF_FORMAT "%d "
#define MAT_ROW_COL MATRIX_INT_ROW_COL
#define MAT_P MATRIX_INT_P
#define C_TO_VALUE INT2FIX
#define C_TO_VALUE2 INT2NUM
#define MAT_ROW_P MATRIX_INT_ROW_P
#define MAT_COL_P MATRIX_INT_COL_P
#define CHECK_MAT CHECK_MATRIX_INT
#define MAT_VIEW_P MATRIX_INT_VIEW_P
#define VEC_ROW_COL VECTOR_INT_ROW_COL
#define VEC_P VECTOR_INT_P
#define VEC_ROW_P VECTOR_INT_ROW_P
#define VEC_COL_P VECTOR_INT_COL_P
#define CHECK_VEC CHECK_VECTOR_INT
#define VEC_VIEW_P VECTOR_INT_VIEW_P
#endif

#ifdef BASE_DOUBLE
VALUE cgsl_poly, cgsl_poly_int;
VALUE cgsl_poly_workspace;
VALUE cgsl_poly_complex_workspace;
#ifdef GSL_1_1_LATER
VALUE cgsl_poly_dd;
VALUE cgsl_poly_taylor;
#endif
#endif

#ifdef BASE_INT
double gsl_poly_int_eval(const BASE c[], const int len, const double x)
{
  int i;
  double ans = (double) c[len-1];
  for(i=len-1; i>0; i--) ans = (double) c[i-1] + x * ans;
  return ans;
}
#endif
#ifdef BASE_DOUBLE
#ifdef HAVE_NARRAY_H
#include "include/rb_gsl_with_narray.h"
#endif
#ifdef GSL_1_11_LATER
static VALUE rb_gsl_complex_poly_complex_eval(VALUE a, VALUE b);
#endif
static VALUE rb_gsl_poly_eval_singleton(VALUE klass, VALUE a, VALUE x)
{
  gsl_vector *v = NULL, *vx, *vnew;
  gsl_matrix *mx, *mnew;
  double rslt;
  int flag = 0;
  size_t i, N, n;
  VALUE val;
  double *ptr0;
  double *ptr1, *ptr2;  
#ifdef HAVE_NARRAY_H
  int shape[1];
#endif  
#ifdef GSL_1_11_LATER
  gsl_complex *z, zz;
  gsl_vector_complex *vz, *vznew;
  if (rb_obj_is_kind_of(a, cgsl_vector_complex)) 
    return rb_gsl_complex_poly_complex_eval(a, x);
#endif
  switch (TYPE(a)) {
  case T_ARRAY:
#ifdef GSL_1_11_LATER
    if (rb_obj_is_kind_of(rb_ary_entry(a, 0), cgsl_complex)) 
      return rb_gsl_complex_poly_complex_eval(a, x);
#endif
    v = make_cvector_from_rarray(a);
    N = v->size;
    ptr0 = v->data;
    flag = 1;
    break;
  default:
    if (rb_obj_is_kind_of(a, cgsl_vector)) {
      Data_Get_Struct(a, gsl_vector, v);
      N = v->size;
      ptr0 = v->data;
#ifdef HAVE_NARRAY_H
    } else if (NA_IsNArray(a)) {
      N = NA_TOTAL(a);
      ptr0 = NA_PTR_TYPE(a, double*);
#endif
    } else {
      rb_raise(rb_eTypeError, "Wrong argument type %s (Array, GSL::Vector or NArray expected)",
        rb_class2name(CLASS_OF(a)));
    }
  }
  switch (TYPE(x)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    val = rb_float_new(gsl_poly_eval(ptr0, N, NUM2DBL(x)));
    break;    
  case T_ARRAY:
    n = RARRAY_LEN(x);
    val = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
      rslt = gsl_poly_eval(ptr0, N, NUM2DBL(rb_ary_entry(x, i)));
      rb_ary_store(val, i, rb_float_new(rslt));
    }
    break;  
  default:
    if (rb_obj_is_kind_of(x, cgsl_vector)) {
      Data_Get_Struct(x, gsl_vector, vx);
      vnew = gsl_vector_alloc(vx->size);
      val = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
      n = vx->size;
      ptr1 = vx->data;
      ptr2 = vnew->data;
    } else if (rb_obj_is_kind_of(x, cgsl_matrix)) {
      Data_Get_Struct(x, gsl_matrix, mx);
      mnew = gsl_matrix_alloc(mx->size1, mx->size2);
      val = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);      
      n = mx->size1*mx->size2;
      ptr1 = mx->data;
      ptr2 = mnew->data;
#ifdef HAVE_NARRAY_H
    } else if (NA_IsNArray(x)) {
      shape[0] = NA_TOTAL(x);
      n = shape[0];
      val = na_make_object(NA_DFLOAT, 1, shape, CLASS_OF(x));
      ptr1 = NA_PTR_TYPE(x, double*);
      ptr2 = NA_PTR_TYPE(val, double*);      
#endif
#ifdef GSL_1_11_LATER
    } else if (rb_obj_is_kind_of(x, cgsl_complex)) {
      Data_Get_Struct(x, gsl_complex, z);
      zz = gsl_poly_complex_eval(ptr0, N, *z);
      z = make_complex(GSL_REAL(zz), GSL_IMAG(zz));
      if (flag == 1) gsl_vector_free(v);
      return Data_Wrap_Struct(cgsl_complex, 0, free, z);
    } else if (rb_obj_is_kind_of(x, cgsl_vector_complex)) {
      Data_Get_Struct(x, gsl_vector_complex, vz);
      vznew = gsl_vector_complex_alloc(vz->size);
      for (i = 0; i < vz->size; i++) {
	zz = gsl_poly_complex_eval(ptr0, N, gsl_vector_complex_get(vz, i));
	gsl_vector_complex_set(vznew, i, zz);
      }
      if (flag == 1) gsl_vector_free(v);
      return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vznew);
#endif
    } else {
      rb_raise(rb_eTypeError, "Wrong argument type %s (A number, Array, GSL::Vector or NArray expected)",
        rb_class2name(CLASS_OF(a)));      
    }
    for (i = 0; i < n; i++) {
      ptr2[i] = gsl_poly_eval(ptr0, N, ptr1[i]);
    }
  }
  if (flag == 1) gsl_vector_free(v);
  return val;
}
#ifdef GSL_1_11_LATER
static VALUE rb_gsl_complex_poly_complex_eval(VALUE a, VALUE b)
{
  gsl_vector_complex *coef = NULL, *zb, *vnew;
  gsl_complex *zc;
  gsl_complex z, *zx, *res;
  VALUE ret;
  size_t i, N;
  int flag = 0;
  if (rb_obj_is_kind_of(a, cgsl_vector_complex)) {
    Data_Get_Struct(a, gsl_vector_complex, coef);
    N = coef->size;
    zc = (gsl_complex*) coef->data;
  } else if (TYPE(a) == T_ARRAY) {
    N = RARRAY_LEN(a);
    zc = (gsl_complex*) malloc(sizeof(gsl_complex));
    flag = 1;
    for (i = 0; i < N; i++) {
      Data_Get_Struct(rb_ary_entry(a, i), gsl_complex, zx);
      zc[i] = *zx;
    }
  } else {
    rb_raise(rb_eTypeError, "rb_gsl_complex_poly_complex_solve: wrong argument type %s (GSL::Vector::Complex or Array expected)\n", rb_class2name(CLASS_OF(a)));
  }

  switch (TYPE(b)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    res = (gsl_complex*) malloc(sizeof(gsl_complex));
    ret = Data_Wrap_Struct(cgsl_complex, 0, free, res);
    GSL_SET_REAL(&z, NUM2DBL(b));
    GSL_SET_IMAG(&z, 0.0);
    *res = gsl_complex_poly_complex_eval(zc, coef->size, z);
    break;
  case T_ARRAY:
    ret = rb_ary_new2(RARRAY_LEN(b));
    for (i = 0; (int) i < RARRAY_LEN(b); i++) {
      Data_Get_Struct(rb_ary_entry(b, i), gsl_complex, zx);
      res = (gsl_complex*) malloc(sizeof(gsl_complex));
      *res = gsl_complex_poly_complex_eval(zc, N, *zx);
      rb_ary_store(ret, i, Data_Wrap_Struct(cgsl_complex, 0, free, res));
    }
    break;
  default:
    if (rb_obj_is_kind_of(b, cgsl_complex)) {
      res = (gsl_complex*) malloc(sizeof(gsl_complex));
      ret = Data_Wrap_Struct(cgsl_complex, 0, free, res);
      Data_Get_Struct(b, gsl_complex, zx);
      *res = gsl_complex_poly_complex_eval(zc, N, *zx);
    } else if (rb_obj_is_kind_of(b, cgsl_vector_complex)) {
      Data_Get_Struct(b, gsl_vector_complex, zb);
      vnew = gsl_vector_complex_alloc(zb->size);
      ret = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
      for (i = 0; i < zb->size; i++) {
	z = gsl_vector_complex_get(zb, i);
	gsl_vector_complex_set(vnew, i, gsl_complex_poly_complex_eval(zc, N, z));
      }
    } else {
      rb_raise(rb_eTypeError, "Wrong argument type %s.\n", rb_class2name(CLASS_OF(b)));
    }
  }
  if (flag == 1) free(zc);
  return ret;
}
#endif
#endif

static VALUE FUNCTION(rb_gsl_poly,eval)(VALUE obj, VALUE xx)
{
  GSL_TYPE(gsl_poly) *p = NULL;
  GSL_TYPE(gsl_vector) *v = NULL;
  GSL_TYPE(gsl_matrix) *m = NULL;
  gsl_vector *vnew = NULL;
  gsl_matrix *mnew = NULL;
  VALUE x, ary;
  size_t i, j;
#ifdef BASE_DOUBLE
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2;
  size_t n;
#endif
#ifdef GSL_1_11_LATER
  gsl_complex *z, zz;
  gsl_vector_complex *vz, *vznew;
#endif
#endif

  Data_Get_Struct(obj, GSL_TYPE(gsl_poly), p);
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new(FUNCTION(gsl_poly,eval)(p->data, p->size, NUM2DBL(xx)));
    break;
  case T_ARRAY:
    ary = rb_ary_new2(RARRAY_LEN(xx));
    for (i = 0; (int) i < RARRAY_LEN(xx); i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(x);
      rb_ary_store(ary, i, rb_float_new(FUNCTION(gsl_poly,eval)(p->data, p->size, NUM2DBL(x))));
    }
    return ary;
    break;
  default:
#ifdef BASE_DOUBLE
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      GetNArray(xx, na);
      ptr1 = (double*) na->ptr;
      n = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary,double*);
      for (i = 0; i < n; i++) ptr2[i] = FUNCTION(gsl_poly,eval)(p->data,p->size,ptr1[i]);
      return ary;
    }
#endif
#endif
    if (VEC_P(xx)) {
      Data_Get_Struct(xx, GSL_TYPE(gsl_vector), v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	gsl_vector_set(vnew, i, FUNCTION(gsl_poly,eval)(p->data, p->size, FUNCTION(gsl_vector,get)(v, i)));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MAT_P(xx)) {
      Data_Get_Struct(xx, GSL_TYPE(gsl_matrix), m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  gsl_matrix_set(mnew, i, j, 
			 FUNCTION(gsl_poly,eval)(p->data, p->size, FUNCTION(gsl_matrix,get)(m, i, j)));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
#ifdef BASE_DOUBLE
#ifdef GSL_1_11_LATER
    } else if (rb_obj_is_kind_of(xx, cgsl_complex)) {
      Data_Get_Struct(xx, gsl_complex, z);
      zz = gsl_poly_complex_eval(p->data, p->size, *z);
      z = make_complex(GSL_REAL(zz), GSL_IMAG(zz));
      return Data_Wrap_Struct(cgsl_complex, 0, free, z);
    } else if (rb_obj_is_kind_of(xx, cgsl_vector_complex)) {
      Data_Get_Struct(xx, gsl_vector_complex, vz);
      vznew = gsl_vector_complex_alloc(vz->size);
      for (i = 0; i < vz->size; i++) {
	zz = gsl_poly_complex_eval(p->data, p->size, gsl_vector_complex_get(vz, i));
	gsl_vector_complex_set(vznew, i, zz);
      }
      return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vznew);
#endif
#endif
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  }
  return Qnil; /* never reach here */
}

/* singleton method */
static VALUE FUNCTION(rb_gsl_poly,eval2)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_poly) *p = NULL;
  GSL_TYPE(gsl_vector) *v = NULL;
  GSL_TYPE(gsl_matrix) *m = NULL;
  gsl_vector  *vnew = NULL;
  gsl_matrix  *mnew = NULL;
  VALUE xx, x, ary;
  size_t i, j, size;
#ifdef BASE_DOUBLE
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2;
#endif
#endif
  switch (argc) {
  case 2:
    Data_Get_Struct(argv[0], GSL_TYPE(gsl_poly), p);
    size = p->size;
    xx = argv[1];
    break;
  case 3:
    Data_Get_Struct(argv[0], GSL_TYPE(gsl_poly), p);
    size = FIX2INT(argv[1]);
    xx = argv[2];
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new(FUNCTION(gsl_poly,eval)(p->data, size, NUM2DBL(xx)));
    break;
  case T_ARRAY:
    ary = rb_ary_new2(RARRAY_LEN(xx));
    for (i = 0; (int) i < RARRAY_LEN(xx); i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(x);
      rb_ary_store(ary, i, rb_float_new(FUNCTION(gsl_poly,eval)(p->data, size, NUM2DBL(x))));
    }
    return ary;
    break;
  default:
#ifdef BASE_DOUBLE
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      GetNArray(xx, na);
      ptr1 = (double*) na->ptr;
      size = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary,double*);
      for (i = 0; i < size; i++) ptr2[i] = FUNCTION(gsl_poly,eval)(p->data,p->size,ptr1[i]);
      return ary;
    }
#endif
#endif
    if (VEC_P(xx)) {
      Data_Get_Struct(xx, GSL_TYPE(gsl_vector), v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	gsl_vector_set(vnew, i, FUNCTION(gsl_poly,eval)(p->data, size, FUNCTION(gsl_vector,get)(v, i)));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MAT_P(xx)) {
      Data_Get_Struct(xx, GSL_TYPE(gsl_matrix), m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  gsl_matrix_set(mnew, i, j, 
			 FUNCTION(gsl_poly,eval)(p->data, size, FUNCTION(gsl_matrix,get)(m, i, j)));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  }
  return Qnil; /* never reach here */
}

/* ax*x + bx + c = 0 */
static VALUE FUNCTION(rb_gsl_poly,solve_quadratic)(int argc, VALUE *argv, VALUE obj)
{
  double x0, x1;
  GSL_TYPE(gsl_poly) *v = NULL;
  gsl_vector *r;
  int n;
  switch (argc) {
  case 3:
    n = gsl_poly_solve_quadratic(NUMCONV2(argv[0]), NUMCONV2(argv[1]), NUMCONV2(argv[2]),
				 &x0, &x1);
    break;
  case 1:
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      n = gsl_poly_solve_quadratic(NUMCONV2(rb_ary_entry(argv[0], 0)), 
				   NUMCONV2(rb_ary_entry(argv[0], 1)),
				   NUMCONV2(rb_ary_entry(argv[0], 2)),
				   &x0, &x1);
      break;
    default:
      CHECK_VEC(argv[0]);
      Data_Get_Struct(argv[0], GSL_TYPE(gsl_poly), v);
      n = gsl_poly_solve_quadratic(FUNCTION(gsl_vector,get)(v, 0),
				   FUNCTION(gsl_vector,get)(v, 1),
				   FUNCTION(gsl_vector,get)(v, 2),
				   &x0, &x1);
      break;
    }
    break;
  default:
    rb_raise(rb_eArgError, 
	     "wrong number of arguments (3 numbers or 1 array or 1 vector)");
    break;
  }
  // If n == 0, we want to return an empty gsl_vector, but gsl_vectors can'y be
  // empty, so just return an empty Array.
  if(n == 0) {
    return rb_ary_new();
  }
  r = gsl_vector_alloc(n);
  switch(n) {
    case 2: gsl_vector_set(r, 1, x1); /* fall through */
    case 1: gsl_vector_set(r, 0, x0);
  }
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, r);
}

static VALUE FUNCTION(rb_gsl_poly,complex_solve_quadratic)(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex z0, z1;
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector_complex *r = NULL;
  int n;
  switch (argc) {
  case 3:
    n = gsl_poly_complex_solve_quadratic(NUMCONV2(argv[0]), 
					 NUMCONV2(argv[1]), NUMCONV2(argv[2]),
					 &z0, &z1);
    break;
  case 1:
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      n = gsl_poly_complex_solve_quadratic(NUMCONV2(rb_ary_entry(argv[0], 0)), 
					   NUMCONV2(rb_ary_entry(argv[0], 1)),
					   NUMCONV2(rb_ary_entry(argv[0], 2)),
					   &z0, &z1);
      break;
    default:
      CHECK_VEC(argv[0]);
      Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), v);
      n = gsl_poly_complex_solve_quadratic(FUNCTION(gsl_vector,get)(v, 0),
					   FUNCTION(gsl_vector,get)(v, 1),
					   FUNCTION(gsl_vector,get)(v, 2),
					   &z0, &z1);
      break;
    }
    break;
  default:
    rb_raise(rb_eArgError, 
	     "wrong number of arguments (3 numbers or 1 array or 1 vector)");
    break;
  }
  // If n == 0, we want to return an empty gsl_vector, but gsl_vectors can'y be
  // empty, so just return an empty Array.
  if(n == 0) {
    return rb_ary_new();
  }
  r = gsl_vector_complex_alloc(n);
  switch(n) {
    case 2: gsl_vector_complex_set(r, 1, z1); /* fall through */
    case 1: gsl_vector_complex_set(r, 0, z0);
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, r);
}

static VALUE FUNCTION(rb_gsl_poly,solve_cubic)(int argc, VALUE *argv, VALUE obj)
{
  double x0, x1, x2;
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector *r = NULL;
  int n;
  switch (argc) {
  case 3:
    n = gsl_poly_solve_cubic(NUMCONV2(argv[0]), NUMCONV2(argv[1]), NUMCONV2(argv[2]),
			     &x0, &x1, &x2);
    break;
  case 1:
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      n = gsl_poly_solve_cubic(NUMCONV2(rb_ary_entry(argv[0], 0)), 
				   NUMCONV2(rb_ary_entry(argv[0], 1)),
				   NUMCONV2(rb_ary_entry(argv[0], 2)),
				   &x0, &x1, &x2);
      break;
    default:
      CHECK_VEC(argv[0]);
      Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), v);
      n = gsl_poly_solve_cubic(FUNCTION(gsl_vector,get)(v, 0),
				   FUNCTION(gsl_vector,get)(v, 1),
				   FUNCTION(gsl_vector,get)(v, 2),
				   &x0, &x1, &x2);
      break;
    }
    break;
  default:
    rb_raise(rb_eArgError, 
	     "wrong number of arguments (3 numbers or 1 array or 1 vector)");
    break;
  }
  r = gsl_vector_alloc(n);
  switch(n) {
    case 3: gsl_vector_set(r, 2, x2); /* fall through */
    case 2: gsl_vector_set(r, 1, x1); /* fall through */
    case 1: gsl_vector_set(r, 0, x0);
  }
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, r);
}

static VALUE FUNCTION(rb_gsl_poly,complex_solve_cubic)(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex z0, z1, z2;
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector_complex *r = NULL;
  int n;
  switch (argc) {
  case 3:
    n = gsl_poly_complex_solve_cubic(NUMCONV2(argv[0]), NUMCONV2(argv[1]), NUMCONV2(argv[2]),
				     &z0, &z1, &z2);
    break;
  case 1:
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      n = gsl_poly_complex_solve_cubic(NUMCONV2(rb_ary_entry(argv[0], 0)), 
				   NUMCONV2(rb_ary_entry(argv[0], 1)),
				   NUMCONV2(rb_ary_entry(argv[0], 2)),
				   &z0, &z1, &z2);
      break;
    default:
      CHECK_VEC(argv[0]);
      Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), v);
      n = gsl_poly_complex_solve_cubic(FUNCTION(gsl_vector,get)(v, 0),
				   FUNCTION(gsl_vector,get)(v, 1),
				   FUNCTION(gsl_vector,get)(v, 2),
				   &z0, &z1, &z2);
      break;
    }
    break;
  default:
    rb_raise(rb_eArgError, 
	     "wrong number of arguments (3 numbers or 1 array or 1 vector)");
    break;
  }
  r = gsl_vector_complex_alloc(n);
  switch(n) {
    case 3: gsl_vector_complex_set(r, 2, z2); /* fall through */
    case 2: gsl_vector_complex_set(r, 1, z1); /* fall through */
    case 1: gsl_vector_complex_set(r, 0, z0);
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, r);
}

#ifdef HAVE_GSL_POLY_SOLVE_QUARTIC
static VALUE FUNCTION(rb_gsl_poly,solve_quartic)(int argc, VALUE *argv, VALUE obj)
{
  double x0, x1, x2, x3;
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector *r = NULL;
  int n;
  switch (argc) {
  case 4:
    n = gsl_poly_solve_quartic(NUMCONV2(argv[0]), NUMCONV2(argv[1]), NUMCONV2(argv[2]),
			       NUMCONV2(argv[3]),
			       &x0, &x1, &x2, &x3);
    break;
  case 1:
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      n = gsl_poly_solve_quartic(NUMCONV2(rb_ary_entry(argv[0], 0)), 
			       NUMCONV2(rb_ary_entry(argv[0], 1)),
			       NUMCONV2(rb_ary_entry(argv[0], 2)),
			       NUMCONV2(rb_ary_entry(argv[0], 3)),
			       &x0, &x1, &x2, &x3);
      break;
    default:
      CHECK_VEC(argv[0]);
      Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), v);
      n = gsl_poly_solve_quartic(FUNCTION(gsl_vector,get)(v, 0),
			       FUNCTION(gsl_vector,get)(v, 1),
			       FUNCTION(gsl_vector,get)(v, 2),
			       FUNCTION(gsl_vector,get)(v, 3),
			       &x0, &x1, &x2, &x3);
      break;
    }
    break;
  default:
    rb_raise(rb_eArgError, 
	     "wrong number of arguments (3 numbers or 1 array or 1 vector)");
    break;
  }
  r = gsl_vector_alloc(n);
  switch(n) {
    case 4: gsl_vector_set(r, 3, x3); /* fall through */
    case 3: gsl_vector_set(r, 2, x2); /* fall through */
    case 2: gsl_vector_set(r, 1, x1); /* fall through */
    case 1: gsl_vector_set(r, 0, x0);
  }
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, r);
}

static VALUE FUNCTION(rb_gsl_poly,complex_solve_quartic)(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex z0, z1, z2, z3;
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector_complex *r = NULL;
  int n;
  switch (argc) {
  case 4:
    n = gsl_poly_complex_solve_quartic(NUMCONV2(argv[0]), NUMCONV2(argv[1]), 
				       NUMCONV2(argv[2]), NUMCONV2(argv[3]),
				       &z0, &z1, &z2, &z3);
    break;
  case 1:
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      n = gsl_poly_complex_solve_quartic(NUMCONV2(rb_ary_entry(argv[0], 0)), 
					 NUMCONV2(rb_ary_entry(argv[0], 1)),
					 NUMCONV2(rb_ary_entry(argv[0], 2)),
					 NUMCONV2(rb_ary_entry(argv[0], 3)),
					 &z0, &z1, &z2, &z3);
      break;
    default:
      CHECK_VEC(argv[0]);
      Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), v);
      n = gsl_poly_complex_solve_quartic(FUNCTION(gsl_vector,get)(v, 0),
					 FUNCTION(gsl_vector,get)(v, 1),
					 FUNCTION(gsl_vector,get)(v, 2),
					 FUNCTION(gsl_vector,get)(v, 3),
					 &z0, &z1, &z2, &z3);
      break;
    }
    break;
  default:
    rb_raise(rb_eArgError, 
	     "wrong number of arguments (3 numbers or 1 array or 1 vector)");
    break;
  }
  r = gsl_vector_complex_alloc(n);
  switch(n) {
    case 4: gsl_vector_complex_set(r, 0, z0); /* fall through */
    case 3: gsl_vector_complex_set(r, 1, z1); /* fall through */
    case 2: gsl_vector_complex_set(r, 2, z2); /* fall through */
    case 1: gsl_vector_complex_set(r, 3, z3);
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, r);
}

#endif
/* singleton method */
VALUE FUNCTION(rb_gsl_poly,complex_solve)(int argc, VALUE *argv, VALUE obj)
{
  int size = -1, size2;
  gsl_poly_complex_workspace *w = NULL;
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector *a = NULL, *z = NULL;
  gsl_complex c;
  gsl_vector_complex *r = NULL;
  int i, flag = 0;
  // local variable "status" declared and set, but never used
  //int status;

  switch (argc) {
  case 1:
    break;
  case 2:
    if (TYPE(argv[1]) == T_FIXNUM) size = FIX2INT(argv[1]);
    break;
  case 3:
    if (TYPE(argv[1]) == T_FIXNUM) size = FIX2INT(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1-3)", argc);
    break;
  }

  switch (TYPE(argv[0])) {
  case T_ARRAY:
    if (size < 0) size = RARRAY_LEN(argv[0]);
    a = gsl_vector_alloc(size);
    for (i = 0; i < size; i++) gsl_vector_set(a, i, NUMCONV2(rb_ary_entry(argv[0], i)));
    break;
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    if (rb_obj_is_kind_of(argv[argc-1], cgsl_poly_workspace)) size = argc - 1;
    else size = argc;
    a = gsl_vector_alloc(size);
    for (i = 0; i < size; i++) gsl_vector_set(a, i, NUMCONV2(argv[i]));
    break;
  default:
    if (rb_obj_is_kind_of(argv[0], GSL_TYPE(cgsl_vector))) {
      Data_Get_Struct(argv[0], GSL_TYPE(gsl_vector), v);
      if (size < 0) size = v->size;
    } else {
      rb_raise(rb_eTypeError, "wrong argument type (Array, Vector, or Numeric expected");
    }
    a = gsl_vector_alloc(v->size);
    for (i = 0; i < size; i++) gsl_vector_set(a, i, FUNCTION(gsl_vector,get)(v, i));
    break;
  }
  size2 = 2*(size - 1);
  z = gsl_vector_alloc(size2);
  if (rb_obj_is_kind_of(argv[argc-1], cgsl_poly_workspace)
      || rb_obj_is_kind_of(argv[argc-1], cgsl_poly_complex_workspace)) {
    Data_Get_Struct(argv[argc-1], gsl_poly_complex_workspace, w);
    flag = 0;
  } else {
    w = gsl_poly_complex_workspace_alloc(size);
    flag = 1;
  }
  /*status =*/ gsl_poly_complex_solve(a->data, size, w, z->data);

  if (flag == 1) gsl_poly_complex_workspace_free(w);
  gsl_vector_free(a);
  r = gsl_vector_complex_alloc(size - 1);
  for (i = 0; i < size - 1; i++) {
    c.dat[0] = gsl_vector_get(z, i*2);
    c.dat[1] = gsl_vector_get(z, i*2+1);
    gsl_vector_complex_set(r, i, c);
  }
  gsl_vector_free(z);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, r);
}


/* a2 x*x + a1 x + a0 = 0 */
static VALUE FUNCTION(rb_gsl_poly,solve_quadratic2)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector *r = NULL;
  gsl_vector_complex *r2 = NULL;
  double a2, a1, a0;
  double d;                /* determinant */
  double x0, x1;
  int n;
  gsl_complex z0, z1;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (v->size < 3) {
    rb_raise(rb_eArgError, "the order of the object is less than 3.");
  }
  a2 = FUNCTION(gsl_vector,get)(v, 2);   /* coefficients */
  a1 = FUNCTION(gsl_vector,get)(v, 1);
  a0 = FUNCTION(gsl_vector,get)(v, 0);
  d = a1*a1 - 4.0*a2*a0;    /* determinant */
  if (d >= 0.0) {
    n = gsl_poly_solve_quadratic(a2, a1, a0, &x0, &x1);
    r = gsl_vector_alloc(n);
    switch(n) {
      case 2: gsl_vector_set(r, 1, x1); /* fall through */
      case 1: gsl_vector_set(r, 0, x0);
    }
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, r);
  } else {
    n = gsl_poly_complex_solve_quadratic(a2, a1, a0, &z0, &z1);
    r2 = gsl_vector_complex_alloc(n);
    switch(n) {
      case 2: gsl_vector_complex_set(r2, 1, z1); /* fall through */
      case 1: gsl_vector_complex_set(r2, 0, z0);
    }
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, r2);
  }

}


static VALUE FUNCTION(rb_gsl_poly,complex_solve_quadratic2)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector_complex *r = NULL;
  double a2, a1, a0;
  int n;
  gsl_complex z0, z1;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (v->size < 3) {
    rb_raise(rb_eArgError, "the order of the object is less than 3.");
  }
  a2 = FUNCTION(gsl_vector,get)(v, 2);   /* coefficients */
  a1 = FUNCTION(gsl_vector,get)(v, 1);
  a0 = FUNCTION(gsl_vector,get)(v, 0);
  n = gsl_poly_complex_solve_quadratic(a2, a1, a0, &z0, &z1);
  r = gsl_vector_complex_alloc(n);
  switch(n) {
    case 2: gsl_vector_complex_set(r, 1, z1); /* fall through */
    case 1: gsl_vector_complex_set(r, 0, z0);
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, r);
}

/* x**3 + a2 x**2 + a1 x + a0 = 0 */
static VALUE FUNCTION(rb_gsl_poly,solve_cubic2)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector *r = NULL;
  double a3, a2, a1, a0;
  double x0, x1, x2;
  int n;

  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (v->size < 4) {
    rb_raise(rb_eArgError, "the order of the object is less than 4.");
  }
  a3 = FUNCTION(gsl_vector,get)(v, 3);   
  a2 = FUNCTION(gsl_vector,get)(v, 2)/a3;   /* coefficients */
  a1 = FUNCTION(gsl_vector,get)(v, 1)/a3;
  a0 = FUNCTION(gsl_vector,get)(v, 0)/a3;
  n = gsl_poly_solve_cubic(a2, a1, a0, &x0, &x1, &x2);
  r = gsl_vector_alloc(n);
  switch(n) {
    case 3: gsl_vector_set(r, 2, x2); /* fall through */
    case 2: gsl_vector_set(r, 1, x1); /* fall through */
    case 1: gsl_vector_set(r, 0, x0);
  }
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, r);
}

static VALUE FUNCTION(rb_gsl_poly,complex_solve_cubic2)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector_complex *r = NULL;
  double a3, a2, a1, a0;
  int n;
  gsl_complex z0, z1, z2;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (v->size < 4) {
    rb_raise(rb_eArgError, "the order of the object is less than 4.");
  }
  a3 = FUNCTION(gsl_vector,get)(v, 3);   
  a2 = FUNCTION(gsl_vector,get)(v, 2)/a3;   /* coefficients */
  a1 = FUNCTION(gsl_vector,get)(v, 1)/a3;
  a0 = FUNCTION(gsl_vector,get)(v, 0)/a3;
  n = gsl_poly_complex_solve_cubic(a2, a1, a0, &z0, &z1, &z2);
  r = gsl_vector_complex_alloc(n);
  switch(n) {
    case 3: gsl_vector_complex_set(r, 2, z2); /* fall through */
    case 2: gsl_vector_complex_set(r, 1, z1); /* fall through */
    case 1: gsl_vector_complex_set(r, 0, z0);
  }
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, r);
}

#ifdef HAVE_GSL_POLY_SOLVE_QUARTIC
/* a4 x**4 + a3 x**3 + a2 x**2 + a1 x + a0 = 0 */
static VALUE FUNCTION(rb_gsl_poly,solve_quartic2)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector *r = NULL;
  double a4, a3, a2, a1, a0;
  double x0, x1, x2, x3;
  int n;

  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (v->size < 5) {
    rb_raise(rb_eArgError, "the order of the object is less than 4.");
  }
  a4 = FUNCTION(gsl_vector,get)(v, 4);   
  a3 = FUNCTION(gsl_vector,get)(v, 3)/a4;   
  a2 = FUNCTION(gsl_vector,get)(v, 2)/a4;   /* coefficients */
  a1 = FUNCTION(gsl_vector,get)(v, 1)/a4;
  a0 = FUNCTION(gsl_vector,get)(v, 0)/a4;
  n = gsl_poly_solve_quartic(a3, a2, a1, a0, &x0, &x1, &x2, &x3);
  r = gsl_vector_alloc(3);
  r = gsl_vector_alloc(4);
  gsl_vector_set(r, 0, x0);
  gsl_vector_set(r, 1, x1);
  gsl_vector_set(r, 2, x2);
  gsl_vector_set(r, 3, x3);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, r);
}

static VALUE FUNCTION(rb_gsl_poly,complex_solve_quartic2)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector_complex *r = NULL;
  double a4, a3, a2, a1, a0;
  int n;
  gsl_complex z0, z1, z2, z3;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  if (v->size < 5) {
    rb_raise(rb_eArgError, "the order of the object is less than 4.");
  }
  a4 = FUNCTION(gsl_vector,get)(v, 4);   
  a3 = FUNCTION(gsl_vector,get)(v, 3)/a4;   
  a2 = FUNCTION(gsl_vector,get)(v, 2)/a4;   /* coefficients */
  a1 = FUNCTION(gsl_vector,get)(v, 1)/a4;
  a0 = FUNCTION(gsl_vector,get)(v, 0)/a4;
  n = gsl_poly_complex_solve_quartic(a3, a2, a1, a0, &z0, &z1, &z2, &z3);
  r = gsl_vector_complex_alloc(4);
  gsl_vector_complex_set(r, 0, z0);
  gsl_vector_complex_set(r, 1, z1);
  gsl_vector_complex_set(r, 2, z2);
  gsl_vector_complex_set(r, 3, z3);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, r);
}

#endif

/* method */
VALUE FUNCTION(rb_gsl_poly,complex_solve2)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  gsl_vector *a = NULL, *z = NULL;
  size_t size, size2, i;
  gsl_poly_complex_workspace *w = NULL;
  gsl_complex c;
  gsl_vector_complex *r = NULL;
  int flag = 0;
  // local variable "status" declared and set, but never used
  //int status;

  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  size = v->size;
  size2 = 2*(size - 1);
  z = gsl_vector_alloc(size2);
  a = gsl_vector_alloc(v->size);
  for (i = 0; i < v->size; i++) {
    gsl_vector_set(a, i, FUNCTION(gsl_vector,get)(v, i));
  }

  if (argc == 1 && rb_obj_is_kind_of(argv[0], cgsl_poly_workspace)) {
    Data_Get_Struct(argv[0], gsl_poly_complex_workspace, w);
    flag = 0;
  } else {
    w = gsl_poly_complex_workspace_alloc(size);
    flag = 1;
  }
  
  /*status =*/ gsl_poly_complex_solve(a->data, size, w, z->data);

  r = gsl_vector_complex_alloc(size - 1);
  for (i = 0; i < size - 1; i++) {
    c.dat[0] = gsl_vector_get(z, i*2);
    c.dat[1] = gsl_vector_get(z, i*2+1);
    gsl_vector_complex_set(r, i, c);
  }
  gsl_vector_free(a);
  gsl_vector_free(z);
  if (flag == 1) gsl_poly_complex_workspace_free(w);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, r);
}

#ifdef BASE_INT
static VALUE rb_gsl_poly_int_to_f(VALUE obj)
{
  gsl_vector *v;
  gsl_vector_int *vi;
  size_t i;
  Data_Get_Struct(obj, gsl_vector_int, vi);
  v = gsl_vector_alloc(vi->size);
  for (i = 0; i < v->size; i++) 
    gsl_vector_set(v, i, (double) gsl_vector_int_get(vi, i));
  return Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, v);
}
#endif

#ifdef BASE_DOUBLE
static VALUE rb_gsl_poly_to_i(VALUE obj)
{
  gsl_vector *v;
  gsl_vector_int *vi;
  size_t i;
  Data_Get_Struct(obj, gsl_vector, v);
  vi = gsl_vector_int_alloc(v->size);
  for (i = 0; i < v->size; i++) 
    gsl_vector_int_set(vi, i, (int) gsl_vector_get(v, i));
  return Data_Wrap_Struct(cgsl_poly_int, 0, gsl_vector_int_free, vi);
}

static VALUE FUNCTION(rb_gsl_poly,workspace_new)(VALUE klass, VALUE n)
{
  gsl_poly_complex_workspace *w = NULL;
  w = gsl_poly_complex_workspace_alloc(FIX2INT(n));
  return Data_Wrap_Struct(klass, 0, gsl_poly_complex_workspace_free, w);
}

#ifdef GSL_1_1_LATER

/* singleton method of the class Poly */
static VALUE rb_gsl_poly_dd_init(VALUE obj, VALUE vxa, VALUE vya)
{
  gsl_vector *xa = NULL, *ya = NULL;
  gsl_poly *dd = NULL;
  Data_Get_Struct(vxa, gsl_vector, xa);
  Data_Get_Struct(vya, gsl_vector, ya);
  dd = gsl_vector_alloc(xa->size);
  gsl_poly_dd_init(dd->data, xa->data, ya->data, xa->size);
  return Data_Wrap_Struct(cgsl_poly_dd, 0, gsl_vector_free, dd);
}

static VALUE rb_gsl_poly_dd_eval(VALUE obj, VALUE xxa, VALUE xx)
{
  gsl_vector *v = NULL;
  gsl_matrix *m = NULL;
  gsl_vector *dd = NULL, *xa, *vnew = NULL;
  gsl_matrix *mnew = NULL;
  VALUE x, ary;
  size_t size, i, j;
  Data_Get_Struct(obj, gsl_vector, dd);
  CHECK_VECTOR(xxa);
  Data_Get_Struct(xxa, gsl_vector, xa);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new(gsl_poly_dd_eval(dd->data, xa->data, dd->size,
					 NUM2DBL(xx)));
    break;
  case T_ARRAY:
    size = RARRAY_LEN(xx);
    ary = rb_ary_new2(size);
    for (i = 0; i < size; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(x);
      rb_ary_store(ary, i, 
		   rb_float_new(gsl_poly_dd_eval(dd->data, xa->data, 
						 dd->size, NUM2DBL(x))));
    }
    return ary;
    break;
  default:
    if (VEC_P(xx)) {
      Data_Get_Struct(xx, GSL_TYPE(gsl_vector), v);
      size = v->size;
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < size; i++) {
	gsl_vector_set(vnew, i, 
		       gsl_poly_dd_eval(dd->data, xa->data, 
					dd->size, FUNCTION(gsl_vector,get)(v, i)));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MAT_P(xx)) {
      Data_Get_Struct(xx, GSL_TYPE(gsl_matrix), m);
      size = m->size1;
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  gsl_matrix_set(mnew, i, j, 
			 gsl_poly_dd_eval(dd->data, xa->data, 
					  dd->size, gsl_matrix_get(m, i, j)));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type");
    }
    break;
  }
  return Qnil; /* never reach here */
}

static VALUE rb_gsl_poly_dd_taylor(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *dd = NULL;
  gsl_vector *xa, *c = NULL, *w = NULL;
  double xp;
  size_t size;
  int flag = 0;
  Data_Get_Struct(obj, gsl_vector, dd);
  switch (argc) {
  case 2:
    size = dd->size;
    xp = NUM2DBL(argv[0]);
    CHECK_VECTOR(argv[1]);
    Data_Get_Struct(argv[1], gsl_vector, xa);
    w = gsl_vector_alloc(size);
    flag = 1;
    break;
  case 3:
    xp = NUM2DBL(argv[0]);
    CHECK_VECTOR(argv[1]);
    Data_Get_Struct(argv[1], gsl_vector, xa);
    if (TYPE(argv[2]) == T_FIXNUM) {
      size = FIX2INT(argv[2]);
      w = gsl_vector_alloc(size);
      flag = 1;
    } else {
      CHECK_VEC(argv[2]);
      Data_Get_Struct(argv[2], GSL_TYPE(gsl_vector), w);
      size = dd->size;
    }
    break;
  case 4:
    Need_Float(argv[0]);
    CHECK_VECTOR(argv[1]);
    CHECK_FIXNUM(argv[2]);
    CHECK_VECTOR(argv[3]);
    xp = NUM2DBL(argv[0]);
    Data_Get_Struct(argv[1], gsl_vector, xa);
    size = FIX2INT(argv[2]);
    Data_Get_Struct(argv[3], GSL_TYPE(gsl_vector), w);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
  }
  c = gsl_vector_alloc(size);
  gsl_poly_dd_taylor(c->data, xp, dd->data, xa->data, size, w->data);
  if (flag == 1) gsl_vector_free(w);
  return Data_Wrap_Struct(cgsl_poly_taylor, 0, gsl_vector_free, c);
}

#endif
#endif

static VALUE FUNCTION(rb_gsl_poly,order)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  return INT2FIX(v->size - 1);
}

int FUNCTION(gsl_poly,conv)(const BASE *a, size_t na, const BASE *b, size_t nb,
		  BASE *c, size_t *nc)
{
  BASE x;
  size_t i, j;
  *nc = na + nb - 1;
  for (i = 0; i < *nc; i++) c[i] = 0;
  for (i = 0; i < *nc; i++) {
    if (i >= na) break;
    x = a[i];
    for (j = 0; j < *nc; j++) {
      if (j >= nb) break;
      else c[i+j] += x*b[j];
    }
  }
  return 0;
}

GSL_TYPE(gsl_vector)* FUNCTION(gsl_poly,conv_vector)(const GSL_TYPE(gsl_vector) *v1, const GSL_TYPE(gsl_vector) *v2)
{
  GSL_TYPE(gsl_vector) *vnew = NULL;
  size_t n, tmp;
  if (v1->size == 1) { 
    vnew = FUNCTION(make_vector,clone)(v2);
    FUNCTION(gsl_vector,scale)(vnew, FUNCTION(gsl_vector,get)(v1, 0));
    return vnew;
  } else if (v2->size == 1) {
    vnew = FUNCTION(make_vector,clone)(v1);
    FUNCTION(gsl_vector,scale)(vnew, FUNCTION(gsl_vector,get)(v2, 0));
    return vnew;
  } else {
    n = v1->size + v2->size - 1;
    vnew = FUNCTION(gsl_vector,calloc)(n);
    FUNCTION(gsl_poly,conv)(v1->data, v1->size, v2->data, v2->size, vnew->data, &tmp);
    return vnew;
  }
}

GSL_TYPE(gsl_vector)* FUNCTION(gsl_poly,reduce)(const GSL_TYPE(gsl_vector) *v)
{
  size_t i, nn = v->size;
  GSL_TYPE(gsl_vector) *vnew = NULL;
  for (i = v->size-1; 0 <= (int) i; i--) {
    if (!gsl_fcmp(FUNCTION(gsl_vector,get)(v, i), 0.0, 1e-10)) {
      nn = i;
      break;
    }
  }
  if (nn == 0) nn = v->size;
  vnew = FUNCTION(gsl_vector,alloc)(nn);
  for (i = 0; i < nn; i++) {
    FUNCTION(gsl_vector,set)(vnew, i, FUNCTION(gsl_vector,get)(v, i));
  }
  return vnew;
}

GSL_TYPE(gsl_vector)* FUNCTION(gsl_poly,deriv)(const GSL_TYPE(gsl_vector) *v)
{
  GSL_TYPE(gsl_vector) *vnew = NULL;
  size_t i;
  vnew = FUNCTION(gsl_vector,alloc)(v->size - 1);
  for (i = 0; i < v->size - 1; i++) {
    FUNCTION(gsl_vector,set)(vnew, i, FUNCTION(gsl_vector,get)(v, i+1)*(i+1));
  }
  return vnew;
}

GSL_TYPE(gsl_vector)* FUNCTION(gsl_poly,integ)(const GSL_TYPE(gsl_vector) *v)
{
  GSL_TYPE(gsl_vector) *vnew = NULL;
  size_t i;
  vnew = FUNCTION(gsl_vector,alloc)(v->size + 1);
  FUNCTION(gsl_vector,set)(vnew, 0, 0.0);
  for (i = 1; i < v->size + 1; i++) {
    FUNCTION(gsl_vector,set)(vnew, i, FUNCTION(gsl_vector,get)(v, i-1)/i);
  }
  return vnew;
}

GSL_TYPE(gsl_vector)* FUNCTION(gsl_poly,deconv_vector)(const GSL_TYPE(gsl_vector) *c, const GSL_TYPE(gsl_vector) *a, 
				   GSL_TYPE(gsl_vector) **r)
{
  GSL_TYPE(gsl_vector) *vnew = NULL, *a2 = NULL, *c2 = NULL, *vtmp = NULL;
  GSL_TYPE(gsl_vector) *rtmp = NULL;
  BASE x, y, z, aa;
  size_t n, i, j, k, jj;
  c2 = FUNCTION(gsl_poly,reduce)(c);
  a2 = FUNCTION(gsl_poly,reduce)(a);
  n = c2->size - a2->size + 1;
  vnew = FUNCTION(gsl_vector,calloc)(n);
  rtmp = FUNCTION(gsl_vector,alloc)(c2->size - 1);
  aa = FUNCTION(gsl_vector,get)(a2, a2->size - 1);
  FUNCTION(gsl_vector,set)(vnew, n-1, FUNCTION(gsl_vector,get)(c2, c2->size-1)/aa);
  for (i = n - 2, k = 1; k < n; i--, k++) {
    x = FUNCTION(gsl_vector,get)(c2, c2->size-1-k);
    for (j = n-1;; j--) {
      z = FUNCTION(gsl_vector,get)(vnew, j);
      jj = c2->size-1-k-j;
      //if (jj > k || jj < 0) continue;
      if (jj <= k) {
      y = FUNCTION(gsl_vector,get)(a2, jj);
      x -= y*z;
      }
      if (j == 0) break;
    }
    FUNCTION(gsl_vector,set)(vnew, i, x/aa);
  }
  vtmp = FUNCTION(gsl_poly,conv_vector)(vnew, a2);
  for (i = 0; i < rtmp->size; i++) {
    x = FUNCTION(gsl_vector,get)(c2, i);
    y = FUNCTION(gsl_vector,get)(vtmp, i);
    FUNCTION(gsl_vector,set)(rtmp, i, x - y);
  }
  *r = FUNCTION(gsl_poly,reduce)(rtmp);
  FUNCTION(gsl_vector,free)(rtmp);
  FUNCTION(gsl_vector,free)(vtmp);
  FUNCTION(gsl_vector,free)(c2);
  FUNCTION(gsl_vector,free)(a2);
  return vnew;
}

GSL_TYPE(gsl_poly)* FUNCTION(get_poly,get)(VALUE obj, int *flag)
{
  GSL_TYPE(gsl_poly) *p = NULL;
  size_t i;
  switch (TYPE(obj)) {
  case T_ARRAY:
    p = FUNCTION(gsl_vector,alloc)(RARRAY_LEN(obj));
    for (i = 0; i < p->size; i++) FUNCTION(gsl_vector,set)(p, i, (BASE) NUM2DBL(rb_ary_entry(obj, i)));
    *flag = 1;
    break;
  case T_FLOAT:
  case T_FIXNUM:
    p = FUNCTION(gsl_vector,alloc)(1);
    FUNCTION(gsl_vector,set)(p, 0, (BASE) NUM2DBL(obj));
    *flag = 1;
    break;
  default:
    CHECK_VEC(obj);
    Data_Get_Struct(obj, GSL_TYPE(gsl_vector), p);
    *flag = 0;
    break;
  }
  return p;
}

static VALUE FUNCTION(rb_gsl_poly,conv)(VALUE obj, VALUE bb)
{
  GSL_TYPE(gsl_vector) *v = NULL, *v2 = NULL, *vnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  switch (TYPE(bb)) {
  case T_FIXNUM:
  case T_FLOAT:
    vnew = FUNCTION(gsl_vector,alloc)(v->size);
    FUNCTION(gsl_vector,memcpy)(vnew, v);
    FUNCTION(gsl_vector,scale)(vnew, (BASE) NUM2DBL(bb));
    break;
  default:
    CHECK_VEC(bb);
    Data_Get_Struct(bb, GSL_TYPE(gsl_vector), v2);
    vnew = FUNCTION(gsl_poly,conv_vector)(v, v2);
    break;
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), vnew);
}

VALUE FUNCTION(rb_gsl_poly,deconv)(VALUE obj, VALUE bb)
{
  GSL_TYPE(gsl_poly) *v = NULL, *v2 = NULL, *vnew = NULL, *r = NULL;
  int flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  switch (TYPE(bb)) {
  case T_ARRAY:
    v2 = FUNCTION(get_poly,get)(bb, &flag);
    break;
  case T_FLOAT:
  case T_FIXNUM:
    v2 = FUNCTION(gsl_vector,alloc)(1);
    FUNCTION(gsl_vector,set)(v2, 0, (BASE) NUM2DBL(bb));
    break;
  default:
    CHECK_VEC(bb);
    Data_Get_Struct(bb, GSL_TYPE(gsl_vector), v2);
    break;
  }
  vnew = FUNCTION(gsl_poly,deconv_vector)(v, v2, &r);
  if (flag == 1) FUNCTION(gsl_vector,free)(v2);
  if (FUNCTION(gsl_vector,isnull)(r)) 
    return Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), vnew);
  else
    return rb_ary_new3(2, Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), vnew),
		       Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), r));
}

static VALUE FUNCTION(rb_gsl_poly,reduce)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_poly,reduce)(v);
  if (vnew == NULL) {
    return Qnil;
  } else if (vnew->size == 0) {
    return Qnil;
  } else if (FUNCTION(gsl_vector,isnull)(vnew)) {
    return INT2FIX(0);
  } else if (vnew->size == 1) {
    return rb_float_new(FUNCTION(gsl_vector,get)(vnew, 0));
  } else {
    return Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), vnew);
  }
  return Qnil; /* never reach here */
}

static VALUE FUNCTION(rb_gsl_poly,deriv)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_poly,deriv)(v);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), vnew);
}

static VALUE FUNCTION(rb_gsl_poly,integ)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(gsl_poly,integ)(v);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), vnew);
}

static VALUE FUNCTION(rb_gsl_poly,conv2)(VALUE klass, VALUE v1, VALUE v2)
{
  GSL_TYPE(gsl_poly) *p1 = NULL, *p2 = NULL, *p3 = NULL;
  int flag1 = 0, flag2 = 0;
  size_t i;
  VALUE ary;
  p1 = FUNCTION(get_poly,get)(v1, &flag1);
  p2 = FUNCTION(get_poly,get)(v2, &flag2);
  p3 = FUNCTION(gsl_poly,conv_vector)(p1, p2);
  if (flag1 == 1) FUNCTION(gsl_vector,free)(p1);
  if (flag2 == 1) FUNCTION(gsl_vector,free)(p2);
  if (flag1 == 1 && flag2 == 1) {	
    ary = rb_ary_new2(p3->size);
    for (i = 0; i < p3->size; i++)
      rb_ary_store(ary, i, C_TO_VALUE2(FUNCTION(gsl_vector,get)(p3, i)));
    FUNCTION(gsl_vector,free)(p3);
    return ary;
  } else {
    return Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), p3);
  }
}

static VALUE FUNCTION(rb_gsl_poly,deconv2)(VALUE klass, VALUE v1, VALUE v2)
{
  GSL_TYPE(gsl_poly) *p1 = NULL, *p2 = NULL;
  GSL_TYPE(gsl_poly) *r = NULL, *vnew = NULL;
  int flag1 = 0, flag2 = 0;
  p1 = FUNCTION(get_poly,get)(v1, &flag1);
  p2 = FUNCTION(get_poly,get)(v2, &flag2);
  vnew = FUNCTION(gsl_poly,deconv_vector)(p1, p2, &r);
  if (flag1 == 1) FUNCTION(gsl_vector,free)(p1);
  if (flag2 == 1) FUNCTION(gsl_vector,free)(p2);
  if (FUNCTION(gsl_vector,isnull)(r)) 
    return Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), vnew);
  else
    return rb_ary_new3(2, Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), vnew),
		       Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), r));
}

GSL_TYPE(gsl_poly)* FUNCTION(gsl_poly,add)(const GSL_TYPE(gsl_poly) *a, 
	const GSL_TYPE(gsl_poly) *b)
{
  GSL_TYPE(gsl_poly) *c = NULL;
  const GSL_TYPE(gsl_poly) *longer;
  size_t i, n;
  if (a->size > b->size) {
    c = FUNCTION(gsl_vector,alloc)(a->size);
    longer = a;
  } else {
    c = FUNCTION(gsl_vector,alloc)(b->size);
      longer = b;
  }
  n = GSL_MIN(a->size, b->size);
  for (i = 0; i < n; i++) {
    FUNCTION(gsl_vector,set)(c, i, FUNCTION(gsl_vector,get)(a, i) + FUNCTION(gsl_vector,get)(b, i));
  }
  for (i = n; i < c->size; i++)
    FUNCTION(gsl_vector,set)(c, i, FUNCTION(gsl_vector,get)(longer, i));
  return c;
}

static VALUE FUNCTION(rb_gsl_poly,add)(VALUE obj, VALUE bb)
{
  GSL_TYPE(gsl_vector) *v = NULL, *vnew = NULL, *vb = NULL;
  BASE b;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  switch (TYPE(bb)) {
  case T_FLOAT:
  case T_FIXNUM:
    b = (BASE) NUM2DBL(bb);
    vnew = FUNCTION(gsl_vector,alloc)(v->size);
    FUNCTION(gsl_vector,memcpy)(vnew, v);
    FUNCTION(gsl_vector,set)(vnew, 0, FUNCTION(gsl_vector,get)(v, 0) + b);
    break;
  default:
    CHECK_VEC(bb);
    Data_Get_Struct(bb, GSL_TYPE(gsl_vector), vb);
    vnew = FUNCTION(gsl_poly,add)(v, vb);
  }
  return Data_Wrap_Struct(CLASS_OF(obj), 0, FUNCTION(gsl_vector,free), vnew);
}

static VALUE rb_gsl_poly_uminus(VALUE obj);
static VALUE rb_gsl_poly_int_uminus(VALUE obj);
static VALUE FUNCTION(rb_gsl_poly,sub)(VALUE obj, VALUE bb)
{
  switch (TYPE(bb)) {
  case T_FLOAT:
  case T_FIXNUM:
    return FUNCTION(rb_gsl_poly,add)(obj, C_TO_VALUE2(-(BASE)NUM2DBL(bb)));
    break;
  default:
    CHECK_VEC(bb);
    return FUNCTION(rb_gsl_poly,add)(obj, FUNCTION(rb_gsl_poly,uminus)(bb));
    break;
  }
}

static VALUE FUNCTION(rb_gsl_poly,uminus)(VALUE obj)
{
  GSL_TYPE(gsl_poly) *p = NULL, *pnew = NULL;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), p);
  pnew = FUNCTION(gsl_vector,alloc)(p->size);
  for (i = 0; i < pnew->size; i++) FUNCTION(gsl_vector,set)(pnew, i, -FUNCTION(gsl_vector,get)(p, i));
  return Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), pnew);
}

static VALUE FUNCTION(rb_gsl_poly,uplus)(VALUE obj)
{
  return obj;
}

static VALUE FUNCTION(rb_gsl_poly,coerce)(VALUE obj, VALUE other)
{
  GSL_TYPE(gsl_vector) *vb;
  switch (TYPE(other)) {
  case T_FLOAT:
  case T_FIXNUM:
    vb = FUNCTION(gsl_vector,calloc)(1);
    FUNCTION(gsl_vector,set)(vb, 0, (BASE) NUM2DBL(other));
    return rb_ary_new3(2, Data_Wrap_Struct(CLASS_OF(obj), 0, FUNCTION(gsl_vector,free), vb),
		       obj);
    break;
  default:
    CHECK_VEC(other);
    return rb_ary_new3(3, other, obj);
    break;
  }
}

static VALUE FUNCTION(rb_gsl_poly,to_gv)(VALUE obj)
{
  GSL_TYPE(gsl_vector) *v, *vnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(gsl_vector), v);
  vnew = FUNCTION(make_vector,clone)(v);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_poly), 0, FUNCTION(gsl_vector,free), vnew);
}

static VALUE FUNCTION(rb_gsl_poly,companion_matrix)(VALUE obj)
{
  GSL_TYPE(gsl_poly) *p = NULL;
  BASE z;
  gsl_matrix *m;
  size_t i, j, size;
  Data_Get_Struct(obj, GSL_TYPE(gsl_poly), p);
  size = p->size - 1;
  m = gsl_matrix_calloc(size, size);
  z = FUNCTION(gsl_vector,get)(p, size);
  for (j = 0; j < size; j++) 
    gsl_matrix_set(m, 0, size-j-1, -FUNCTION(gsl_vector,get)(p, j)/z);
  for (i = 1; i < size; i++) {
   gsl_matrix_set(m, i, i-1, 1.0);
  }
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, m);
}

static VALUE FUNCTION(rb_gsl_poly,info)(VALUE obj)
{
  GSL_TYPE(gsl_poly) *v;
  char buf[256];
  Data_Get_Struct(obj, GSL_TYPE(gsl_poly), v);
  sprintf(buf, "Class:      %s\n", rb_class2name(CLASS_OF(obj)));
  sprintf(buf, "%sSuperClass: %s\n", buf, rb_class2name(RCLASS_SUPER(CLASS_OF(obj))));
  sprintf(buf, "%sOrder:      %d\n", buf, (int) v->size-1);
  return rb_str_new2(buf);
}

#ifdef BASE_DOUBLE
#include "include/rb_gsl_fit.h"
/* singleton */
static VALUE rb_gsl_poly_fit(int argc, VALUE *argv, VALUE obj)
{
  gsl_multifit_linear_workspace *space = NULL;
  gsl_matrix *X = NULL, *cov = NULL;
  gsl_vector *x, *y = NULL, *c = NULL;
  gsl_vector_view xx, yy;
  size_t order, i, j;
  double chisq, val;
  int status, flag = 0;
  VALUE vc, vcov;
  if (argc != 3 && argc != 4) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 3 or 4)", argc);
  x = &xx.vector;
  y = &yy.vector;
  Data_Get_Vector(argv[0], x);
  Data_Get_Vector(argv[1], y);
  order = NUM2INT(argv[2]);
  if (argc == 4) {
    Data_Get_Struct(argv[3], gsl_multifit_linear_workspace, space);
  } else {
    space = gsl_multifit_linear_alloc(x->size, order + 1);
    flag = 1;
  }
  cov = gsl_matrix_alloc(order + 1, order + 1);
  c = gsl_vector_alloc(order + 1);
  X = gsl_matrix_alloc(x->size, order + 1);
  for (i = 0; i < x->size; i++) {
    val = 1.0;
    gsl_matrix_set(X, i, 0, val);
    for (j = 1; j <= order; j++) {
      val *= gsl_vector_get(x, i);
      gsl_matrix_set(X, i, j, val);
    }
  }
  status = gsl_multifit_linear(X, y, c, cov, &chisq, space);
  if (flag == 1) gsl_multifit_linear_free(space);
  vc = Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, c);
  vcov = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, cov);
  gsl_matrix_free(X);
  return rb_ary_new3(4, vc, vcov, rb_float_new(chisq), INT2FIX(status));
}

static VALUE rb_gsl_poly_wfit(int argc, VALUE *argv, VALUE obj)
{
  gsl_multifit_linear_workspace *space = NULL;
  gsl_matrix *X = NULL, *cov = NULL;
  gsl_vector *x, *y = NULL, *w, *c = NULL;
  size_t order, i, j;
  double chisq, val;
  int status, flag = 0;
  VALUE vc, vcov;
  if (argc != 4 && argc != 5) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 4 or 5)", argc);
  Data_Get_Vector(argv[0], x);
  Data_Get_Vector(argv[1], w);
  Data_Get_Vector(argv[2], y);
  order = NUM2INT(argv[3]);
  if (argc == 5) {
    Data_Get_Struct(argv[4], gsl_multifit_linear_workspace, space);
  } else {
    space = gsl_multifit_linear_alloc(x->size, order + 1);
    flag = 1;
  }
  cov = gsl_matrix_alloc(order + 1, order + 1);
  c = gsl_vector_alloc(order + 1);
  X = gsl_matrix_alloc(x->size, order + 1);
  for (i = 0; i < x->size; i++) {
    val = 1.0;
    gsl_matrix_set(X, i, 0, val);
    for (j = 1; j <= order; j++) {
      val *= gsl_vector_get(x, i);
      gsl_matrix_set(X, i, j, val);
    }
  }
  status = gsl_multifit_wlinear(X, w, y, c, cov, &chisq, space);
  if (flag == 1) gsl_multifit_linear_free(space);
  vc = Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, c);
  vcov = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, cov);
  gsl_matrix_free(X);
  return rb_ary_new3(4, vc, vcov, rb_float_new(chisq), INT2FIX(status));
}
#endif

#ifdef BASE_DOUBLE
#ifdef GSL_1_13_LATER
static VALUE rb_gsl_poly_eval_derivs_singleton(int argc, VALUE *argv, VALUE klass)
{
  VALUE ary;
  gsl_vector *v = NULL, *v2 = NULL;
  size_t i, lenc, lenres;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2;
  int shape[1];
#endif

  if (argc < 2) rb_raise(rb_eArgError, "Wrong number of arguments (%d for >= 2)", argc);
  if (rb_obj_is_kind_of(argv[0], rb_cArray)) {
    v = gsl_vector_alloc(RARRAY_LEN(argv[0]));
    lenc = v->size;
    for (i = 0; i < lenc; i++) {
      gsl_vector_set(v, i, NUM2DBL(rb_ary_entry(argv[0], i)));
    }
    if (argc == 2) lenres = lenc + 1;
    else lenres = FIX2INT(argv[2]);
    v2 = gsl_vector_alloc(lenres);
    gsl_poly_eval_derivs(v->data, lenc, NUM2DBL(argv[1]), v2->data, lenres);
    ary = rb_ary_new2(lenres);
    for (i = 0; i < lenres; i++) {
      rb_ary_store(ary, i, rb_float_new(gsl_vector_get(v2, i)));
    }
    gsl_vector_free(v2);
    gsl_vector_free(v);
    return ary;
  }
  if (rb_obj_is_kind_of(argv[0], cgsl_vector)) {
    Data_Get_Struct(argv[0], gsl_vector, v);
    lenc = v->size;
    if (argc == 2) lenres = lenc + 1;
    else lenres = FIX2INT(argv[2]);
    v2 = gsl_vector_alloc(lenres);
    gsl_poly_eval_derivs(v->data, lenc, NUM2DBL(argv[1]), v2->data, lenres);
    return Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, v2);
  }
#ifdef HAVE_NARRAY_H
  if (NA_IsNArray(argv[0])) {
    GetNArray(argv[0], na);
    ptr1 = (double*) na->ptr;
    lenc = na->total;
    if (argc == 2) lenres = lenc + 1;
    else lenres = FIX2INT(argv[2]);
    shape[0] = lenres;
    ary = na_make_object(NA_DFLOAT, na->rank, shape, CLASS_OF(argv[0]));
    ptr2 = NA_PTR_TYPE(ary,double*);
    gsl_poly_eval_derivs(ptr1, lenc, NUM2DBL(argv[1]), ptr2, lenres);
    return ary;
  }
#endif
  return Qnil;  // Never comes here
}
static VALUE rb_gsl_poly_eval_derivs(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *v, *v2;
  size_t lenc, lenres;
  Data_Get_Struct(obj, gsl_vector, v);
  lenc = v->size;
  switch (argc) {
  case 1:
    lenres = lenc + 1;
    break;
  case 2:
    lenres = FIX2INT(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for > 1)", argc);
  }
  v2 = gsl_vector_alloc(lenres);
  gsl_poly_eval_derivs(v->data, lenc, NUM2DBL(argv[0]), v2->data, lenres);
  return Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, v2);
}
#endif
#endif

void FUNCTION(Init_gsl_poly,init)(VALUE module)
{
#ifdef BASE_DOUBLE
  VALUE mgsl_poly_complex;
  cgsl_poly = rb_define_class_under(module, "Poly", cgsl_vector);
  cgsl_poly_int = rb_define_class_under(cgsl_poly, "Int", cgsl_vector_int);

  cgsl_poly_workspace = rb_define_class_under(cgsl_poly, "Workspace", cGSL_Object);
  mgsl_poly_complex = rb_define_module_under(cgsl_poly, "Complex");
  cgsl_poly_complex_workspace = rb_define_class_under(mgsl_poly_complex, 
						      "Workspace", cGSL_Object);
  rb_define_singleton_method(cgsl_poly_workspace, "alloc", 
			     rb_gsl_poly_workspace_new, 1);
  rb_define_singleton_method(cgsl_poly_complex_workspace, "alloc",
			     rb_gsl_poly_workspace_new, 1);

  rb_define_singleton_method(mgsl_poly_complex, "solve_quadratic", 
			     FUNCTION(rb_gsl_poly,complex_solve_quadratic), -1);
  rb_define_singleton_method(mgsl_poly_complex, "solve_cubic",
			     FUNCTION(rb_gsl_poly,complex_solve_cubic), -1);
#ifdef HAVE_GSL_POLY_SOLVE_QUARTIC
  rb_define_singleton_method(mgsl_poly_complex, "solve_quartic",
			     FUNCTION(rb_gsl_poly,complex_solve_quartic), -1);
#endif

  rb_define_singleton_method(mgsl_poly_complex, "solve",
			     FUNCTION(rb_gsl_poly,complex_solve), -1);
  rb_define_singleton_method(mgsl_poly_complex, "roots", 
			     FUNCTION(rb_gsl_poly,complex_solve), -1);
#endif

  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "solve_quadratic", 
			     FUNCTION(rb_gsl_poly,solve_quadratic), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "solve_cubic", 
			     FUNCTION(rb_gsl_poly,solve_cubic), -1);

  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "complex_solve_quadratic", 
			     FUNCTION(rb_gsl_poly,complex_solve_quadratic), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "complex_solve_cubic", 
			     FUNCTION(rb_gsl_poly,complex_solve_cubic), -1);
#ifdef HAVE_GSL_POLY_SOLVE_QUARTIC
  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "solve_quartic",
			     FUNCTION(rb_gsl_poly,solve_quartic), -1);

  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "complex_solve_quartic",
			     FUNCTION(rb_gsl_poly,complex_solve_quartic), -1);
#endif
  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "complex_solve",
			     FUNCTION(rb_gsl_poly,complex_solve), -1);

  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "solve",
			     FUNCTION(rb_gsl_poly,complex_solve), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "roots", 
			     FUNCTION(rb_gsl_poly,complex_solve), -1);

  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "eval", 
			     FUNCTION(rb_gsl_poly,eval2), -1);

  rb_define_method(GSL_TYPE(cgsl_poly), "eval",
		   FUNCTION(rb_gsl_poly,eval), 1);
  rb_define_alias(GSL_TYPE(cgsl_poly), "at", "eval");

  rb_define_method(GSL_TYPE(cgsl_poly), "solve_quadratic", 
		   FUNCTION(rb_gsl_poly,solve_quadratic2), 0);
  rb_define_method(GSL_TYPE(cgsl_poly), "complex_solve_quadratic", 
		   FUNCTION(rb_gsl_poly,complex_solve_quadratic2), 0);

  rb_define_method(GSL_TYPE(cgsl_poly), "solve_cubic", 
		   FUNCTION(rb_gsl_poly,solve_cubic2), 0);
  rb_define_method(GSL_TYPE(cgsl_poly), "complex_solve_cubic", 
		   FUNCTION(rb_gsl_poly,complex_solve_cubic2), 0);

#ifdef HAVE_GSL_POLY_SOLVE_QUARTIC
  rb_define_method(GSL_TYPE(cgsl_poly), "solve_quartic", 
		   FUNCTION(rb_gsl_poly,solve_quartic2), 0);
  rb_define_method(GSL_TYPE(cgsl_poly), "complex_solve_quartic", 
		   FUNCTION(rb_gsl_poly,complex_solve_quartic2), 0);
#endif

  rb_define_method(GSL_TYPE(cgsl_poly), "complex_solve", 
		   FUNCTION(rb_gsl_poly,complex_solve2), -1);
  rb_define_alias(GSL_TYPE(cgsl_poly), "solve", "complex_solve");
  rb_define_alias(GSL_TYPE(cgsl_poly), "roots", "complex_solve");

#ifdef BASE_INT
  rb_define_method(cgsl_poly_int, "to_f", rb_gsl_poly_int_to_f, 0);
#endif

#ifdef BASE_DOUBLE
  //  rb_define_singleton_method(cgsl_poly, "eval", rb_gsl_poly_eval_singleton, 2);
  rb_define_method(cgsl_poly, "to_i", rb_gsl_poly_to_i, 0);
#ifdef GSL_1_11_LATER
  rb_define_singleton_method(cgsl_poly, "complex_eval", rb_gsl_poly_eval_singleton, 2);
  rb_define_method(cgsl_vector_complex, "eval", rb_gsl_complex_poly_complex_eval, 1);
#endif
#ifdef GSL_1_1_LATER
  cgsl_poly_dd = rb_define_class_under(cgsl_poly, "DividedDifference", cgsl_poly);
  cgsl_poly_taylor = rb_define_class_under(cgsl_poly, "Taylor", cgsl_poly);
  rb_define_singleton_method(cgsl_poly, "dd_init", rb_gsl_poly_dd_init, 2);

  rb_define_method(cgsl_poly_dd, "eval",rb_gsl_poly_dd_eval, 2);
  rb_define_method(cgsl_poly_dd, "taylor", rb_gsl_poly_dd_taylor, -1);
#endif
#endif

  rb_define_method(GSL_TYPE(cgsl_poly), "order", FUNCTION(rb_gsl_poly,order), 0);

/*****/
  rb_define_method(GSL_TYPE(cgsl_poly), "conv", FUNCTION(rb_gsl_poly,conv), 1);
  rb_define_alias(GSL_TYPE(cgsl_poly), "*", "conv");
  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "conv",
			     FUNCTION(rb_gsl_poly,conv2), 2);

  rb_define_method(GSL_TYPE(cgsl_poly), "deconv",
		   FUNCTION(rb_gsl_poly,deconv), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "deconv",
			     FUNCTION(rb_gsl_poly,deconv2), 2);

  rb_define_method(GSL_TYPE(cgsl_poly), "reduce", 
		   FUNCTION(rb_gsl_poly,reduce), 1);
  rb_define_method(GSL_TYPE(cgsl_poly), "deriv", FUNCTION(rb_gsl_poly,deriv), 1);
  rb_define_method(GSL_TYPE(cgsl_poly), "integ", FUNCTION(rb_gsl_poly,integ), 1);

/*****/

  rb_define_method(GSL_TYPE(cgsl_poly), "add", FUNCTION(rb_gsl_poly,add), 1);
  rb_define_alias(GSL_TYPE(cgsl_poly), "+", "add");
  rb_define_method(GSL_TYPE(cgsl_poly), "sub", FUNCTION(rb_gsl_poly,sub), 1);
  rb_define_alias(GSL_TYPE(cgsl_poly), "-", "sub");

  rb_define_method(GSL_TYPE(cgsl_poly), "-@", FUNCTION(rb_gsl_poly,uminus), 0);
  rb_define_method(GSL_TYPE(cgsl_poly), "+@", FUNCTION(rb_gsl_poly,uplus), 0);

  rb_define_method(GSL_TYPE(cgsl_poly), "coerce", 
		   FUNCTION(rb_gsl_poly,coerce), 1);
  rb_define_method(GSL_TYPE(cgsl_poly), "to_gv", FUNCTION(rb_gsl_poly,to_gv), 0);
  rb_define_alias(GSL_TYPE(cgsl_poly), "to_v", "to_gv");

  rb_define_method(GSL_TYPE(cgsl_poly), "companion_matrix", 
		   FUNCTION(rb_gsl_poly,companion_matrix), 0);
  rb_define_alias(GSL_TYPE(cgsl_poly), "compan", "companion_matrix");

  /*****/
  rb_define_method(GSL_TYPE(cgsl_poly), "info", 
		   FUNCTION(rb_gsl_poly,info), 0);

#ifdef BASE_DOUBLE
  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "fit", 
		   FUNCTION(rb_gsl_poly,fit), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_poly), "wfit", 
		   FUNCTION(rb_gsl_poly,wfit), -1);

#ifdef GSL_1_13_LATER
  rb_define_singleton_method(cgsl_poly, "eval_derivs", rb_gsl_poly_eval_derivs_singleton, -1);
  rb_define_method(cgsl_vector, "eval_derivs", rb_gsl_poly_eval_derivs, -1);
#endif

#endif
}

#undef NUMCONV
#undef NUMCONV2
#undef PRINTF_FORMAT
#undef VEC_ROW_COL
#undef VEC_P
#undef C_TO_VALUE
#undef C_TO_VALUE2
#undef VEC_COL_P
#undef VEC_ROW_P
#undef CHECK_VEC
#undef VEC_VIEW_P

#undef MAT_ROW_COL
#undef MAT_P
#undef C_TO_VALUE
#undef MAT_COL_P
#undef MAT_ROW_P
#undef CHECK_MAT
#undef MAT_VIEW_P
