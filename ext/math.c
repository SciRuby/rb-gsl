/*
  math.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/

#include "rb_gsl_math.h"
#include "rb_gsl_complex.h"
#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif

static void rb_gsl_define_const(VALUE module);

static void rb_gsl_define_const(VALUE module)
{
  rb_define_const(module, "M_E", rb_float_new(M_E));
  rb_define_const(module, "M_LOG2E", rb_float_new(M_LOG2E));
  rb_define_const(module, "M_LOG10E", rb_float_new(M_LOG10E));
  rb_define_const(module, "M_SQRT2", rb_float_new(M_SQRT2));
  rb_define_const(module, "M_SQRT1_2", rb_float_new(M_SQRT1_2));
  rb_define_const(module, "M_SQRT3", rb_float_new(M_SQRT3));
  rb_define_const(module, "M_PI", rb_float_new(M_PI));
  rb_define_const(module, "M_PI_2", rb_float_new(M_PI_2));
  rb_define_const(module, "M_PI_4", rb_float_new(M_PI_4));
  rb_define_const(module, "M_SQRTPI", rb_float_new(M_SQRTPI));
  rb_define_const(module, "M_2_SQRTPI", rb_float_new(M_2_SQRTPI));
  rb_define_const(module, "M_1_PI", rb_float_new(M_1_PI));
  rb_define_const(module, "M_2_PI", rb_float_new(M_2_PI));
  rb_define_const(module, "M_LN10", rb_float_new(M_LN10));
  rb_define_const(module, "M_LN2", rb_float_new(M_LN2));
  rb_define_const(module, "M_LNPI", rb_float_new(M_LNPI));
  rb_define_const(module, "M_EULER", rb_float_new(M_EULER));
  rb_define_const(module, "POSINF", rb_float_new(GSL_POSINF));
  rb_define_const(module, "NEGINF", rb_float_new(GSL_NEGINF));
  rb_define_const(module, "NAN", rb_float_new(GSL_NAN));
}

static VALUE rb_GSL_POSINF(VALUE obj)
{
  return rb_float_new(GSL_POSINF);
}

static VALUE rb_GSL_NEGINF(VALUE obj)
{
  return rb_float_new(GSL_NEGINF);
}

static VALUE rb_GSL_NAN(VALUE obj)
{
  return rb_float_new(GSL_NAN);
}

static VALUE rb_gsl_isnan(VALUE obj, VALUE x)
{
  Need_Float(x);
  return INT2FIX(gsl_isnan(NUM2DBL(x)));
}

static VALUE rb_gsl_isnan2(VALUE obj, VALUE x)
{
  Need_Float(x);
  if (gsl_isnan(NUM2DBL(x))) return Qtrue;
  else return Qfalse;
}

static VALUE rb_gsl_isinf(VALUE obj, VALUE x)
{
  Need_Float(x);
  return INT2FIX(gsl_isinf(NUM2DBL(x)));
}

static VALUE rb_gsl_isinf2(VALUE obj, VALUE x)
{
  //  Need_Float(x);
  if (gsl_isinf(NUM2DBL(x))) return Qtrue;
  else return Qfalse;
}

static VALUE rb_gsl_finite(VALUE obj, VALUE x)
{
  Need_Float(x);
  return INT2FIX(gsl_finite(NUM2DBL(x)));
}

static VALUE rb_gsl_finite2(VALUE obj, VALUE x)
{
  Need_Float(x);
  if (gsl_finite(NUM2DBL(x))) return Qtrue;
  else return Qfalse;
}

/*****/

static VALUE rb_gsl_math_eval(double (*func)(const double), VALUE xx);

VALUE rb_gsl_math_complex_eval(gsl_complex (*func)(gsl_complex), VALUE obj)
{
  gsl_complex *z, *znew;
  gsl_vector_complex *v, *vnew;
  gsl_matrix_complex *m, *mnew;
  size_t i, j;
  if (COMPLEX_P(obj)) {
    Data_Get_Struct(obj, gsl_complex, z);
    znew = xmalloc(sizeof(gsl_complex));
    *znew = (*func)(*z);
    return Data_Wrap_Struct(cgsl_complex, 0, free, znew);
  } else if (VECTOR_COMPLEX_P(obj)) {
    Data_Get_Struct(obj, gsl_vector_complex, v);
    vnew = gsl_vector_complex_alloc(v->size);
    for (i = 0; i < v->size; i++) {
      z = GSL_COMPLEX_AT(v, i);
      gsl_vector_complex_set(vnew, i, (*func)(*z));
    }
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vnew);
  } else if (MATRIX_COMPLEX_P(obj)) {
    Data_Get_Struct(obj, gsl_matrix_complex, m);
    mnew = gsl_matrix_complex_alloc(m->size1, m->size2);
    for (i = 0; i < m->size1; i++) {
      for (j = 0; j < m->size2; j++) {
	gsl_matrix_complex_set(mnew, i, j, (*func)(gsl_matrix_complex_get(m, i, j)));
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

static VALUE rb_gsl_math_eval(double (*func)(const double), VALUE xx)
{
  VALUE x, ary;
  size_t i, size;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2;
#endif
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new((*func)(NUM2DBL(xx)));
    break;
  case T_ARRAY:
    size = RARRAY_LEN(xx);
    ary = rb_ary_new2(size);
    for (i = 0; i < size; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(x);
      //      rb_ary_store(ary, i, rb_float_new((*func)(RFLOAT(x)->value)));
      rb_ary_store(ary, i, rb_float_new((*func)(NUM2DBL(x))));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      GetNArray(xx, na);
      ptr1 = (double*) na->ptr;
      size = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < size; i++) ptr2[i] = (*func)(ptr1[i]);
      return ary;
    }
#endif
    if (VECTOR_P(xx)) {
      return vector_eval_create(xx, func);
    } else if (MATRIX_P(xx)) {
      return matrix_eval_create(xx, func);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (Array or Vector or Matrix expected)", rb_class2name(CLASS_OF(xx)));
    }
    break;
  }
  /* never reach here */
  return Qnil;
}

static VALUE rb_gsl_math_eval2(double (*func)(const double, const double), VALUE xx,
			       VALUE yy);
static VALUE rb_gsl_math_eval2(double (*func)(const double, const double), VALUE xx,
			       VALUE yy)
{
  VALUE x, y, ary;
  size_t i, j, size;
  gsl_vector *v = NULL, *v2 = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *m2 = NULL, *mnew = NULL;
#ifdef HAVE_NARRAY_H
  struct NARRAY *nax, *nay;
  double *ptr1, *ptr2, *ptr3;
#endif
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    Need_Float(yy);
    return rb_float_new((*func)(NUM2DBL(xx), NUM2DBL(yy)));
    break;
  case T_ARRAY:
    Check_Type(yy, T_ARRAY);
    size = RARRAY_LEN(xx);
    //    if (size != RARRAY(yy)->len) rb_raise(rb_eRuntimeError, "array sizes are different.");
    if ((int) size != RARRAY_LEN(yy)) rb_raise(rb_eRuntimeError, "array sizes are different.");
    ary = rb_ary_new2(size);
    for (i = 0; i < size; i++) {
      x = rb_ary_entry(xx, i);
      y = rb_ary_entry(yy, i);
      Need_Float(x); Need_Float(y);
      //      rb_ary_store(ary, i, rb_float_new((*func)(RFLOAT(x)->value, RFLOAT(y)->value)));
      rb_ary_store(ary, i, rb_float_new((*func)(NUM2DBL(x), NUM2DBL(y))));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      GetNArray(xx, nax);
      GetNArray(yy, nay);
      ptr1 = (double*) nax->ptr;
      ptr2 = (double*) nay->ptr;
      size = nax->total;
      ary = na_make_object(NA_DFLOAT, nax->rank, nax->shape, CLASS_OF(xx));
      ptr3 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < size; i++) ptr3[i] = (*func)(ptr1[i], ptr2[i]);
      return ary;
    }
#endif
    if (VECTOR_P(xx)) {
      CHECK_VECTOR(yy);
      Data_Get_Struct(xx, gsl_vector, v);
      Data_Get_Struct(yy, gsl_vector, v2);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	gsl_vector_set(vnew, i, (*func)(gsl_vector_get(v, i), gsl_vector_get(v2, i)));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MATRIX_P(xx)) {
      CHECK_MATRIX(yy);
      Data_Get_Struct(xx, gsl_matrix, m);
      Data_Get_Struct(yy, gsl_matrix, m2);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  gsl_matrix_set(mnew, i, j, (*func)(gsl_matrix_get(m, i, j), gsl_matrix_get(m2, i, j)));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else {
      rb_raise(rb_eTypeError, 
	       "wrong argument type %s "
	       "(Array or Vector or Matrix expected)", rb_class2name(CLASS_OF(xx)));
    }
    break;
  }
  /* never reach here */
  return Qnil;
}

static VALUE rb_gsl_sqrt(VALUE obj, VALUE x)
{
  return rb_gsl_math_eval(sqrt, x);
}

static VALUE rb_gsl_log1p(VALUE obj, VALUE x)
{
  return rb_gsl_math_eval(gsl_log1p, x);
}

static VALUE rb_gsl_expm1(VALUE obj, VALUE x)
{
  return rb_gsl_math_eval(gsl_expm1, x);
}

static VALUE rb_gsl_hypot(VALUE obj, VALUE x, VALUE y)
{
  return rb_gsl_math_eval2(gsl_hypot, x, y);
}
#ifdef GSL_1_10_LATER
static VALUE rb_gsl_hypot3(VALUE obj, VALUE x, VALUE y, VALUE z)
{
  Need_Float(x);
  Need_Float(y);
  Need_Float(z);    
  return rb_float_new(gsl_hypot3(NUM2DBL(x), NUM2DBL(y), NUM2DBL(z)));
}
#endif
static VALUE rb_gsl_acosh(VALUE obj, VALUE x)
{
  if (COMPLEX_P(x) || VECTOR_COMPLEX_P(x) || MATRIX_COMPLEX_P(x)) 
    return rb_gsl_math_complex_eval(gsl_complex_arccosh, x);
  return rb_gsl_math_eval(gsl_acosh, x);
}

static VALUE rb_gsl_asinh(VALUE obj, VALUE x)
{
  if (COMPLEX_P(x) || VECTOR_COMPLEX_P(x) || MATRIX_COMPLEX_P(x)) 
    return rb_gsl_math_complex_eval(gsl_complex_arcsinh, x);
  return rb_gsl_math_eval(gsl_asinh, x);
}

static VALUE rb_gsl_atanh(VALUE obj, VALUE x)
{
  if (COMPLEX_P(x) || VECTOR_COMPLEX_P(x) || MATRIX_COMPLEX_P(x)) 
    return rb_gsl_math_complex_eval(gsl_complex_arctanh, x);
  return rb_gsl_math_eval(gsl_atanh, x);
}

#include <math.h>
/* xx: Numeric, Complex, Vector, Matrix
   nn: Numeric, Complex
*/
VALUE rb_gsl_pow(VALUE obj, VALUE xx, VALUE nn)
{
  VALUE x, ary, argv[2];
  size_t i, j, size;
  double n;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2;
#endif
  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new(pow(NUM2DBL(xx), NUM2DBL(nn)));
    break;
  case T_ARRAY:
    n = NUM2DBL(nn);
    size = RARRAY_LEN(xx);
    ary = rb_ary_new2(size);
    for (i = 0; i < size; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(x);
      rb_ary_store(ary, i, rb_float_new(pow(NUM2DBL(x), n)));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      n = NUM2DBL(nn);
      GetNArray(xx, na);
      ptr1 = (double*) na->ptr;
      size = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < size; i++) ptr2[i] = pow(ptr1[i], n);
      return ary;
    }
#endif
    if (VECTOR_P(xx)) {
      n = NUM2DBL(nn);
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	gsl_vector_set(vnew, i, pow(gsl_vector_get(v, i), n));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } 
    if (MATRIX_P(xx)) {
      n = NUM2DBL(nn);
      Data_Get_Struct(xx, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  gsl_matrix_set(mnew, i, j, pow(gsl_matrix_get(m, i, j), n));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } 
    if (COMPLEX_P(xx) || VECTOR_COMPLEX_P(xx) || MATRIX_COMPLEX_P(xx)) {
      argv[0] = xx;
      argv[1] = nn;
      return rb_gsl_complex_pow(2, argv, obj);
    }
    rb_raise(rb_eTypeError, "wrong argument type %s (Array or Vector or Matrix expected)", rb_class2name(CLASS_OF(xx)));
    break;
  }
  /* never reach here */
  return Qnil;
}

static VALUE rb_gsl_pow_int(VALUE obj, VALUE xx, VALUE nn)
{
  VALUE x, ary, argv[2];
  size_t i, j, size;
  int n;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix *m = NULL, *mnew = NULL;
#ifdef HAVE_NARRAY_H
  struct NARRAY *na;
  double *ptr1, *ptr2;
#endif

  if (CLASS_OF(xx) == rb_cRange) xx = rb_gsl_range2ary(xx);
  switch (TYPE(xx)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return rb_float_new(gsl_pow_int(NUM2DBL(xx), FIX2INT(nn)));
    break;
  case T_ARRAY:
    CHECK_FIXNUM(nn);
    n = FIX2INT(nn);
    size = RARRAY_LEN(xx);
    ary = rb_ary_new2(size);
    for (i = 0; i < size; i++) {
      x = rb_ary_entry(xx, i);
      Need_Float(x);
      //      rb_ary_store(ary, i, rb_float_new(gsl_pow_int(RFLOAT(x)->value, n)));
      rb_ary_store(ary, i, rb_float_new(gsl_pow_int(NUM2DBL(x), n)));
    }
    return ary;
    break;
  default:
#ifdef HAVE_NARRAY_H
    if (NA_IsNArray(xx)) {
      CHECK_FIXNUM(nn);
      n = FIX2INT(nn);
      GetNArray(xx, na);
      ptr1 = (double*) na->ptr;
      size = na->total;
      ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
      ptr2 = NA_PTR_TYPE(ary, double*);
      for (i = 0; i < size; i++) ptr2[i] = gsl_pow_int(ptr1[i], n);
      return ary;
    }
#endif
    if (VECTOR_P(xx)) {
      CHECK_FIXNUM(nn);
      n = FIX2INT(nn);
      Data_Get_Struct(xx, gsl_vector, v);
      vnew = gsl_vector_alloc(v->size);
      for (i = 0; i < v->size; i++) {
	gsl_vector_set(vnew, i, gsl_pow_int(gsl_vector_get(v, i), n));
      }
      return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vnew);
    } else if (MATRIX_P(xx)) {
      CHECK_FIXNUM(nn);
      n = FIX2INT(nn);
      Data_Get_Struct(xx, gsl_matrix, m);
      mnew = gsl_matrix_alloc(m->size1, m->size2);
      for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
	  gsl_matrix_set(mnew, i, j, gsl_pow_int(gsl_matrix_get(m, i, j), n));
	}
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else if (COMPLEX_P(xx) || VECTOR_COMPLEX_P(xx) || MATRIX_COMPLEX_P(xx)) {
      argv[0] = xx;
      argv[1] = nn;
      return rb_gsl_complex_pow_real(2, argv, obj);
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s (Array or Vector or Matrix expected)", rb_class2name(CLASS_OF(xx)));
    }
    break;
  }
  /* never reach here */
  return Qnil;
}

static VALUE rb_gsl_pow_2(VALUE obj, VALUE x)
{
  return rb_gsl_math_eval(gsl_pow_2, x);
}

static VALUE rb_gsl_pow_3(VALUE obj, VALUE x)
{
  return rb_gsl_math_eval(gsl_pow_3, x);
}

static VALUE rb_gsl_pow_4(VALUE obj, VALUE x)
{
  return rb_gsl_math_eval(gsl_pow_4, x);
}

static VALUE rb_gsl_pow_5(VALUE obj, VALUE x)
{
  return rb_gsl_math_eval(gsl_pow_5, x);
}

static VALUE rb_gsl_pow_6(VALUE obj, VALUE x)
{
  return rb_gsl_math_eval(gsl_pow_6, x);
}

static VALUE rb_gsl_pow_7(VALUE obj, VALUE x)
{
  return rb_gsl_math_eval(gsl_pow_7, x);
}

static VALUE rb_gsl_pow_8(VALUE obj, VALUE x)
{
  return rb_gsl_math_eval(gsl_pow_8, x);
}

static VALUE rb_gsl_pow_9(VALUE obj, VALUE x)
{
  return rb_gsl_math_eval(gsl_pow_9, x);
}

/*****/
static VALUE rb_GSL_SIGN(VALUE obj, VALUE x)
{
  return INT2FIX(GSL_SIGN(NUM2DBL(x)));
}

static VALUE rb_GSL_IS_ODD(VALUE obj, VALUE n)
{
  CHECK_FIXNUM(n);
  return INT2FIX(GSL_IS_ODD(FIX2INT(n)));
}

static VALUE rb_GSL_IS_ODD2(VALUE obj, VALUE n)
{
  CHECK_FIXNUM(n);
  if (GSL_IS_ODD(FIX2INT(n))) return Qtrue;
  else return Qfalse;
}

static VALUE rb_GSL_IS_EVEN(VALUE obj, VALUE n)
{
  CHECK_FIXNUM(n);
  return INT2FIX(GSL_IS_EVEN(FIX2INT(n)));
}

static VALUE rb_GSL_IS_EVEN2(VALUE obj, VALUE n)
{
  CHECK_FIXNUM(n);
  if (GSL_IS_EVEN(FIX2INT(n))) return Qtrue;
  else return Qfalse;
}

static VALUE rb_GSL_MAX(VALUE obj, VALUE aa, VALUE bb)
{
  double a, b;
  double max;
  /*  Need_Float(aa); Need_Float(bb);*/
  a = NUM2DBL(aa);
  b = NUM2DBL(bb);
  max = GSL_MAX_DBL(a, b);
  if (gsl_fcmp(max, a, 1.0e-10) == 0) return aa;
  else return bb;

}

static VALUE rb_GSL_MIN(VALUE obj, VALUE aa, VALUE bb)
{
  double a, b;
  double min;
  /*  Need_Float(aa); Need_Float(bb);*/
  a = NUM2DBL(aa);
  b = NUM2DBL(bb);
  min = GSL_MIN_DBL(a, b);
  if (gsl_fcmp(min, a, 1.0e-10) == 0) return aa;
  else return bb;
}

static VALUE rb_GSL_MAX_DBL(VALUE obj, VALUE aa, VALUE bb)
{
  Need_Float(aa); Need_Float(bb);
  return rb_float_new(GSL_MAX_DBL(NUM2DBL(aa), NUM2DBL(bb)));
}

static VALUE rb_GSL_MIN_DBL(VALUE obj, VALUE aa, VALUE bb)
{
  Need_Float(aa); Need_Float(bb);
  return rb_float_new(GSL_MAX_DBL(NUM2DBL(aa), NUM2DBL(bb)));
}

static VALUE rb_GSL_MAX_INT(VALUE obj, VALUE aa, VALUE bb)
{
  if (TYPE(aa) != T_FIXNUM || TYPE(bb) != T_FIXNUM)
    return rb_GSL_MAX(obj, aa, bb);
  else 
    return INT2FIX(GSL_MAX_INT(FIX2INT(aa), FIX2INT(bb)));
}

static VALUE rb_GSL_MIN_INT(VALUE obj, VALUE aa, VALUE bb)
{
  if (TYPE(aa) != T_FIXNUM || TYPE(bb) != T_FIXNUM)
    return rb_GSL_MIN(obj, aa, bb);
  return 
    INT2FIX(GSL_MIN_INT(FIX2INT(aa), FIX2INT(bb)));
}

#ifdef GSL_1_3_LATER
static VALUE rb_gsl_ldexp(VALUE obj, VALUE x, VALUE e)
{
  return rb_float_new(gsl_ldexp(NUM2DBL(x), FIX2INT(e)));
}

static VALUE rb_gsl_frexp(VALUE obj, VALUE x)
{
  int e;
  double val;
  Need_Float(x);
  val = gsl_frexp(NUM2DBL(x), &e);
  return rb_ary_new3(2, rb_float_new(val), INT2FIX(e));
}
#endif

static VALUE rb_gsl_fcmp(int argc, VALUE *argv, VALUE obj)
{
  double a, b, epsilon = 1e-10;
  switch (argc) {
  case 3:
    epsilon = NUM2DBL(argv[2]);
    /* no break, do next */
  case 2:
    a = NUM2DBL(argv[0]);
    b = NUM2DBL(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  return INT2FIX(gsl_fcmp(a, b, epsilon));
}

static VALUE rb_gsl_equal(int argc, VALUE *argv, VALUE obj)
{
  double a, b, epsilon = 1e-10;
  int retval;

  switch (argc) {
  case 3:
    epsilon = NUM2DBL(argv[2]);
    /* no break, do next */
  case 2:
    a = NUM2DBL(argv[0]);
    b = NUM2DBL(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  retval = gsl_fcmp(a, b, epsilon);
  if (retval == 0) return Qtrue;
  else return Qfalse;
}

void Init_gsl_math(VALUE module)
{
  rb_gsl_define_const(module);

  rb_define_module_function(module, "posinf", rb_GSL_POSINF, 0);
  rb_define_module_function(module, "neginf", rb_GSL_NEGINF, 0);
  rb_define_module_function(module, "nan", rb_GSL_NAN, 0);
  rb_define_module_function(module, "isnan", rb_gsl_isnan, 1);
  rb_define_module_function(module, "isnan?", rb_gsl_isnan2, 1);
  rb_define_module_function(module, "isinf", rb_gsl_isinf, 1);
  rb_define_module_function(module, "isinf?", rb_gsl_isinf2, 1);
  rb_define_module_function(module, "finite", rb_gsl_finite, 1);
  rb_define_module_function(module, "finite?", rb_gsl_finite2, 1);
			    
  rb_define_module_function(module, "sqrt", rb_gsl_sqrt, 1);
  rb_define_module_function(module, "log1p", rb_gsl_log1p, 1);
  rb_define_module_function(module, "expm1", rb_gsl_expm1, 1);
  rb_define_module_function(module, "hypot", rb_gsl_hypot, 2);
#ifdef GSL_1_10_LATER
  rb_define_module_function(module, "hypot3", rb_gsl_hypot3, 3);  
#endif
  rb_define_module_function(module, "acosh", rb_gsl_acosh, 1);
  rb_define_module_function(module, "asinh", rb_gsl_asinh, 1);
  rb_define_module_function(module, "atanh", rb_gsl_atanh, 1);
  rb_define_module_function(module, "pow", rb_gsl_pow, 2);
  rb_define_module_function(module, "pow_int", rb_gsl_pow_int, 2);
  rb_define_module_function(module, "pow_2", rb_gsl_pow_2, 1);
  rb_define_module_function(module, "pow_3", rb_gsl_pow_3, 1);
  rb_define_module_function(module, "pow_4", rb_gsl_pow_4, 1);
  rb_define_module_function(module, "pow_5", rb_gsl_pow_5, 1);
  rb_define_module_function(module, "pow_6", rb_gsl_pow_6, 1);
  rb_define_module_function(module, "pow_7", rb_gsl_pow_7, 1);
  rb_define_module_function(module, "pow_8", rb_gsl_pow_8, 1);
  rb_define_module_function(module, "pow_9", rb_gsl_pow_9, 1);
  rb_define_module_function(module, "sign", rb_GSL_SIGN, 1);
  rb_define_module_function(module, "SIGN", rb_GSL_SIGN, 1);
  
  rb_define_module_function(module, "is_odd", rb_GSL_IS_ODD, 1);
  rb_define_module_function(module, "IS_ODD", rb_GSL_IS_ODD, 1);
  rb_define_module_function(module, "is_odd?", rb_GSL_IS_ODD2, 1);
  rb_define_module_function(module, "IS_ODD?", rb_GSL_IS_ODD2, 1);
  rb_define_module_function(module, "is_even", rb_GSL_IS_EVEN, 1);
  rb_define_module_function(module, "IS_EVEN", rb_GSL_IS_EVEN, 1);
  rb_define_module_function(module, "is_even?", rb_GSL_IS_EVEN2, 1);
  rb_define_module_function(module, "IS_EVEN?", rb_GSL_IS_EVEN2, 1);
  rb_define_module_function(module, "max", rb_GSL_MAX, 2);
  rb_define_module_function(module, "MAX", rb_GSL_MAX, 2);
  rb_define_module_function(module, "min", rb_GSL_MIN, 2);
  rb_define_module_function(module, "MIN", rb_GSL_MIN, 2);

  rb_define_module_function(module, "max_dbl", rb_GSL_MAX_DBL, 2);
  rb_define_module_function(module, "MAX_DBL", rb_GSL_MAX_DBL, 2);
  rb_define_module_function(module, "min_dbl", rb_GSL_MIN_DBL, 2);
  rb_define_module_function(module, "MIN_DBL", rb_GSL_MIN_DBL, 2);

  rb_define_module_function(module, "max_int", rb_GSL_MAX_INT, 2);
  rb_define_module_function(module, "MAX_INT", rb_GSL_MAX_INT, 2);
  rb_define_module_function(module, "min_int", rb_GSL_MIN_INT, 2);
  rb_define_module_function(module, "MIN_INT", rb_GSL_MIN_INT, 2);

  rb_define_module_function(module, "fcmp", rb_gsl_fcmp, -1);
  rb_define_singleton_method(module, "equal?", rb_gsl_equal, -1);
#ifdef GSL_1_3_LATER
  rb_define_module_function(module, "ldexp", rb_gsl_ldexp, 2);
  rb_define_module_function(module, "frexp", rb_gsl_frexp, 1);
#endif
}
