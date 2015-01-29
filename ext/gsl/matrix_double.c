/*
  matrix_double.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada
  
  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_array.h"
#include "include/rb_gsl_complex.h"
#ifdef HAVE_NARRAY_H
#include "include/rb_gsl_with_narray.h"
#endif

enum {
  GSL_MATRIX_ADD,
  GSL_MATRIX_SUB,
  GSL_MATRIX_MUL,
  GSL_MATRIX_DIV,
};

static VALUE rb_gsl_matrix_arithmetics(int flag, VALUE obj, VALUE bb);
VALUE rb_gsl_linalg_LU_solve(int argc, VALUE *argv, VALUE obj);

static VALUE rb_gsl_matrix_arithmetics(int flag, VALUE obj, VALUE bb)
{
  gsl_matrix *m = NULL, *mb = NULL, *mnew = NULL;
  gsl_matrix_complex *cm = NULL, *cmb = NULL, *cmnew = NULL;
  gsl_complex *c = NULL;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_vector_complex *cv = NULL, *cvnew = NULL;
  Data_Get_Struct(obj, gsl_matrix, m);
  switch (TYPE(bb)) {
  case T_FLOAT:
  case T_FIXNUM:
    switch (flag) {
    case GSL_MATRIX_ADD:
      mnew = make_matrix_clone(m);
      gsl_matrix_add_constant(mnew, NUM2DBL(bb));
      break;
    case GSL_MATRIX_SUB:
      mnew = make_matrix_clone(m);
      gsl_matrix_add_constant(mnew, -NUM2DBL(bb));
      break;
    case GSL_MATRIX_MUL:
      mnew = make_matrix_clone(m);
      gsl_matrix_scale(mnew, NUM2DBL(bb));
      break;
    case GSL_MATRIX_DIV:
      mnew = make_matrix_clone(m);
      gsl_matrix_scale(mnew, 1.0/NUM2DBL(bb));
      break;
    default:
      rb_raise(rb_eRuntimeError, "operation not defined");
      break;
    }
    return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    break;
  default:
    if (MATRIX_INT_P(bb)) bb = rb_gsl_matrix_int_to_f(bb);
    if (VECTOR_INT_P(bb)) bb = rb_gsl_vector_int_to_f(bb);

    if (rb_obj_is_kind_of(bb, cgsl_matrix)) {
      Data_Get_Struct(bb, gsl_matrix, mb);
      switch (flag) {
      case GSL_MATRIX_ADD:
  mnew = make_matrix_clone(m);
  gsl_matrix_add(mnew, mb);
  break;
      case GSL_MATRIX_SUB:
  mnew = make_matrix_clone(m);
  gsl_matrix_sub(mnew,mb);
  break;
      case GSL_MATRIX_MUL:
  mnew = make_matrix_clone(m);
  gsl_matrix_mul_elements(mnew, mb);
  break;
      case GSL_MATRIX_DIV:
  mnew = make_matrix_clone(m);
  gsl_matrix_div_elements(mnew, mb);
  break;
      default:
  rb_raise(rb_eRuntimeError, "operation not defined");
  break;
      }
      return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
    } else if (rb_obj_is_kind_of(bb, cgsl_matrix_complex)) {
      Data_Get_Struct(bb, gsl_matrix_complex, cmb);
      cmnew = matrix_to_complex(m);
      switch (flag) {
      case GSL_MATRIX_ADD:
  gsl_matrix_complex_add(cmnew, cmb);
  break;
      case GSL_MATRIX_SUB:
  gsl_matrix_complex_sub(cmnew,cmb);
  break;
      case GSL_MATRIX_MUL:
  gsl_matrix_complex_mul_elements(cmnew, cmb);
  break;
      case GSL_MATRIX_DIV:
  gsl_matrix_complex_div_elements(cmnew, cmb);
  break;
      default:
  rb_raise(rb_eRuntimeError, "operation not defined");
  break;
      }
      return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
    } else if (rb_obj_is_kind_of(bb, cgsl_complex)) {
      Data_Get_Struct(bb, gsl_complex, c);
      cmnew = matrix_to_complex(m);
      switch (flag) {
      case GSL_MATRIX_ADD:
  gsl_matrix_complex_add_constant(cmnew, *c);
  break;
      case GSL_MATRIX_SUB:
  gsl_matrix_complex_add_constant(cmnew, gsl_complex_negative(*c));
  break;
      case GSL_MATRIX_MUL:
  gsl_matrix_complex_scale(cmnew, *c);
  break;
      case GSL_MATRIX_DIV:
  gsl_matrix_complex_scale(cmnew, gsl_complex_inverse(*c));
  break;
      default:
  rb_raise(rb_eRuntimeError, "operation not defined");
  break;
      }
      return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
    } else if (rb_obj_is_kind_of(bb, cgsl_vector)) {
      if (!VECTOR_COL_P(bb))
  rb_raise(rb_eTypeError, 
     "Operation with %s is not defined (GSL::Vector::Col expected)", 
     rb_class2name(CLASS_OF(bb)));
      Data_Get_Struct(bb, gsl_vector, v);
      switch (flag) {
      case GSL_MATRIX_MUL:
  //  vnew = gsl_vector_alloc(v->size);
  vnew = gsl_vector_alloc(m->size1);
  if (vnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
  gsl_matrix_mul_vector(vnew, m, v);
  return Data_Wrap_Struct(cgsl_vector_col, 0, gsl_vector_free, vnew);
  break;
      case GSL_MATRIX_DIV:
  return rb_gsl_linalg_LU_solve(1, &bb, obj);
  break;
      default:
  rb_raise(rb_eRuntimeError, "operation is not defined %s and Matrix",
     rb_class2name(CLASS_OF(bb)));
  break;
      }
    } else if (rb_obj_is_kind_of(bb, cgsl_vector_complex)) {
      Data_Get_Struct(bb, gsl_vector_complex, cv);
      switch (flag) {
      case GSL_MATRIX_MUL:
  cm = matrix_to_complex(m);
  //  cvnew = gsl_vector_complex_alloc(cv->size);
  cvnew = gsl_vector_complex_alloc(m->size1);
  if (cvnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
  gsl_matrix_complex_mul_vector(cvnew, cm, cv);
  gsl_matrix_complex_free(cm);
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cvnew);
  break;
      default:
  rb_raise(rb_eRuntimeError, "operation is not defined %s and Matrix",
     rb_class2name(CLASS_OF(bb)));
  break;
      }
    } else {
      rb_raise(rb_eTypeError, "wrong argument type %s", rb_class2name(CLASS_OF(bb)));
    }
  }
  /* never reach here */
  return Qnil;
}

static VALUE rb_gsl_matrix_add(VALUE obj, VALUE bb)
{
  return rb_gsl_matrix_arithmetics(GSL_MATRIX_ADD, obj, bb);
}

static VALUE rb_gsl_matrix_sub(VALUE obj, VALUE bb)
{
  return rb_gsl_matrix_arithmetics(GSL_MATRIX_SUB, obj, bb);
}

VALUE rb_gsl_matrix_mul_elements(VALUE obj, VALUE bb)
{
  return rb_gsl_matrix_arithmetics(GSL_MATRIX_MUL, obj, bb);
}

/* matrix multiplication */
static VALUE rb_gsl_matrix_mul(VALUE obj, VALUE bb)
{
  gsl_matrix *m = NULL, *b = NULL, *mnew = NULL;
  gsl_vector *v = NULL, *vnew = NULL;
  gsl_matrix_complex *mc, *mcb, *mcnew;
  gsl_vector_complex *vc, *vcnew;
  gsl_complex za, zb;
  Data_Get_Struct(obj, gsl_matrix, m);
  if (VECTOR_INT_P(bb)) bb = rb_gsl_vector_int_to_f(bb);
  if (MATRIX_P(bb)) {
    Data_Get_Struct(bb, gsl_matrix, b);
    mnew = gsl_matrix_alloc(m->size1, b->size2);
    gsl_linalg_matmult(m, b, mnew);
    return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
  } else if (VECTOR_P(bb)) {
    Data_Get_Struct(bb, gsl_vector, v);
    //    vnew = gsl_vector_alloc(v->size);
    //    printf("%d %d\n", m->size1, m->size2);
    vnew = gsl_vector_alloc(m->size1);
    if (vnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
    gsl_matrix_mul_vector(vnew, m, v);
    return Data_Wrap_Struct(VECTOR_ROW_COL(bb), 0, gsl_vector_free, vnew);
  } else if (MATRIX_COMPLEX_P(bb)) {
    Data_Get_Struct(bb, gsl_matrix_complex, mcb);
    mc = matrix_to_complex(m);
    mcnew = gsl_matrix_complex_alloc(m->size1, m->size2);
    gsl_matrix_complex_mul(mcnew, mc, mcb);
    gsl_matrix_complex_free(mc);
    return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, mcnew);
  } else if (VECTOR_COMPLEX_P(bb)) {
    Data_Get_Struct(bb, gsl_vector_complex, vc);
    vcnew = gsl_vector_complex_calloc(m->size1);
    mc = matrix_to_complex(m);
    za.dat[0] = 1.0; za.dat[1] = 0.0;
    zb.dat[0] = 0.0; zb.dat[1] = 0.0;
    gsl_blas_zgemv(CblasNoTrans, za, mc, vc, zb, vcnew);
    gsl_matrix_complex_free(mc);
    return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, vcnew);
  } else {
    switch (TYPE(bb)) {
    case T_FIXNUM:    case T_BIGNUM:    case T_FLOAT:
      return rb_gsl_matrix_mul_elements(obj, bb);
      break;
    default:
      rb_raise(rb_eTypeError, 
         "wrong argument type %s", rb_class2name(CLASS_OF(bb)));
      break;
    }
  }
  return Qnil;
}

static VALUE rb_gsl_matrix_mul_bang(VALUE obj, VALUE bb)
{
  gsl_matrix *m = NULL, *b = NULL, *mtmp = NULL;
  gsl_vector *v = NULL, *vnew = NULL;
  Data_Get_Struct(obj, gsl_matrix, m);
  if (MATRIX_P(bb)) {
    Data_Get_Struct(bb, gsl_matrix, b);
    mtmp = gsl_matrix_alloc(m->size1, b->size2);
    gsl_linalg_matmult(m, b, mtmp);
    gsl_matrix_memcpy(m, mtmp);
    gsl_matrix_free(mtmp);
    return obj;
  } else if (VECTOR_P(bb)) {
    Data_Get_Struct(bb, gsl_vector, v);
    vnew = gsl_vector_alloc(v->size);
    if (vnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_alloc failed");
    gsl_matrix_mul_vector(vnew, m, v);
    return Data_Wrap_Struct(VECTOR_ROW_COL(bb), 0, gsl_vector_free, vnew);
  } else {
    switch (TYPE(bb)) {
    case T_FIXNUM:
    case T_BIGNUM:
    case T_FLOAT:
      gsl_matrix_scale(m, NUM2DBL(bb));
      return obj;
      break;
    default:
      rb_raise(rb_eTypeError, 
         "wrong argument type %s", rb_class2name(CLASS_OF(bb)));
      break;
    }
  }
}

static VALUE rb_gsl_matrix_div_elements(VALUE obj, VALUE bb)
{
  return rb_gsl_matrix_arithmetics(GSL_MATRIX_DIV, obj, bb);
}

static VALUE rb_gsl_matrix_to_complex(VALUE obj)
{
  gsl_matrix *m = NULL;
  gsl_matrix_complex *cm = NULL;
  Data_Get_Struct(obj, gsl_matrix, m);
  cm = matrix_to_complex(m);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cm);
}

gsl_matrix_view* gsl_matrix_view_alloc()
{
  gsl_matrix_view *mv = NULL;
  mv = ALLOC(gsl_matrix_view);
  if (mv == NULL) rb_raise(rb_eRuntimeError, "malloc failed");
  return mv;
}

void gsl_matrix_view_free(gsl_matrix_view * mv);
void gsl_matrix_view_free(gsl_matrix_view * mv)
{
  free((gsl_matrix_view *) mv);
}

static VALUE rb_gsl_matrix_coerce(VALUE obj, VALUE other)
{
  gsl_matrix *m = NULL, *mnew = NULL;
  gsl_matrix_complex *cm = NULL;
  double x;
  gsl_complex *z = NULL;
  VALUE vcm;
  Data_Get_Struct(obj, gsl_matrix, m);
  switch (TYPE(other)) {
  case T_FLOAT:
  case T_FIXNUM:
    mnew = gsl_matrix_alloc(m->size1, m->size2);
    if (mnew == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_alloc failed");
    x = NUM2DBL(other);
    gsl_matrix_set_all(mnew, x);
    return rb_ary_new3(2, Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew), obj);
    break;
  default:
    if (rb_obj_is_kind_of(other, cgsl_complex)) {
      Data_Get_Struct(other, gsl_complex, z);
      cm = gsl_matrix_complex_alloc(m->size1, m->size2);
      if (cm == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
      gsl_matrix_complex_set_all(cm, *z);
      vcm =  Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cm);
      return rb_ary_new3(2, vcm, rb_gsl_matrix_to_complex(obj));
    } else if (rb_obj_is_kind_of(other, cgsl_matrix_complex)) {
      cm = matrix_to_complex(m);
      vcm = Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cm);
      return rb_ary_new3(2, other, vcm);
    } else {
      rb_raise(rb_eTypeError, "cannot coerce %s to Matrix", 
         rb_class2name(CLASS_OF(other)));
    }
    break;
  }
}

static VALUE rb_gsl_matrix_to_f(VALUE obj)
{
  return obj;
}

VALUE rb_gsl_matrix_to_i(VALUE obj)
{
  gsl_matrix *m = NULL;
  gsl_matrix_int *mi = NULL;
  size_t i, j;
  Data_Get_Struct(obj, gsl_matrix, m);
  mi = gsl_matrix_int_alloc(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      gsl_matrix_int_set(mi, i, j, (int) gsl_matrix_get(m, i, j));
    }
  }
  return Data_Wrap_Struct(cgsl_matrix_int, 0, gsl_matrix_int_free, mi);
}

static VALUE rb_gsl_matrix_clean(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m = NULL, *mnew = NULL;
  double eps = 1e-10;
  size_t n, i;
  switch (argc) {
  case 0:
    /* do nothing */
    break;
  case 1:
    Need_Float(argv[0]);
    eps = NUM2DBL(argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments (%d for 0 or 1)", argc);
    break;
  }
  Data_Get_Struct(obj, gsl_matrix, m);
  mnew = make_matrix_clone(m);
  n = m->size1*m->size2;
  for (i = 0; i < n; i++) if (fabs(mnew->data[i]) < eps) mnew->data[i] = 0.0;
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);;
}

static VALUE rb_gsl_matrix_clean_bang(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *m = NULL;
  double eps = 1e-10;
  size_t n, i;
  switch (argc) {
  case 0:
    /* do nothing */
    break;
  case 1:
    Need_Float(argv[0]);
    eps = NUM2DBL(argv[0]);
    break;
  default:
    rb_raise(rb_eArgError, "too many arguments (%d for 0 or 1)", argc);
    break;
  }
  Data_Get_Struct(obj, gsl_matrix, m);
  n = m->size1*m->size2;
  for (i = 0; i < n; i++) if (fabs(m->data[i]) < eps) m->data[i] = 0.0;
  return obj;
}


static VALUE rb_gsl_matrix_op_inplace(VALUE mm1, VALUE mm2, 
              int (*f)(gsl_matrix*, const gsl_matrix*))
{
  gsl_matrix *m1, *m2;
  Data_Get_Struct(mm1, gsl_matrix, m1);
  Data_Get_Struct(mm2, gsl_matrix, m2);
  (*f)(m1, m2);
  return mm1;
}

static VALUE rb_gsl_matrix_add_inplace(VALUE mm1, VALUE mm2)
{
  return rb_gsl_matrix_op_inplace(mm1, mm2, gsl_matrix_add);
}

static VALUE rb_gsl_matrix_sub_inplace(VALUE mm1, VALUE mm2)
{
  return rb_gsl_matrix_op_inplace(mm1, mm2, gsl_matrix_sub);
}

VALUE rb_gsl_sf_eval1(double (*func)(double), VALUE argv);

static VALUE rb_gsl_matrix_sin(VALUE obj)
{
  return rb_gsl_sf_eval1(sin, obj);
}

static VALUE rb_gsl_matrix_cos(VALUE obj)
{
  return rb_gsl_sf_eval1(cos, obj);
}

static VALUE rb_gsl_matrix_tan(VALUE obj)
{
  return rb_gsl_sf_eval1(tan, obj);
}

static VALUE rb_gsl_matrix_exp(VALUE obj)
{
  return rb_gsl_sf_eval1(exp, obj);
}

static VALUE rb_gsl_matrix_log(VALUE obj)
{
  return rb_gsl_sf_eval1(log, obj);
}

static VALUE rb_gsl_matrix_log10(VALUE obj)
{
  return rb_gsl_sf_eval1(log10, obj);
}

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "include/rb_gsl_rng.h"
static VALUE rb_gsl_matrix_randx(int argc, VALUE *argv, VALUE klass,
         double (*f)(const gsl_rng*));
static VALUE rb_gsl_matrix_rand(int argc, VALUE *argv, VALUE klass)
{
  return rb_gsl_matrix_randx(argc, argv, klass, gsl_rng_uniform);
}

static VALUE rb_gsl_matrix_randn(int argc, VALUE *argv, VALUE klass)
{
  return rb_gsl_matrix_randx(argc, argv, klass, gsl_ran_ugaussian);
}

static VALUE rb_gsl_matrix_randx(int argc, VALUE *argv, VALUE klass,
         double (*f)(const gsl_rng*))
{
  gsl_matrix *m;
  gsl_rng *rng;
  size_t size1, size2;
  size_t i, j;
  switch (argc) {
  case 3:
    if (!rb_obj_is_kind_of(argv[2], cgsl_rng)) {
      rb_raise(rb_eTypeError, 
         "Wrong argument type (GSL::Rng expected)");
    }
    Data_Get_Struct(argv[2], gsl_rng, rng);
    size1 = FIX2INT(argv[0]);
    size2 = FIX2INT(argv[1]);
    break;
  case 2:
    size1 = FIX2INT(argv[0]);
    size2 = FIX2INT(argv[1]);
    rng = gsl_rng_alloc(gsl_rng_default);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 2 or 3)", argc);
  }
  m = gsl_matrix_alloc(size1, size2);
  for (i = 0; i < size1; i++) {
    for (j = 0; j < size2; j++) {
      gsl_matrix_set(m, i, j, (*f)(rng));
    }
  }
  if (argc == 2) gsl_rng_free(rng);
  return Data_Wrap_Struct(klass, 0, gsl_matrix_free, m);
}

void Init_gsl_matrix_init(VALUE module);
void Init_gsl_matrix(VALUE module)
{
  Init_gsl_matrix_init(module);

  rb_define_method(cgsl_matrix, "add", rb_gsl_matrix_add, 1);
  rb_define_alias(cgsl_matrix, "+", "add");
  rb_define_method(cgsl_matrix, "sub", rb_gsl_matrix_sub, 1);
  rb_define_alias(cgsl_matrix, "-", "sub");
  rb_define_method(cgsl_matrix, "mul_elements", rb_gsl_matrix_mul_elements, 1);
  rb_define_method(cgsl_matrix, "div_elements", rb_gsl_matrix_div_elements, 1);
  rb_define_alias(cgsl_matrix, "/", "div_elements");

  rb_define_method(cgsl_matrix, "mul", rb_gsl_matrix_mul, 1);
  rb_define_alias(cgsl_matrix, "*", "mul");
  rb_define_method(cgsl_matrix, "mul!", rb_gsl_matrix_mul_bang, 1);
  /***/
 
  rb_define_method(cgsl_matrix, "to_complex", rb_gsl_matrix_to_complex, 0);

  rb_define_method(cgsl_matrix, "coerce", rb_gsl_matrix_coerce, 1);

  rb_undef_method(cgsl_matrix_view_ro, "set");

  rb_define_method(cgsl_matrix, "to_i", rb_gsl_matrix_to_i, 0);
  rb_define_method(cgsl_matrix, "to_f", rb_gsl_matrix_to_f, 0);

  rb_define_method(cgsl_matrix, "clean", rb_gsl_matrix_clean, -1);
  rb_define_method(cgsl_matrix, "clean!", rb_gsl_matrix_clean_bang, -1);

  /*****/
  rb_define_method(cgsl_matrix, "add!", rb_gsl_matrix_add_inplace, 1);
  rb_define_method(cgsl_matrix, "sub!", rb_gsl_matrix_sub_inplace, 1);

  rb_define_method(cgsl_matrix, "sin", rb_gsl_matrix_sin, 0);
  rb_define_method(cgsl_matrix, "cos", rb_gsl_matrix_cos, 0);
  rb_define_method(cgsl_matrix, "tan", rb_gsl_matrix_tan, 0);
  rb_define_method(cgsl_matrix, "exp", rb_gsl_matrix_exp, 0);
  rb_define_method(cgsl_matrix, "log", rb_gsl_matrix_log, 0);
  rb_define_method(cgsl_matrix, "log10", rb_gsl_matrix_log10, 0);

  rb_define_singleton_method(cgsl_matrix, "rand", rb_gsl_matrix_rand, -1);
  rb_define_singleton_method(cgsl_matrix, "randn", rb_gsl_matrix_randn, -1);

}
