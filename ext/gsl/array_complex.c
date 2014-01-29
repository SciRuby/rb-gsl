/*
  array_complex.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_complex.h"
#include "include/rb_gsl_array.h"

enum {
  GSL_COMPLEX_ADD,
  GSL_COMPLEX_SUB,
  GSL_COMPLEX_MUL,
  GSL_COMPLEX_DIV,
};

static VALUE rb_gsl_complex_arithmetics5(int flag, VALUE obj, VALUE bb);

static VALUE rb_gsl_complex_arithmetics5(int flag, VALUE obj, VALUE bb)
{
  gsl_complex *a = NULL, *b = NULL, *c = NULL, tmp, tmp2;
  gsl_matrix *m = NULL;
  gsl_matrix_complex *cm = NULL, *cmself = NULL;
  gsl_vector *v = NULL;
  gsl_vector_complex *cv = NULL, *cvnew = NULL;
  gsl_complex (*func1)(gsl_complex, gsl_complex);
  // local variables "func2" iand "func3" declared and set, but never used
  //int (*func2)(gsl_matrix_complex*, const gsl_matrix_complex*);
  //int (*func3)(gsl_matrix_complex*, const gsl_complex);
  int flagcm = 0;
  switch (flag) {
  case GSL_COMPLEX_ADD:
    func1 = gsl_complex_add;
    //func2 = gsl_matrix_complex_add;
    //func3 = gsl_matrix_complex_add_constant;
    break;
  case GSL_COMPLEX_SUB:
    func1 = gsl_complex_sub;
    //func2 = gsl_matrix_complex_sub;
    //func3 = gsl_matrix_complex_add_constant;
    break;
  case GSL_COMPLEX_MUL:
    func1 = gsl_complex_mul;
    //func2 = gsl_matrix_complex_mul_elements;
    //func3 = gsl_matrix_complex_scale;
    break;
  case GSL_COMPLEX_DIV:
    func1 = gsl_complex_div;
    //func2 = gsl_matrix_complex_div_elements;
    //func3 = gsl_matrix_complex_scale;
    break;
  default:
    rb_raise(rb_eRuntimeError, "undefined operation");
  }

  CHECK_COMPLEX(obj);
  Data_Get_Struct(obj, gsl_complex, a);
  switch (TYPE(bb)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    tmp2 = gsl_complex_rect(NUM2DBL(bb), 0.0);
    b = &tmp2;
    tmp = (*func1)(*a, *b);
    switch (flag) {
    case GSL_COMPLEX_ADD:
    case GSL_COMPLEX_SUB:
    case GSL_COMPLEX_MUL:
    case GSL_COMPLEX_DIV:
      c = ALLOC(gsl_complex);
      *c = tmp;
      return Data_Wrap_Struct(cgsl_complex, 0, free, c); 
      break;
    }
    break;
  default:
    if (COMPLEX_P(bb)) {
      Data_Get_Struct(bb, gsl_complex, b);
      tmp = (*func1)(*a, *b);
      switch (flag) {
      case GSL_COMPLEX_ADD:
      case GSL_COMPLEX_SUB:
      case GSL_COMPLEX_MUL:
      case GSL_COMPLEX_DIV:
	c = ALLOC(gsl_complex);
	*c = tmp;
	return Data_Wrap_Struct(cgsl_complex, 0, free, c); 
	break;
      }
    } else {
      if (VECTOR_P(bb)) {
	Data_Get_Struct(bb, gsl_vector, v);
	cv = vector_to_complex(v);
	cvnew = gsl_vector_complex_alloc(v->size);
	if (cvnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
	gsl_vector_complex_set_all(cvnew, *a);
	switch (flag) {
	case GSL_COMPLEX_ADD:
	  gsl_vector_complex_add(cvnew, cv);
	  break;
	case GSL_COMPLEX_SUB:
	  gsl_vector_complex_sub(cvnew, cv);
	  break;
	case GSL_COMPLEX_MUL:
	  gsl_vector_complex_mul(cvnew, cv);
	  break;
	case GSL_COMPLEX_DIV:
	  gsl_vector_complex_add(cvnew, cv);
	  break;
	}
	gsl_vector_complex_free(cv);
	return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cvnew);
      }
      if (VECTOR_COMPLEX_P(bb)) {
	Data_Get_Struct(bb, gsl_vector_complex, cv);
	cvnew = gsl_vector_complex_alloc(v->size);
	if (cvnew == NULL) rb_raise(rb_eNoMemError, "gsl_vector_complex_alloc failed");
	gsl_vector_complex_set_all(cvnew, *a);
	switch (flag) {
	case GSL_COMPLEX_ADD:
	  gsl_vector_complex_add(cvnew, cv);
	  break;
	case GSL_COMPLEX_SUB:
	  gsl_vector_complex_sub(cvnew, cv);
	  break;
	case GSL_COMPLEX_MUL:
	  gsl_vector_complex_mul(cvnew, cv);
	  break;
	case GSL_COMPLEX_DIV:
	  gsl_vector_complex_add(cvnew, cv);
	  break;
	}
	return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free, cvnew);
      }


      if (MATRIX_P(bb)) {
	Data_Get_Struct(bb, gsl_matrix, m);
	cm = matrix_to_complex(m);
	flagcm = 1;
      }	else if (MATRIX_COMPLEX_P(bb)) {
	Data_Get_Struct(bb, gsl_matrix_complex, cm);
      } else {
	rb_raise(rb_eTypeError, "wrong argument type %s", rb_class2name(CLASS_OF(bb)));
      }
      cmself = gsl_matrix_complex_alloc(m->size1, m->size2);
      if (cmself == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
      gsl_matrix_complex_set_all(cmself, *a);
      switch (flag) {
      case GSL_COMPLEX_ADD:
	gsl_matrix_complex_add(cmself, cm);
	break;
      case GSL_COMPLEX_SUB:
	gsl_matrix_complex_sub(cmself, cm);
	break;
      case GSL_COMPLEX_MUL:
	gsl_matrix_complex_mul_elements(cmself, cm);
	break;
      case GSL_COMPLEX_DIV:
	gsl_matrix_complex_div_elements(cmself, cm);
	break;
      }
      if (flagcm == 1) gsl_matrix_complex_free(cm);
      return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmself);
    }
  }
  /* never reach here */
  return Qnil;
}

static VALUE rb_gsl_complex_add(VALUE obj, VALUE bb)
{
  return rb_gsl_complex_arithmetics5(GSL_COMPLEX_ADD, obj, bb);
}

static VALUE rb_gsl_complex_sub(VALUE obj, VALUE bb)
{
  return rb_gsl_complex_arithmetics5(GSL_COMPLEX_SUB, obj, bb);
}

static VALUE rb_gsl_complex_mul(VALUE obj, VALUE bb)
{
  return rb_gsl_complex_arithmetics5(GSL_COMPLEX_MUL, obj, bb);
}

static VALUE rb_gsl_complex_div(VALUE obj, VALUE bb)
{
  return rb_gsl_complex_arithmetics5(GSL_COMPLEX_DIV, obj, bb);
}

static VALUE rb_gsl_complex_coerce(VALUE obj, VALUE other)
{
  gsl_complex *c = NULL;
  gsl_matrix *m = NULL;
  gsl_matrix_complex *cmnew = NULL, *cmself = NULL;
  VALUE vcmself, vcmnew;
  double x;
  switch (TYPE(other)) {
  case T_FLOAT:  case T_FIXNUM:  case T_BIGNUM:
    x = NUM2DBL(other);
    c = ALLOC(gsl_complex);
    *c = gsl_complex_rect(x, 0.0);
    return rb_ary_new3(2, Data_Wrap_Struct(cgsl_complex, 0, free, c),
		       obj);
    break;
  default:
    if (MATRIX_P(other)) {
      Data_Get_Struct(other, gsl_matrix, m);
      cmnew = matrix_to_complex(m);
      vcmnew = Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmnew);
      cmself = gsl_matrix_complex_alloc(m->size1, m->size2);
      if (cmself == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
      Data_Get_Struct(obj, gsl_complex, c);
      gsl_matrix_complex_set_all(cmself, *c);
      vcmself = Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmself);
      return rb_ary_new3(2, vcmself, vcmnew);
    } 
    if (MATRIX_COMPLEX_P(other)) {
      Data_Get_Struct(other, gsl_matrix_complex, cmnew);
      cmself = gsl_matrix_complex_alloc(cmnew->size1, cmnew->size2);
      if (cmself == NULL) rb_raise(rb_eNoMemError, "gsl_matrix_complex_alloc failed");
      vcmself = Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, cmself);
      return rb_ary_new3(2, vcmself, other);
    } else {
      rb_raise(rb_eTypeError, "cannot coerce to GSL::Complex");
    }
  }
}

void Init_gsl_array_complex(VALUE mgsl)
{
  rb_define_method(cgsl_complex, "coerce", rb_gsl_complex_coerce, 1);

  rb_define_method(cgsl_complex, "add", rb_gsl_complex_add, 1);
  rb_define_alias(cgsl_complex, "+", "add");
  rb_define_method(cgsl_complex, "sub", rb_gsl_complex_sub, 1);
  rb_define_alias(cgsl_complex, "-", "sub");
  rb_define_method(cgsl_complex, "mul", rb_gsl_complex_mul, 1);
  rb_define_alias(cgsl_complex, "*", "mul");
  rb_define_method(cgsl_complex, "div", rb_gsl_complex_div, 1);
  rb_define_alias(cgsl_complex, "/", "div");
}
