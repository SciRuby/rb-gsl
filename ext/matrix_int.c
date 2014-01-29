/*
  matrix_int.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2006 by Yoshiki Tsunesada
  
  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "rb_gsl_array.h"
#include "rb_gsl_complex.h"
#ifdef HAVE_NARRAY_H
#include "rb_gsl_with_narray.h"
#endif

int gsl_linalg_matmult_int(const gsl_matrix_int *A, 
			   const gsl_matrix_int *B, gsl_matrix_int *C);


VALUE rb_gsl_matrix_to_i(VALUE obj);

static VALUE rb_gsl_matrix_int_to_i(VALUE obj)
{
  return obj;
}

VALUE rb_gsl_matrix_int_to_f(VALUE obj)
{
  gsl_matrix_int *m;
  gsl_matrix *mnew;
  size_t i, j;
  Data_Get_Struct(obj, gsl_matrix_int, m);
  mnew = gsl_matrix_alloc(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      gsl_matrix_set(mnew, i, j, (double) gsl_matrix_int_get(m, i, j));
    }
  }
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, mnew);
}

static VALUE rb_gsl_matrix_int_to_complex(VALUE obj)
{
  gsl_matrix_int *m;
  gsl_matrix_complex *mnew;
  gsl_complex z;
  size_t i, j;
  Data_Get_Struct(obj, gsl_matrix_int, m);
  mnew = gsl_matrix_complex_alloc(m->size1, m->size2);
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++) {
      GSL_SET_REAL(&z, (double) gsl_matrix_int_get(m, i, j));
      GSL_SET_IMAG(&z, 0.0);
      gsl_matrix_complex_set(mnew, i, j, z);
    }
  }
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, mnew);
}

static VALUE rb_gsl_matrix_int_coerce(VALUE obj, VALUE other)
{
  return rb_ary_new3(2, other, rb_gsl_matrix_int_to_f(obj));
}

enum {
  GSL_MATRIX_INT_ADD,
  GSL_MATRIX_INT_SUB,
  GSL_MATRIX_INT_MUL,
  GSL_MATRIX_INT_DIV,
};

static VALUE rb_gsl_matrix_int_operation1(VALUE obj, VALUE other, int flag)
{
  gsl_matrix_int *a, *anew, *b;
  gsl_vector_int *vi, *vinew;
  double bval;
  // local variable "result" declared and set, but never used
  //int result;
  Data_Get_Struct(obj, gsl_matrix_int, a);
  switch (TYPE(other)) {
  case T_FIXNUM:
  case T_FLOAT:
    bval = NUM2INT(other);
    anew = make_matrix_int_clone(a);
    switch (flag) {
    case GSL_MATRIX_INT_ADD: 
      /*result =*/ gsl_matrix_int_add_constant(anew, bval);
      break;
    case GSL_MATRIX_INT_SUB: 
      /*result =*/ gsl_matrix_int_add_constant(anew, -bval);
      break;
    case GSL_MATRIX_INT_MUL: 
      /*result =*/ gsl_matrix_int_scale(anew, bval);
      break;
    case GSL_MATRIX_INT_DIV: 
      /*result =*/ gsl_matrix_int_scale(anew, 1.0/bval);
      break;
    default:
      break;
    }
    break;
  default:
    if (MATRIX_P(other)) other = rb_gsl_matrix_to_i(other);
    if (VECTOR_P(other)) other = rb_gsl_vector_to_i(other);
    if (MATRIX_INT_P(other)) {
      anew = make_matrix_int_clone(a);
      Data_Get_Struct(other, gsl_matrix_int, b);
      switch (flag) {
      case GSL_MATRIX_INT_ADD:
	/*result =*/ gsl_matrix_int_add(anew, b);
	break;
      case GSL_MATRIX_INT_SUB:
	/*result =*/ gsl_matrix_int_sub(anew, b);
	break;
      case GSL_MATRIX_INT_MUL:
	/*result =*/ gsl_matrix_int_mul_elements(anew, b);
	break;
      case GSL_MATRIX_INT_DIV:
	/*result =*/ gsl_matrix_int_div_elements(anew, b);
	break;
      default:
	break;
      }
    } else if (VECTOR_INT_COL_P(other)) {
      switch (flag) {
      case GSL_MATRIX_INT_MUL:
	Data_Get_Struct(other, gsl_vector_int, vi);
	vinew = gsl_vector_int_alloc(vi->size);
	gsl_matrix_int_mul_vector(vinew, a, vi);
	return Data_Wrap_Struct(cgsl_vector_int_col, 0, gsl_vector_int_free, vinew);
	break;
      default:
	rb_raise(rb_eRuntimeError, "Operation not defined");
      }
    } else {
      rb_raise(rb_eTypeError, "Operation not defined with %s",
	       rb_class2name(CLASS_OF(other)));
    }
    break;
  }
  return Data_Wrap_Struct(cgsl_matrix_int, 0, gsl_matrix_int_free, anew);
}

static VALUE rb_gsl_matrix_int_add(VALUE obj, VALUE other)
{
  return rb_gsl_matrix_int_operation1(obj, other, GSL_MATRIX_INT_ADD);
}

static VALUE rb_gsl_matrix_int_sub(VALUE obj, VALUE other)
{
  return rb_gsl_matrix_int_operation1(obj, other, GSL_MATRIX_INT_SUB);
}

static VALUE rb_gsl_matrix_int_mul(VALUE obj, VALUE other)
{
  return rb_gsl_matrix_int_operation1(obj, other, GSL_MATRIX_INT_MUL);
}

static VALUE rb_gsl_matrix_int_div(VALUE obj, VALUE other)
{
  return rb_gsl_matrix_int_operation1(obj, other, GSL_MATRIX_INT_DIV);
}

static VALUE rb_gsl_matrix_int_matrix_mul(VALUE obj, VALUE bb)
{
  gsl_matrix_int *m = NULL, *b = NULL, *mnew = NULL;
  gsl_vector_int *vi, *vinew;
  Data_Get_Struct(obj, gsl_matrix_int, m);
  if (MATRIX_INT_P(bb)) {
    Data_Get_Struct(bb, gsl_matrix_int, b);
    mnew = gsl_matrix_int_alloc(m->size1, b->size2);
    gsl_linalg_matmult_int(m, b, mnew);
    return Data_Wrap_Struct(cgsl_matrix_int, 0, gsl_matrix_int_free, mnew);
  } else {
    if (VECTOR_INT_COL_P(bb)) {
      Data_Get_Struct(bb, gsl_vector_int, vi);
      vinew = gsl_vector_int_alloc(vi->size);
      gsl_matrix_int_mul_vector(vinew, m, vi);
      return Data_Wrap_Struct(cgsl_vector_int_col, 0, gsl_vector_int_free, vinew);
    }
    switch (TYPE(bb)) {
    case T_FIXNUM:
      return rb_gsl_matrix_int_mul(obj, bb);
	/*      return rb_gsl_matrix_int_power(obj, bb);*/
      break;
    default:
      rb_raise(rb_eTypeError, "wrong argument type %s (Matrix::Int, Vector::Int::Col or Fixnum expected)",
	       rb_class2name(CLASS_OF(bb)));
      break;
    }
  }
}

int gsl_linalg_matmult_int(const gsl_matrix_int *A, 
			   const gsl_matrix_int *B, gsl_matrix_int *C)
{
  if (A->size2 != B->size1 || A->size1 != C->size1 || B->size2 != C->size2)
    {
      GSL_ERROR ("matrix sizes are not conformant", GSL_EBADLEN);
    }
  else
    {
      int a, b;
      int temp;
      size_t i, j, k;
      
      for (i = 0; i < C->size1; i++)
        {
          for (j = 0; j < C->size2; j++)
            {
              a = gsl_matrix_int_get(A, i, 0);
              b = gsl_matrix_int_get(B, 0, j);
              temp = a * b;
              for (k = 1; k < A->size2; k++)
                {
                  a = gsl_matrix_int_get(A, i, k);
                  b = gsl_matrix_int_get(B, k, j);
                  temp += a * b;
                }
              gsl_matrix_int_set(C, i, j, temp);
            }
        }

      return GSL_SUCCESS;
    }
}

void Init_gsl_matrix_int_init(VALUE module);
void Init_gsl_matrix_int(VALUE module)
{
  Init_gsl_matrix_int_init(module);
  /*****/

  rb_define_method(cgsl_matrix_int, "to_f", rb_gsl_matrix_int_to_f, 0);
  rb_define_method(cgsl_matrix_int, "to_i", rb_gsl_matrix_int_to_i, 0);
  rb_define_method(cgsl_matrix_int, "to_complex", rb_gsl_matrix_int_to_complex, 0);

  /*****/
  rb_define_method(cgsl_matrix_int, "coerce", rb_gsl_matrix_int_coerce, 1);
  /*****/
  rb_define_method(cgsl_matrix_int, "add", rb_gsl_matrix_int_add, 1);
  rb_define_alias(cgsl_matrix_int, "+", "add");
  rb_define_method(cgsl_matrix_int, "sub", rb_gsl_matrix_int_sub, 1);
  rb_define_alias(cgsl_matrix_int, "-", "sub");
  rb_define_method(cgsl_matrix_int, "mul", rb_gsl_matrix_int_mul, 1);
  /*  rb_define_alias(cgsl_matrix_int, "*", "mul");*/
  rb_define_method(cgsl_matrix_int, "div", rb_gsl_matrix_int_div, 1);
  rb_define_alias(cgsl_matrix_int, "/", "div");

  rb_define_method(cgsl_matrix_int, "matrix_mul", rb_gsl_matrix_int_matrix_mul, 1);
  rb_define_alias(cgsl_matrix_int, "*", "matrix_mul");
  /*****/

}
