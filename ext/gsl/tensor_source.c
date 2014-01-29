/*
  tensor_source.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

/*
  The tensor package is developed by Jordi Burguet-Caltell, and
  distributed separately as an add-on package.
  http://savannah.nongnu.org/projects/tensor/
 */

#ifdef HAVE_TENSOR_TENSOR_H

#include "rb_gsl_tensor.h"
#include "rb_gsl_common.h"

#ifdef BASE_DOUBLE
VALUE cgsl_tensor, cgsl_tensor_int;
VALUE cgsl_tensor_view, cgsl_tensor_int_view;
#define NUMCONV(x) NUM2DBL(x)
#define C_TO_VALUE(x) rb_float_new(x)
#define CHECK_TEN(x) CHECK_TENSOR(x)
#define TEN_P(x) TENSOR_P(x)
#define VEC_P(x) VECTOR_P(x)
#define MAT_P(x) MATRIX_P(x)
#else defined(BASE_INT)
#define NUMCONV(x) FIX2INT(x)
#define C_TO_VALUE(x) INT2FIX(x)
#define CHECK_TEN(x) CHECK_TENSOR_INT(x)
#define TEN_P(x) TENSOR_INT_P(x)
#define VEC_P(x) VECTOR_INT_P(x)
#define MAT_P(x) MATRIX_INT_P(x)
#endif

GSL_TYPE(rbgsl_tensor)* FUNCTION(rbgsl_tensor,alloc)(const unsigned int rank, 
						     const size_t dimension)
{
  GSL_TYPE(rbgsl_tensor) *t;
  t = ALLOC(GSL_TYPE(rbgsl_tensor));
  t->tensor = FUNCTION(tensor,alloc)(rank, dimension);
  if (rank == 0) 
    t->indices = gsl_permutation_alloc(1);
  else
    t->indices = gsl_permutation_alloc(rank);
  return t;
}

GSL_TYPE(rbgsl_tensor)* FUNCTION(rbgsl_tensor,calloc)(const unsigned int rank, 
						      const size_t dimension)
{
  GSL_TYPE(rbgsl_tensor) *t;
  t = ALLOC(GSL_TYPE(rbgsl_tensor));
  t->tensor = FUNCTION(tensor,calloc)(rank, dimension);
  if (rank == 0) 
    t->indices = gsl_permutation_alloc(1);
  else
    t->indices = gsl_permutation_alloc(rank);
  return t;
}

GSL_TYPE(rbgsl_tensor)* FUNCTION(rbgsl_tensor,copy)(const GSL_TYPE(rbgsl_tensor) *t)
{
  GSL_TYPE(rbgsl_tensor) *tnew;
  tnew = ALLOC(GSL_TYPE(rbgsl_tensor));
  if (t->tensor->rank == 0) 
    tnew->indices = gsl_permutation_alloc(1);
  else
    tnew->indices = gsl_permutation_alloc(t->tensor->rank);
  tnew->tensor = FUNCTION(tensor,copy)(t->tensor);
  return tnew;
}

void FUNCTION(rbgsl_tensor,free)(GSL_TYPE(rbgsl_tensor) *t)
{
  gsl_permutation_free(t->indices);
  FUNCTION(tensor,free)(t->tensor);
  free((GSL_TYPE(rbgsl_tensor) *) t);
}

void FUNCTION(rbgsl_tensor,free2)(GSL_TYPE(rbgsl_tensor) *t)
{
  gsl_permutation_free(t->indices);
  free((GSL_TYPE(tensor)*) t->tensor);
  free((GSL_TYPE(rbgsl_tensor) *) t);
}

/* singleton methods */
static VALUE FUNCTION(rb_tensor,new)(int argc, VALUE *argv, VALUE klass)
{
  unsigned int rank;
  size_t dim;
  GSL_TYPE(rbgsl_tensor) *t = NULL;
  switch (argc) {
  case 2:
    rank = FIX2UINT(argv[0]);
    dim = (size_t) FIX2UINT(argv[1]);
    t = FUNCTION(rbgsl_tensor,alloc)(rank, dim);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2, rank and dimension)",
	     argc);
    break;
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), t);
}

static VALUE FUNCTION(rb_tensor,calloc)(VALUE klass, VALUE r, VALUE s)
{
  unsigned int rank;
  size_t dim;
  GSL_TYPE(rbgsl_tensor) *t = NULL;
  rank = FIX2UINT(r);
  dim = (size_t) FIX2UINT(s);
  t = FUNCTION(rbgsl_tensor,calloc)(rank, dim);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), t);
}

static VALUE FUNCTION(rb_tensor,copy_singleton)(VALUE klass, VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t = NULL, *tnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  tnew = FUNCTION(rbgsl_tensor,copy(t));
  return Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), tnew);
}

static VALUE FUNCTION(rb_tensor,memcpy_singleton)(VALUE klass, VALUE a, VALUE b)
{
  GSL_TYPE(rbgsl_tensor) *dst, *src;
  CHECK_TEN(b);
  Data_Get_Struct(a, GSL_TYPE(rbgsl_tensor), dst);
  Data_Get_Struct(b, GSL_TYPE(rbgsl_tensor), src);
  return INT2FIX(FUNCTION(tensor,memcpy)(dst->tensor, src->tensor));
}

static VALUE FUNCTION(rb_tensor,swap_singleton)(VALUE klass, VALUE a, VALUE b)
{
  GSL_TYPE(rbgsl_tensor) *t1, *t2;
  CHECK_TEN(b);
  Data_Get_Struct(a, GSL_TYPE(rbgsl_tensor), t1);
  Data_Get_Struct(b, GSL_TYPE(rbgsl_tensor), t2);
  return INT2FIX(FUNCTION(tensor,swap)(t1->tensor, t2->tensor));
}

/*****/
static VALUE FUNCTION(rb_tensor,copy)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t = NULL, *tnew = NULL;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  tnew = FUNCTION(rbgsl_tensor,copy)(t);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), tnew);
}

static VALUE FUNCTION(rb_tensor,set_zero)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t = NULL;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  FUNCTION(tensor,set_zero)(t->tensor);
  return obj;
}

static VALUE FUNCTION(rb_tensor,set_all)(VALUE obj, VALUE xx)
{
  GSL_TYPE(rbgsl_tensor) *t = NULL;
  BASE x;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  x = NUMCONV(xx);
  FUNCTION(tensor,set_all)(t->tensor, x);
  return obj;
}

static void rb_tensor_get_indices_array(tensor_indices *v, VALUE ary);
static void rbgsl_tensor_get_indices(int argc, VALUE *argv, tensor_indices *indices, 
				     size_t *n);
#ifdef BASE_DOUBLE
static void rb_tensor_get_indices_array(tensor_indices *v, VALUE ary)
{
  size_t i, nn;
  //  nn = (size_t) GSL_MIN_INT((int) v->size, (int) RARRAY(ary)->len);
  nn = (size_t) GSL_MIN_INT((int) v->size, (int) RARRAY_LEN(ary));
  for (i = 0; i < nn; i++)
    v->data[i] = FIX2UINT(rb_ary_entry(ary, i));
}

static void rbgsl_tensor_get_indices(int argc, VALUE *argv, 
				     tensor_indices *indices, size_t *n)
{
  size_t i;
  for (i = 0; i < indices->size; i++) indices->data[i] = 0;
  switch (argc) {
  case 1:
    switch (TYPE(argv[0])) {
    case T_ARRAY:
      //      *n = (size_t) GSL_MIN_INT((int) indices->size, (int) RARRAY(argv[0])->len);
      *n = (size_t) GSL_MIN_INT((int) indices->size, (int) RARRAY_LEN(argv[0]));
      rb_tensor_get_indices_array(indices, argv[0]);
      break;
    case T_FIXNUM:
      *n = 1;
      indices->data[0] = FIX2INT(argv[0]);
      break;
    default:
      rb_raise(rb_eTypeError, "wrong argument type %s (Array expected)",
	       rb_class2name(CLASS_OF(argv[0])));
      break;
    }
    break;
  default:
    *n = (size_t) GSL_MIN_INT(argc, (int) indices->size);
    for (i = 0; i < *n; i++) {
      CHECK_FIXNUM(argv[i]);
      indices->data[i] = FIX2INT(argv[i]);
    }
    break;
  }
}
#endif

size_t FUNCTION(tensor,position)(const size_t * indices,
				     const GSL_TYPE(tensor) * t);
static VALUE FUNCTION(rb_tensor,position)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t = NULL;
  size_t n, position;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);  
  rbgsl_tensor_get_indices(argc, argv, t->indices, &n);  
  position = (size_t) FUNCTION(tensor,position)(t->indices->data,t->tensor);
  return INT2FIX(position);
}
static VALUE FUNCTION(rb_tensor,subtensor)(int argc, VALUE *argv, VALUE obj);
static VALUE FUNCTION(rb_tensor,get)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t = NULL;
  BASE x;
  size_t n;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);  
  rbgsl_tensor_get_indices(argc, argv, t->indices, &n);  
  if (n < t->tensor->rank) {
    return FUNCTION(rb_tensor,subtensor)(argc, argv, obj);
  } else {
    x = FUNCTION(tensor,get)(t->tensor, t->indices->data);
    return C_TO_VALUE(x);
  }
  return Qnil;
}

static VALUE FUNCTION(rb_tensor,set)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t = NULL;
  size_t n;
  BASE x;
  if (argc < 2) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for >= 2)", argc);
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);  
  rbgsl_tensor_get_indices(argc-1, argv, t->indices, &n);
  x = NUMCONV(argv[argc-1]);
  FUNCTION(tensor,set)(t->tensor, t->indices->data, x);
  return obj;
}

static VALUE FUNCTION(rb_tensor,fread)(VALUE obj, VALUE io)
{
  GSL_TYPE(rbgsl_tensor) *t = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  f = rb_gsl_open_readfile(io, &flag);
  status = FUNCTION(tensor,fread)(f, t->tensor);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_tensor,fwrite)(VALUE obj, VALUE io)
{
  GSL_TYPE(rbgsl_tensor) *t = NULL;
  FILE *f = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  f = rb_gsl_open_writefile(io, &flag);
  status = FUNCTION(tensor,fwrite)(f, t->tensor);
  if (flag == 1) fclose(f);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_tensor,fprintf)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;
  if (argc != 1 && argc != 2) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), h);
  fp = rb_gsl_open_writefile(argv[0], &flag);
  switch (argc) {
  case 2:
    if (TYPE(argv[1]) == T_STRING)
      status = FUNCTION(tensor,fprintf)(fp, h->tensor, STR2CSTR(argv[1]));
    else
      rb_raise(rb_eTypeError, "argv 2 String expected");
    break;
  default:
    status = FUNCTION(tensor,fprintf)(fp, h->tensor, OUT_FORMAT);
    break;
  }
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_tensor,printf)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *h = NULL;
  int status;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), h);
  if (argc == 1) {
    if (TYPE(argv[0]) != T_STRING) 
      rb_raise(rb_eTypeError, "String expected");
    else
      status = FUNCTION(tensor,fprintf)(stdout, h->tensor, STR2CSTR(argv[0]));
  } else {
    status = FUNCTION(tensor,fprintf)(stdout, h->tensor, OUT_FORMAT);
  }
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_tensor,fscanf)(VALUE obj, VALUE io)
{
  GSL_TYPE(rbgsl_tensor) *h = NULL;
  FILE *fp = NULL;
  int status, flag = 0;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), h);
  fp = rb_gsl_open_readfile(io, &flag);
  status = FUNCTION(tensor,fscanf)(fp, h->tensor);
  if (flag == 1) fclose(fp);
  return INT2FIX(status);
}

static VALUE FUNCTION(rb_tensor,swap_indices)(VALUE obj, VALUE ii, VALUE jj)
{
  GSL_TYPE(rbgsl_tensor) *t, *tnew;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  tnew = ALLOC(GSL_TYPE(rbgsl_tensor));
  if (t->tensor->rank == 0) 
    tnew->indices = gsl_permutation_alloc(1);
  else
    tnew->indices = gsl_permutation_alloc(t->tensor->rank);
  tnew->tensor = FUNCTION(tensor,swap_indices)(t->tensor, FIX2INT(ii), FIX2INT(jj));
  return Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), tnew);
}

static VALUE FUNCTION(rb_tensor,max)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  return C_TO_VALUE(FUNCTION(tensor,max)(t->tensor));
}

static VALUE FUNCTION(rb_tensor,min)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  return C_TO_VALUE(FUNCTION(tensor,min)(t->tensor));
}

static VALUE FUNCTION(rb_tensor,minmax)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  BASE min, max;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  FUNCTION(tensor,minmax)(t->tensor, &min, &max);
  return rb_ary_new3(2, C_TO_VALUE(min), C_TO_VALUE(max));
}

static VALUE FUNCTION(rb_tensor,max_index)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  FUNCTION(tensor,max_index)(t->tensor, t->indices->data);
  return Data_Wrap_Struct(cgsl_index, 0, NULL, t->indices);
}

static VALUE FUNCTION(rb_tensor,min_index)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  FUNCTION(tensor,min_index)(t->tensor, t->indices->data);
  return Data_Wrap_Struct(cgsl_index, 0, NULL, t->indices);
}

static VALUE FUNCTION(rb_tensor,minmax_index)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  gsl_permutation *min, *max;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  if (t->tensor->rank == 0) {
    min = gsl_permutation_alloc(1);
    max = gsl_permutation_alloc(1);
  } else {
    min = gsl_permutation_alloc(t->tensor->rank);
    max = gsl_permutation_alloc(t->tensor->rank);
  }
  FUNCTION(tensor,minmax_index)(t->tensor, min->data, max->data);
  return rb_ary_new3(2,
		     Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, min), 
		     Data_Wrap_Struct(cgsl_index, 0, gsl_permutation_free, max));
}

static VALUE FUNCTION(rb_tensor,isnull)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  return INT2FIX(FUNCTION(tensor,isnull)(t->tensor));
}

static VALUE FUNCTION(rb_tensor,isnull2)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  int status;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  status = FUNCTION(tensor,isnull)(t->tensor);
  if (status) return Qtrue;
  else return Qfalse;
}

static VALUE FUNCTION(rb_tensor,oper)(VALUE obj, VALUE bb,
					   int flag)
{
  GSL_TYPE(rbgsl_tensor) *a, *b, *anew;
  BASE x;
  int (*f)(GSL_TYPE(tensor)*, const GSL_TYPE(tensor)*);
  int (*f2)(GSL_TYPE(tensor)*, const double);
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), a);
  anew = FUNCTION(rbgsl_tensor,copy)(a);
  if (TEN_P(bb)) {
    Data_Get_Struct(bb, GSL_TYPE(rbgsl_tensor), b);
    switch (flag) {
    case TENSOR_ADD: f = FUNCTION(&tensor,add); break;
    case TENSOR_SUB: f = FUNCTION(&tensor,sub); break;
    case TENSOR_MUL_ELEMENTS: f = FUNCTION(&tensor,mul_elements); break;
    case TENSOR_DIV_ELEMENTS: f = FUNCTION(&tensor,div_elements); break;
    default: rb_raise(rb_eRuntimeError, "unknown operation"); break;
    } 
    (*f)(anew->tensor, b->tensor);
  } else {
    switch (flag) {
    case TENSOR_ADD:
    case TENSOR_ADD_CONSTANT: 
      x = NUMCONV(bb);
      f2 = FUNCTION(&tensor,add_constant); 
      break;
    case TENSOR_SUB:
      x = -NUMCONV(bb);
      f2 = FUNCTION(&tensor,add_constant); 
      break;
    case TENSOR_ADD_DIAGONAL: f2 = FUNCTION(&tensor,add_diagonal); break;
    case TENSOR_MUL_ELEMENTS:
    case TENSOR_SCALE:
      x = NUMCONV(bb); 
      f2 = FUNCTION(&tensor,scale); 
      break;
    case TENSOR_DIV_ELEMENTS: 
      x = 1.0/NUMCONV(bb); 
      f2 = FUNCTION(&tensor,scale); 
      break;
    default: rb_raise(rb_eRuntimeError, "unknown operation"); break;
    }
    (*f2)(anew->tensor, x);
  }
  return Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), anew);
}


static VALUE FUNCTION(rb_tensor,oper_bang)(VALUE obj, VALUE bb, int flag)
{
  GSL_TYPE(rbgsl_tensor) *a, *b;
  BASE x;
  int (*f)(GSL_TYPE(tensor)*, const GSL_TYPE(tensor)*);
  int (*f2)(GSL_TYPE(tensor)*, const double);
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), a);
  if (TEN_P(bb)) {
    Data_Get_Struct(bb, GSL_TYPE(rbgsl_tensor), b);
    switch (flag) {
    case TENSOR_ADD: f = FUNCTION(&tensor,add); break;
    case TENSOR_SUB: f = FUNCTION(&tensor,sub); break;
    case TENSOR_MUL_ELEMENTS: f = FUNCTION(&tensor,mul_elements); break;
    case TENSOR_DIV_ELEMENTS: f = FUNCTION(&tensor,div_elements); break;
    default: rb_raise(rb_eRuntimeError, "unknown operation"); break;
    } 
    (*f)(a->tensor, b->tensor);
  } else {
    switch (flag) {
    case TENSOR_ADD:
    case TENSOR_ADD_CONSTANT: 
      x = NUMCONV(bb);
      f2 = FUNCTION(&tensor,add_constant); 
      break;
    case TENSOR_SUB:
      x = -NUMCONV(bb);
      f2 = FUNCTION(&tensor,add_constant); 
      break;
    case TENSOR_ADD_DIAGONAL: f2 = FUNCTION(&tensor,add_diagonal); break;
    case TENSOR_MUL_ELEMENTS:
    case TENSOR_SCALE:
      x = NUMCONV(bb); 
      f2 = FUNCTION(&tensor,scale); 
      break;
    case TENSOR_DIV_ELEMENTS: 
      x = 1.0/NUMCONV(bb); 
      f2 = FUNCTION(&tensor,scale); 
      break;
    default: rb_raise(rb_eRuntimeError, "unknown operation"); break;
    }
    (*f2)(a->tensor, x);
  }
  return obj;
}

static VALUE FUNCTION(rb_tensor,add)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper)(obj, bb, TENSOR_ADD);
}

static VALUE FUNCTION(rb_tensor,sub)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper)(obj, bb, TENSOR_SUB);
}

static VALUE FUNCTION(rb_tensor,mul_elements)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper)(obj, bb, TENSOR_MUL_ELEMENTS);
}

static VALUE FUNCTION(rb_tensor,div_elements)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper)(obj, bb, TENSOR_DIV_ELEMENTS);
}

static VALUE FUNCTION(rb_tensor,add_constant)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper)(obj, bb, TENSOR_ADD_CONSTANT);
}

static VALUE FUNCTION(rb_tensor,add_diagonal)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper)(obj, bb, TENSOR_ADD_DIAGONAL);
}

static VALUE FUNCTION(rb_tensor,scale)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper)(obj, bb, TENSOR_SCALE);
}

/***/
static VALUE FUNCTION(rb_tensor,add_bang)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper_bang)(obj, bb, TENSOR_ADD);
}

static VALUE FUNCTION(rb_tensor,sub_bang)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper_bang)(obj, bb, TENSOR_SUB);
}

static VALUE FUNCTION(rb_tensor,mul_elements_bang)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper_bang)(obj, bb, TENSOR_MUL_ELEMENTS);
}

static VALUE FUNCTION(rb_tensor,div_elements_bang)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper_bang)(obj, bb, TENSOR_DIV_ELEMENTS);
}

static VALUE FUNCTION(rb_tensor,add_constant_bang)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper_bang)(obj, bb, TENSOR_ADD_CONSTANT);
}

static VALUE FUNCTION(rb_tensor,add_diagonal_bang)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper_bang)(obj, bb, TENSOR_ADD_DIAGONAL);
}

static VALUE FUNCTION(rb_tensor,scale_bang)(VALUE obj, VALUE bb)
{
  return FUNCTION(rb_tensor,oper_bang)(obj, bb, TENSOR_SCALE);
}

/*****/
static VALUE FUNCTION(rb_tensor,product_singleton)(VALUE obj, VALUE aa, VALUE bb)
{
  GSL_TYPE(rbgsl_tensor) *a, *b, *c;
  CHECK_TEN(aa);
  CHECK_TEN(bb);
  Data_Get_Struct(aa, GSL_TYPE(rbgsl_tensor), a);
  switch (TYPE(bb)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return FUNCTION(rb_tensor,mul_elements)(aa, bb);
    break;
  default:
    Data_Get_Struct(bb, GSL_TYPE(rbgsl_tensor), b);
    c = ALLOC(GSL_TYPE(rbgsl_tensor));
    c->tensor = FUNCTION(tensor,product(a->tensor, b->tensor));
    if (c->tensor->rank == 0)
      c->indices = gsl_permutation_alloc(1);
    else
      c->indices = gsl_permutation_alloc(c->tensor->rank);
    return Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), c);
    break;
  }
}

static VALUE FUNCTION(rb_tensor,product)(VALUE obj, VALUE bb)
{
  GSL_TYPE(rbgsl_tensor) *a, *b, *c;
  switch (TYPE(bb)) {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    return FUNCTION(rb_tensor,mul_elements)(obj, bb);
    break;
  default:
    CHECK_TEN(bb);
    Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), a);
    Data_Get_Struct(bb, GSL_TYPE(rbgsl_tensor), b);
    c = ALLOC(GSL_TYPE(rbgsl_tensor));
    c->tensor = FUNCTION(tensor,product(a->tensor, b->tensor));
    if (c->tensor->rank == 0)
      c->indices = gsl_permutation_alloc(1);
    else
      c->indices = gsl_permutation_alloc(c->tensor->rank);
    return Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), c);
    break;
  }
}

static VALUE FUNCTION(rb_tensor,contract)(VALUE obj, VALUE ii, VALUE jj)
{
  GSL_TYPE(rbgsl_tensor) *t, *tnew;
  size_t rank;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  tnew = ALLOC(GSL_TYPE(rbgsl_tensor));
  tnew->tensor = FUNCTION(tensor,contract)(t->tensor, FIX2INT(ii), FIX2INT(jj));
  if (tnew->tensor->rank == 0) rank = 1;
  else rank = tnew->tensor->rank;
  tnew->indices = gsl_permutation_alloc(rank);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), tnew);
}

static VALUE FUNCTION(rb_tensor,size)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  return INT2FIX(t->tensor->size);
}

static VALUE FUNCTION(rb_tensor,rank)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  return INT2FIX(t->tensor->rank);
}

static VALUE FUNCTION(rb_tensor,dimension)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  return INT2FIX(t->tensor->dimension);
}

static VALUE FUNCTION(rb_tensor,data)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  QUALIFIED_VIEW(gsl_vector,view) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  v = FUNCTION(rb_gsl_make_vector,view)(t->tensor->data, t->tensor->size, 1);
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,view), 0, free, v);
}

static VALUE FUNCTION(rb_tensor,2matrix)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  GSL_TYPE(gsl_matrix) *m = NULL;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  m = (GSL_TYPE(gsl_matrix)*) FUNCTION(tensor,2matrix)(t->tensor);
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_matrix,view), 0, FUNCTION(gsl_matrix,free), m);
}

static VALUE FUNCTION(rb_tensor,2vector)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  GSL_TYPE(gsl_vector) *v = NULL;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  v = (GSL_TYPE(gsl_vector) *) FUNCTION(tensor,2vector)(t->tensor);
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_vector,view), 0, FUNCTION(gsl_vector,free), v);
}

static VALUE FUNCTION(rb_tensor,to_v)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  GSL_TYPE(gsl_vector) *v;

  v = FUNCTION(gsl_vector,alloc)(t->tensor->size);
  memcpy(v->data, t->tensor->data, sizeof(BASE)*v->size);
  return Data_Wrap_Struct(GSL_TYPE(cgsl_vector), 0, FUNCTION(gsl_vector,free), v);
}

/*
  Creates a subtensor slicing the existing tensor.
  NOTE: no new data region is malloced.
  t: Tensor
  rank: rank of the tensor created
*/
GSL_TYPE(tensor) FUNCTION(tensor,subtensor)(const GSL_TYPE(tensor) *t,
						     const unsigned int rank,
						     size_t *indices)
{
  GSL_TYPE(tensor) tnew;
  size_t position;
  tnew.rank = rank;
  tnew.dimension = t->dimension;
  tnew.size =  quick_pow(t->dimension, rank);
  position = FUNCTION(tensor,position)(indices, t);
  if (position >= t->size)
    rb_raise(rb_eRangeError, "wrong indices given");
  tnew.data = t->data + position;
  return tnew;
}

static VALUE FUNCTION(rb_tensor,subtensor)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t, *tnew;
  unsigned int rank;
  size_t n;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  /* n: number of indices given */
  rbgsl_tensor_get_indices(argc, argv, t->indices, &n);  
  rank = t->tensor->rank - n; 
  tnew = ALLOC(GSL_TYPE(rbgsl_tensor));
  tnew->tensor = (GSL_TYPE(tensor)*) malloc(sizeof(GSL_TYPE(tensor)));
  *(tnew->tensor) = FUNCTION(tensor,subtensor)(t->tensor, rank, t->indices->data);
  if (rank == 0) 
    tnew->indices = gsl_permutation_alloc(1);
  else
    tnew->indices = gsl_permutation_alloc(rank);
  return Data_Wrap_Struct(QUALIFIED_VIEW(cgsl_tensor,view), 0, FUNCTION(rbgsl_tensor,free2), tnew);	       
}

#ifdef BASE_DOUBLE
#define SHOW_ELM 6
#define PRINTF_FORMAT "%4.3e "
#else
#define SHOW_ELM 15
#define PRINTF_FORMAT "%d "
#endif
static VALUE FUNCTION(rb_tensor,to_s)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  QUALIFIED_VIEW(gsl_matrix,view) matrix;
  QUALIFIED_VIEW(gsl_vector,view) vector;
  GSL_TYPE(gsl_matrix) *m;
  GSL_TYPE(gsl_vector) *v;
  char buf[16];
  size_t i, j;
  VALUE str;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  str = rb_str_new2("[ ");
  switch (t->tensor->rank) {
  case 2:
    matrix.matrix.data = t->tensor->data;
    matrix.matrix.size1 = t->tensor->dimension;
    matrix.matrix.size2 = t->tensor->dimension;
    matrix.matrix.tda = t->tensor->dimension;
    matrix.matrix.block = 0;
    matrix.matrix.owner = 0;
    m = &(matrix.matrix);
    for (i = 0; i < m->size1; i++) {
      if (i != 0) {
	strcpy(buf, "  ");
	rb_str_cat(str, buf, strlen(buf));
      }
      for (j = 0; j < m->size2; j++) {
	sprintf(buf, PRINTF_FORMAT, FUNCTION(gsl_matrix,get)(m, i, j));
	rb_str_cat(str, buf, strlen(buf));
	if (j == SHOW_ELM) {
	  strcpy(buf, "... ");
	  rb_str_cat(str, buf, strlen(buf));
	  break;
	}
      }
      if (i == 6) {
	strcpy(buf, "\n  ... ]");
	rb_str_cat(str, buf, strlen(buf));
	break;
      }
      if (i == m->size1 - 1) {
	strcpy(buf, "]");
	rb_str_cat(str, buf, strlen(buf));
      } else {
	strcpy(buf, "\n");
	rb_str_cat(str, buf, strlen(buf));
      }
    }
    return str;
    break;
  default:
    vector.vector.data = t->tensor->data;
    vector.vector.stride = 1;
    vector.vector.size = t->tensor->size;
    vector.vector.owner = 0;
    vector.vector.block = 0;
    v = &(vector.vector);
    sprintf(buf,  PRINTF_FORMAT, FUNCTION(gsl_vector,get)(v, 0));
    rb_str_cat(str, buf, strlen(buf));
    for (i = 1; i < v->size; i++) {
      sprintf(buf,  PRINTF_FORMAT, FUNCTION(gsl_vector,get)(v, i));
      rb_str_cat(str, buf, strlen(buf));
      if (i == SHOW_ELM && i != v->size-1) {
        strcpy(buf, "... ");
        rb_str_cat(str, buf, strlen(buf));
        break;
      }
    }
    sprintf(buf, "]");
    rb_str_cat(str, buf, strlen(buf));
    return str;
    break;
  }
}
#undef SHOW_ELM
#undef PRINTF_FORMAT

static VALUE FUNCTION(rb_tensor,inspect)(VALUE obj)
{
  VALUE str;
  char buf[64];
  sprintf(buf, "%s\n", rb_class2name(CLASS_OF(obj)));
  str = rb_str_new2(buf);
  return rb_str_concat(str, FUNCTION(rb_tensor,to_s)(obj));
}

VALUE FUNCTION(rb_tensor,equal)(int argc, VALUE *argv, VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *a, *b;
  GSL_TYPE(gsl_matrix) *m;
  GSL_TYPE(gsl_vector) *v;
  VALUE other;
  double eps = 1e-10;
  size_t i;
  switch (argc) {
  case 2:
    other = argv[0];
    eps = NUM2DBL(argv[1]);
    break;
  case 1:
    other = argv[0];
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  }
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), a);
  if (TEN_P(other)) {
    Data_Get_Struct(other, GSL_TYPE(rbgsl_tensor), b);
    if (a->tensor->rank != b->tensor->rank) return Qfalse;
    if (a->tensor->dimension != b->tensor->dimension) return Qfalse;
    if (a->tensor->size != b->tensor->size) return Qfalse;
    for (i = 0; i < a->tensor->size; i++)
      if (fabs(a->tensor->data[i]-b->tensor->data[i]) > eps) 
	return Qfalse;
    return Qtrue;
  } else if (MAT_P(other)) {
    if (a->tensor->rank != 2) return Qfalse;
    Data_Get_Struct(other, GSL_TYPE(gsl_matrix), m);
    if (a->tensor->dimension != m->size1 || a->tensor->dimension != m->size2)
      return Qfalse;
    for (i = 0; i < a->tensor->size; i++)
      if (fabs(a->tensor->data[i]-m->data[i]) > eps) 
	return Qfalse;
    return Qtrue;
  } else if (VEC_P(other)) {
    Data_Get_Struct(other, GSL_TYPE(gsl_vector), v);
    if (a->tensor->size != v->size) return Qfalse;
    for (i = 0; i < a->tensor->size; i++)
      if (fabs(a->tensor->data[i]-v->data[i]) > eps) 
	return Qfalse;
    return Qtrue;
  } else {
    rb_raise(rb_eTypeError, "wrong argument type %s (Tensor, Matrix or Vector expected)", 
	     rb_class2name(CLASS_OF(other)));
  }
}

static VALUE FUNCTION(rb_tensor,uplus)(VALUE obj)
{
  return obj;
}

static VALUE FUNCTION(rb_tensor,uminus)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t, *tnew;
  size_t i;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  tnew = FUNCTION(rbgsl_tensor,copy)(t);
  for (i = 0; i < tnew->tensor->size; i++)
    tnew->tensor->data[i] *= -1;
  return Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), tnew);
}

static VALUE FUNCTION(rb_tensor,coerce)(VALUE obj, VALUE other)
{
  GSL_TYPE(rbgsl_tensor) *t, *tnew;
  VALUE tt;
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  switch (TYPE(other)) {
  case T_FLOAT:
  case T_FIXNUM:
  case T_BIGNUM:
    tnew = FUNCTION(rbgsl_tensor,alloc)(t->tensor->rank, t->tensor->dimension);
    FUNCTION(tensor,set_all)(tnew->tensor, NUMCONV(other));
    tt = Data_Wrap_Struct(GSL_TYPE(cgsl_tensor), 0, FUNCTION(rbgsl_tensor,free), tnew);
    return rb_ary_new3(2, tt, obj);
    break;
  default:
    rb_raise(rb_eRuntimeError, "undefined operation with %s", 
	     rb_class2name(CLASS_OF(other)));
    break;
  }
}

static VALUE FUNCTION(rb_tensor,info)(VALUE obj)
{
  GSL_TYPE(rbgsl_tensor) *t;
  char buf[256];
  Data_Get_Struct(obj, GSL_TYPE(rbgsl_tensor), t);
  sprintf(buf, "Class:      %s\n", rb_class2name(CLASS_OF(obj)));
  sprintf(buf, "%sSuperClass: %s\n", buf, rb_class2name(RCLASS_SUPER(CLASS_OF(obj))));
  sprintf(buf, "%sRank:       %d\n", buf, (int) t->tensor->rank);
  sprintf(buf, "%sDimension:  %d\n", buf, (int) t->tensor->dimension);
  sprintf(buf, "%sSize:       %d\n", buf, (int) t->tensor->size);
  return rb_str_new2(buf);
}

void FUNCTION(Init_tensor,init)(VALUE module)
{
#ifdef BASE_DOUBLE
  cgsl_tensor = rb_define_class_under(module, "Tensor", cGSL_Object);
  cgsl_tensor_int = rb_define_class_under(cgsl_tensor, "Int", cGSL_Object);
  cgsl_tensor_view = rb_define_class_under(cgsl_tensor, "View", cgsl_tensor);
  cgsl_tensor_int_view = rb_define_class_under(cgsl_tensor_int, "View",
					       cgsl_tensor_int);
  /*
  cgsl_index = rb_define_class_under(cgsl_tensor, "Index", 
  cgsl_permutation);*/
#endif

  rb_define_singleton_method(GSL_TYPE(cgsl_tensor), "new",
			     FUNCTION(rb_tensor,new), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_tensor), "[]",
			     FUNCTION(rb_tensor,new), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_tensor), "alloc",
			     FUNCTION(rb_tensor,new), -1);
  rb_define_singleton_method(GSL_TYPE(cgsl_tensor), "calloc",
			     FUNCTION(rb_tensor,calloc), 2);
  rb_define_singleton_method(GSL_TYPE(cgsl_tensor), "copy",
			     FUNCTION(rb_tensor,copy_singleton), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_tensor), "memcpy",
			     FUNCTION(rb_tensor,memcpy_singleton), 2);
  rb_define_singleton_method(GSL_TYPE(cgsl_tensor), "swap",
			     FUNCTION(rb_tensor,swap_singleton), 2);

  /*****/

  rb_define_method(GSL_TYPE(cgsl_tensor), "copy",
			     FUNCTION(rb_tensor,copy), 0);
  rb_define_alias(GSL_TYPE(cgsl_tensor), "clone", "copy");
  rb_define_alias(GSL_TYPE(cgsl_tensor), "duplicate", "copy");
  rb_define_method(GSL_TYPE(cgsl_tensor), "set_zero",
			     FUNCTION(rb_tensor,set_zero), 0);
  rb_define_method(GSL_TYPE(cgsl_tensor), "set_all",
			     FUNCTION(rb_tensor,set_all), 1);

  rb_define_method(GSL_TYPE(cgsl_tensor), "position",
		   FUNCTION(rb_tensor,position), -1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "get",
			     FUNCTION(rb_tensor,get), -1);
  rb_define_alias(GSL_TYPE(cgsl_tensor), "[]", "get");
  rb_define_method(GSL_TYPE(cgsl_tensor), "set",
			     FUNCTION(rb_tensor,set), -1);
  rb_define_alias(GSL_TYPE(cgsl_tensor), "[]=", "set");

  rb_define_method(GSL_TYPE(cgsl_tensor), "fread",
			     FUNCTION(rb_tensor,fread), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "fwrite",
			     FUNCTION(rb_tensor,fwrite), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "fprintf",
			     FUNCTION(rb_tensor,fprintf), -1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "printf",
			     FUNCTION(rb_tensor,printf), -1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "fscanf",
			     FUNCTION(rb_tensor,fscanf), 1);

  rb_define_method(GSL_TYPE(cgsl_tensor), "swap_indices",
			     FUNCTION(rb_tensor,swap_indices), 2);

  rb_define_method(GSL_TYPE(cgsl_tensor), "max",
			     FUNCTION(rb_tensor,max), 0);
  rb_define_method(GSL_TYPE(cgsl_tensor), "min",
			     FUNCTION(rb_tensor,min), 0);
  rb_define_method(GSL_TYPE(cgsl_tensor), "minmax",
			     FUNCTION(rb_tensor,minmax), 0);
  rb_define_method(GSL_TYPE(cgsl_tensor), "max_index",
			     FUNCTION(rb_tensor,max_index), 0);
  rb_define_method(GSL_TYPE(cgsl_tensor), "min_index",
			     FUNCTION(rb_tensor,min_index), 0);
  rb_define_method(GSL_TYPE(cgsl_tensor), "minmax_index",
			     FUNCTION(rb_tensor,minmax_index), 0);

  rb_define_method(GSL_TYPE(cgsl_tensor), "isnull",
			     FUNCTION(rb_tensor,isnull), 0);
  rb_define_method(GSL_TYPE(cgsl_tensor), "isnull?",
			     FUNCTION(rb_tensor,isnull2), 0);

  rb_define_method(GSL_TYPE(cgsl_tensor), "add",
			     FUNCTION(rb_tensor,add), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "sub",
			     FUNCTION(rb_tensor,sub), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "mul_elements",
			     FUNCTION(rb_tensor,mul_elements), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "div_elements",
			     FUNCTION(rb_tensor,div_elements), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "add_constant",
			     FUNCTION(rb_tensor,add_constant), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "add_diagonal",
			     FUNCTION(rb_tensor,add_diagonal), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "scale",
			     FUNCTION(rb_tensor,scale), 1);
  rb_define_singleton_method(GSL_TYPE(cgsl_tensor), "product",
			     FUNCTION(rb_tensor,product_singleton), 2);
  rb_define_method(GSL_TYPE(cgsl_tensor), "product",
			     FUNCTION(rb_tensor,product), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "contract",
			     FUNCTION(rb_tensor,contract), 2);

  rb_define_alias(GSL_TYPE(cgsl_tensor), "+", "add");
  rb_define_alias(GSL_TYPE(cgsl_tensor), "-", "sub");
  /*  rb_define_alias(GSL_TYPE(cgsl_tensor), "*", "mul_elements");*/
  rb_define_alias(GSL_TYPE(cgsl_tensor), "/", "div_elements");
  rb_define_alias(GSL_TYPE(cgsl_tensor), "*", "product");

  rb_define_method(GSL_TYPE(cgsl_tensor), "add!",
			     FUNCTION(rb_tensor,add_bang), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "sub!",
			     FUNCTION(rb_tensor,sub_bang), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "mul_elements!",
			     FUNCTION(rb_tensor,mul_elements_bang), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "div_elements!",
			     FUNCTION(rb_tensor,div_elements_bang), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "add_constant!",
			     FUNCTION(rb_tensor,add_constant_bang), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "add_diagonal!",
			     FUNCTION(rb_tensor,add_diagonal_bang), 1);
  rb_define_method(GSL_TYPE(cgsl_tensor), "scale!",
			     FUNCTION(rb_tensor,scale_bang), 1);

  rb_define_method(GSL_TYPE(cgsl_tensor), "+@",
			     FUNCTION(rb_tensor,uplus), 0);
  rb_define_method(GSL_TYPE(cgsl_tensor), "-@",
			     FUNCTION(rb_tensor,uminus), 0);

  rb_define_method(GSL_TYPE(cgsl_tensor), "size",
			     FUNCTION(rb_tensor,size), 0);
  rb_define_method(GSL_TYPE(cgsl_tensor), "rank",
			     FUNCTION(rb_tensor,rank), 0);
  rb_define_method(GSL_TYPE(cgsl_tensor), "dimension",
			     FUNCTION(rb_tensor,dimension), 0);
  rb_define_alias(GSL_TYPE(cgsl_tensor), "dim", "dimension");
  rb_define_method(GSL_TYPE(cgsl_tensor), "data",
			     FUNCTION(rb_tensor,data), 0);

  rb_define_method(GSL_TYPE(cgsl_tensor), "to_v",
			     FUNCTION(rb_tensor,to_v), 0);
  rb_define_alias(GSL_TYPE(cgsl_tensor), "to_gv", "to_v");

  rb_define_method(GSL_TYPE(cgsl_tensor), "to_vector",
			     FUNCTION(rb_tensor,2vector), 0);
  rb_define_method(GSL_TYPE(cgsl_tensor), "to_matrix",
			     FUNCTION(rb_tensor,2matrix), 0);

  rb_define_method(GSL_TYPE(cgsl_tensor), "subtensor",
			     FUNCTION(rb_tensor,subtensor), -1);
  rb_define_alias(GSL_TYPE(cgsl_tensor), "view", "subtensor");

  rb_define_method(GSL_TYPE(cgsl_tensor), "to_s",
			     FUNCTION(rb_tensor,to_s), 0);
  rb_define_method(GSL_TYPE(cgsl_tensor), "inspect",
			     FUNCTION(rb_tensor,inspect), 0);

  rb_define_method(GSL_TYPE(cgsl_tensor), "equal?",
			     FUNCTION(rb_tensor,equal), -1);
  rb_define_alias(GSL_TYPE(cgsl_tensor), "==", "equal?");

  rb_define_method(GSL_TYPE(cgsl_tensor), "coerce",
			     FUNCTION(rb_tensor,coerce), 1);

  rb_define_method(GSL_TYPE(cgsl_tensor), "info", 
		   FUNCTION(rb_tensor,info), 0);
}

#undef NUMCONV
#undef C_TO_VALUE
#undef CHECK_TEN
#undef TEN_P
#undef VEC_P
#undef MAT_P

#endif
