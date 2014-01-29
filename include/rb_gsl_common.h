/*
  rb_gsl_common.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/
#ifndef ___RB_GSL_COMMON_H___
#define ___RB_GSL_COMMON_H___

#include "rb_gsl_config.h"

#include <ruby.h>
#ifdef HAVE_RUBY_IO_H
#include <ruby/io.h>
#else
#include <rubyio.h>
#endif

#include <ctype.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_version.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_ieee_utils.h>
#include "rb_gsl_with_narray.h"

EXTERN ID rb_gsl_id_beg, rb_gsl_id_end, rb_gsl_id_excl, rb_gsl_id_to_a;

#ifndef CHECK_FIXNUM
#define CHECK_FIXNUM(x) if(!FIXNUM_P(x))rb_raise(rb_eTypeError,"Fixnum expected");
#endif

#ifndef Need_Float
#define Need_Float(x) (x) = rb_Float(x)
#endif

#ifndef Need_Float2
#define Need_Float2(x,y) do {\
    Need_Float(x);\
    Need_Float(y);} while (0)
#endif

#ifndef COMPLEX_P
#define COMPLEX_P(x) (rb_obj_is_kind_of(x,cgsl_complex))
#endif

#ifndef CHECK_RNG
#define CHECK_RNG(x) if(!rb_obj_is_kind_of(x,cgsl_rng))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Rng expected)");
#endif

#ifndef CHECK_COMPLEX
#define CHECK_COMPLEX(x) if(!rb_obj_is_kind_of(x,cgsl_complex))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Complex expected)");
#endif

#ifndef POLY_P
#define POLY_P(x) (rb_obj_is_kind_of(x,cgsl_poly))
#endif

#ifndef CHECK_POLY
#define CHECK_POLY(x) if(!rb_obj_is_kind_of(x,cgsl_poly))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Poly expected)");
#endif

#ifndef CHECK_POLY_INT
#define CHECK_POLY_INT(x) if(!rb_obj_is_kind_of(x,cgsl_poly_int))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Poly::Int expected)");
#endif

#ifndef RATIONAL_P
#define RATIONAL_P(x) (rb_obj_is_kind_of(x,cgsl_rational))
#endif

#ifndef CHECK_RATIONAL
#define CHECK_RATIONAL(x) if(!rb_obj_is_kind_of(x,cgsl_rational))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Rational expected)");
#endif

/*****/
#ifndef BLOCK_P
#define BLOCK_P(x) (rb_obj_is_kind_of(x,cgsl_block))
#endif

#ifndef CHECK_BLOCK
#define CHECK_BLOCK(x) if(!rb_obj_is_kind_of(x,cgsl_block))\
    rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Block expected)", rb_class2name(CLASS_OF(x)));
#endif
#ifndef BLOCK_INT_P
#define BLOCK_INT_P(x) (rb_obj_is_kind_of(x,cgsl_block_int))
#endif

#ifndef CHECK_BLOCK_INT
#define CHECK_BLOCK_INT(x) if(!rb_obj_is_kind_of(x,cgsl_block_int))\
    rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Block::Int expected)", rb_class2name(CLASS_OF(x)));
#endif
#ifndef BLOCK_UCHAR_P
#define BLOCK_UCHAR_P(x) (rb_obj_is_kind_of(x,cgsl_block_uchar))
#endif

#ifndef CHECK_BLOCK_UCHAR
#define CHECK_BLOCK_UCHAR(x) if(!rb_obj_is_kind_of(x,cgsl_block_uchar))\
    rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Block::Byte expected)", rb_class2name(CLASS_OF(x)));
#endif
/*****/

#ifndef VECTOR_P
#define VECTOR_P(x) (rb_obj_is_kind_of(x,cgsl_vector))
#endif

#ifndef VECTOR_VIEW_P
#define VECTOR_VIEW_P(x) ((CLASS_OF(x)==cgsl_vector_view||CLASS_OF(x)==cgsl_vector_col_view||CLASS_OF(x)==cgsl_vector_view_ro||CLASS_OF(x)==cgsl_vector_col_view_ro))
#endif

#ifndef VECTOR_ROW_P
#define VECTOR_ROW_P(x) ((CLASS_OF(x)==cgsl_vector||CLASS_OF(x)==cgsl_vector_view||CLASS_OF(x)==cgsl_vector_view_ro))
#endif

#ifndef VECTOR_COL_P
#define VECTOR_COL_P(x) ((CLASS_OF(x)==cgsl_vector_col||CLASS_OF(x)==cgsl_vector_col_view||CLASS_OF(x)==cgsl_vector_col_view_ro))
#endif

#ifndef VECTOR_ROW_COL
#define VECTOR_ROW_COL(x) ((rb_obj_is_kind_of(x,cgsl_vector_col)||rb_obj_is_kind_of(x,cgsl_vector_int_col))?cgsl_vector_col:cgsl_vector)
#endif

#ifndef CHECK_VECTOR
#define CHECK_VECTOR(x) if(!rb_obj_is_kind_of(x,cgsl_vector))\
    rb_raise(rb_eTypeError, "wrong argument type %s (GSL::Vector expected)", rb_class2name(CLASS_OF(x)));
#endif

#ifdef HAVE_NARRAY_H
#define Data_Get_Vector(obj,sval) do {\
    if (NA_IsNArray(obj)) {\
      /* Convert obj to GSL::Vector::View */\
      obj = rb_gsl_na_to_gsl_vector_view_method(obj);\
    }\
    CHECK_VECTOR(obj);\
    Data_Get_Struct(obj,gsl_vector,sval);\
} while (0)
#else
#define Data_Get_Vector(obj,sval) do {\
    CHECK_VECTOR(obj);\
    Data_Get_Struct(obj,gsl_vector,sval);\
} while (0)
#endif

/******/
#ifndef VECTOR_INT_P
#define VECTOR_INT_P(x) (rb_obj_is_kind_of(x,cgsl_vector_int))
#endif

#ifndef VECTOR_INT_VIEW_P
#define VECTOR_INT_VIEW_P(x) ((CLASS_OF(x)==cgsl_vector_int_view||CLASS_OF(x)==cgsl_vector_int_col_view||CLASS_OF(x)==cgsl_vector_int_view_ro||CLASS_OF(x)==cgsl_vector_int_col_view_ro))
#endif

#ifndef VECTOR_INT_ROW_P
#define VECTOR_INT_ROW_P(x) ((CLASS_OF(x)==cgsl_vector_int||CLASS_OF(x)==cgsl_vector_int_view||CLASS_OF(x)==cgsl_vector_int_view_ro))
#endif

#ifndef VECTOR_INT_COL_P
#define VECTOR_INT_COL_P(x) ((CLASS_OF(x)==cgsl_vector_int_col||CLASS_OF(x)==cgsl_vector_int_col_view||CLASS_OF(x)==cgsl_vector_int_col_view_ro))
#endif

#ifndef VECTOR_INT_ROW_COL
#define VECTOR_INT_ROW_COL(x) (VECTOR_INT_ROW_P(x)?cgsl_vector_int:cgsl_vector_int_col)
#endif

#ifndef CHECK_VECTOR_INT
#define CHECK_VECTOR_INT(x) if(!rb_obj_is_kind_of(x,cgsl_vector_int))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Vector::Int expected)");
#endif

/******/
#ifndef VECTOR_COMPLEX_P
#define VECTOR_COMPLEX_P(x) (rb_obj_is_kind_of(x,cgsl_vector_complex))
#endif

#ifndef VECTOR_COMPLEX_ROW_P
#define VECTOR_COMPLEX_ROW_P(x) ((CLASS_OF(x)==cgsl_vector_complex||CLASS_OF(x)==cgsl_vector_complex_view))
#endif

#ifndef VECTOR_COMPLEX_COL_P
#define VECTOR_COMPLEX_COL_P(x) ((CLASS_OF(x)==cgsl_vector_complex_col||CLASS_OF(x)==cgsl_vector_complex_col_view))
#endif

#ifndef VECTOR_COMPLEX_ROW_COL
#define VECTOR_COMPLEX_ROW_COL(x) (VECTOR_COMPLEX_ROW_P(x)?cgsl_vector_complex:cgsl_vector_complex_col)
#endif

#ifndef CHECK_VECTOR_COMPLEX
#define CHECK_VECTOR_COMPLEX(x) if(!rb_obj_is_kind_of(x,cgsl_vector_complex))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Vector::Complex expected)");
#endif

#ifndef MATRIX_P
#define MATRIX_P(x) (rb_obj_is_kind_of(x,cgsl_matrix))
#endif

#ifndef CHECK_MATRIX
#define CHECK_MATRIX(x) if(!rb_obj_is_kind_of(x,cgsl_matrix))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Matrix expected)");
#endif

#define Data_Get_Matrix(obj,sval) do {\
    CHECK_MATRIX(obj);\
    Data_Get_Struct(obj,gsl_matrix,sval);\
} while (0)


#ifndef MATRIX_INT_P
#define MATRIX_INT_P(x) (rb_obj_is_kind_of(x,cgsl_matrix_int))
#endif

#ifndef CHECK_MATRIX_INT
#define CHECK_MATRIX_INT(x) if(!rb_obj_is_kind_of(x,cgsl_matrix_int))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Matrix::Int expected)");
#endif

#ifndef MATRIX_COMPLEX_P
#define MATRIX_COMPLEX_P(x) (rb_obj_is_kind_of(x,cgsl_matrix_complex))
#endif

#ifndef CHECK_MATRIX_COMPLEX
#define CHECK_MATRIX_COMPLEX(x) if(!rb_obj_is_kind_of(x,cgsl_matrix_complex))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Matrix::Complex expected)");
#endif

#ifndef TENSOR_P
#define TENSOR_P(x) ((CLASS_OF(x)==cgsl_tensor))
#endif

#ifndef CHECK_TENSOR
#define CHECK_TENSOR(x) if(CLASS_OF(x)!=cgsl_tensor)\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Tensor expected)");
#endif

#ifndef TENSOR_INT_P
#define TENSOR_INT_P(x) ((CLASS_OF(x)==cgsl_tensor_int))
#endif

#ifndef CHECK_TENSOR_INT
#define CHECK_TENSOR_INT(x) if(CLASS_OF(x)!=cgsl_tensor_int)\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Tensor::Int expected)");
#endif

#ifndef PERMUTATION_P
#define PERMUTATION_P(x) (rb_obj_is_kind_of(x,cgsl_permutation))
#endif

#ifndef CHECK_PERMUTATION
#define CHECK_PERMUTATION(x) if(!rb_obj_is_kind_of(x,cgsl_permutation))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Permutation expected)");
#endif

#ifndef PROC_P
#define PROC_P(x) (rb_obj_is_kind_of(x,rb_cProc))
#endif

#ifndef CHECK_PROC
#define CHECK_PROC(x) if(!rb_obj_is_kind_of(x,rb_cProc))\
    rb_raise(rb_eTypeError, "wrong argument type (Proc expected)");
#endif

#ifndef FUNCTION_P
#define FUNCTION_P(x) (rb_obj_is_kind_of(x,cgsl_function))
#endif

#ifndef CHECK_FUNCTION
#define CHECK_FUNCTION(x) if(!rb_obj_is_kind_of(x,cgsl_function))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Function expected)");
#endif

#ifndef FUNCTION_FDF_P
#define FUNCTION_FDF_P(x) (rb_obj_is_kind_of(x,cgsl_function_fdf))
#endif

#ifndef CHECK_FUNCTION_FDF
#define CHECK_FUNCTION_FDF(x) if(!rb_obj_is_kind_of(x,cgsl_function_fdf))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Function_fdf expected)");
#endif

#ifndef HISTOGRAM_P
#define HISTOGRAM_P(x) (rb_obj_is_kind_of(x,cgsl_histogram))
#endif

#ifndef CHECK_HISTOGRAM
#define CHECK_HISTOGRAM(x) if(!rb_obj_is_kind_of(x,cgsl_histogram))\
    rb_raise(rb_eTypeError, "wrong argument type (GSL::Histogram expected)");
#endif

#ifndef RBGSL_SET_CLASS
#ifdef RB_OBJ_WRITE
#define RBGSL_SET_CLASS0(obj0, cls) RB_OBJ_WRITE(obj0, &(RBASIC_CLASS(obj0)), cls)
#else
#define RBGSL_SET_CLASS0(obj0, cls) RBASIC(obj0)->klass = cls
#endif
#define RBGSL_SET_CLASS(obj, cls) do { \
  VALUE _obj_ = (obj); \
  RBGSL_SET_CLASS0(_obj_, cls); \
} while (0)
#endif

void rb_gsl_error_handler(const char *reason, const char *file,
				 int line, int gsl_errno);

FILE* rb_gsl_open_writefile(VALUE io, int *flag);
FILE* rb_gsl_open_readfile(VALUE io, int *flag);

VALUE rb_gsl_obj_read_only(int argc, VALUE *argv, VALUE obj);

int str_tail_grep(const char *s0, const char *s1);
int str_head_grep(const char *s0, const char *s1);

#ifndef STR2CSTR
#define STR2CSTR StringValuePtr
#endif

#ifndef RCLASS_SUPER
#define RCLASS_SUPER(cls) RCLASS(cls)->super
#endif

void make_graphcommand(char *command, VALUE hash);
int rbgsl_complex_equal(const gsl_complex *z1, const gsl_complex *z2, double eps);

gsl_vector* mygsl_vector_down(gsl_vector *p);
void mygsl_vector_up2(gsl_vector *pnew, gsl_vector *p);
gsl_vector* mygsl_vector_up(gsl_vector *p);

gsl_vector_int* mygsl_vector_int_down(gsl_vector_int *p);
void mygsl_vector_int_up2(gsl_vector_int *pnew, gsl_vector_int *p);
gsl_vector_int* mygsl_vector_int_up(gsl_vector_int *p);

#ifndef GSL_1_3_LATER
int gsl_fcmp(const double x1, const double x2, const double epsilon);
#endif

size_t count_columns(const char *str);
char* str_scan_double(const char *str, double *val);
char* str_scan_int(const char *str, int *val);
double* get_ptr_double3(VALUE obj, size_t *size, size_t *stride, int *flag);
gsl_complex ary2complex(VALUE obj);
VALUE vector_eval_create(VALUE obj, double (*func)(double));
VALUE matrix_eval_create(VALUE obj, double (*func)(double));
VALUE rb_gsl_ary_eval1(VALUE ary, double (*f)(double));
#ifdef HAVE_NARRAY_H
VALUE rb_gsl_nary_eval1(VALUE ary, double (*f)(double));
#endif

EXTERN VALUE cGSL_Object;
#endif
