/*
  rb_gsl_array.h
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2001-2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY
*/

#ifndef ___RB_GSL_ARRAY_H___
#define ___RB_GSL_ARRAY_H___

#include <math.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_blas.h>
#include "rb_gsl_common.h"

typedef gsl_permutation gsl_index;

#ifdef HAVE_NARRAY_H
EXTERN VALUE cNArray;
#endif

EXTERN VALUE cgsl_block, cgsl_block_int;
EXTERN VALUE cgsl_block_uchar;
EXTERN VALUE cgsl_block_complex;
EXTERN VALUE cgsl_vector, cgsl_vector_complex;
EXTERN VALUE cgsl_vector_col;
EXTERN VALUE cgsl_vector_col_view;
EXTERN VALUE cgsl_vector_complex_col;
EXTERN VALUE cgsl_vector_complex_col_view;
EXTERN VALUE cgsl_vector_view, cgsl_vector_complex_view;
EXTERN VALUE cgsl_vector_view_ro, cgsl_vector_col_view_ro;
EXTERN VALUE cgsl_vector_complex_view_ro;

EXTERN VALUE cgsl_vector_int, cgsl_vector_int_col;
EXTERN VALUE cgsl_vector_int_view, cgsl_vector_int_col_view;
EXTERN VALUE cgsl_vector_int_view_ro, cgsl_vector_int_col_view_ro;

EXTERN VALUE cgsl_matrix, cgsl_matrix_complex;
EXTERN VALUE cgsl_matrix_view_ro;
EXTERN VALUE cgsl_matrix_complex_view_ro;
EXTERN VALUE cgsl_matrix_view, cgsl_matrix_complex_view;
EXTERN VALUE cgsl_matrix_int, cgsl_matrix_int_view;
EXTERN VALUE cgsl_matrix_int_view_ro;
EXTERN VALUE cgsl_permutation;
EXTERN VALUE cgsl_index;
EXTERN VALUE cgsl_function;
EXTERN VALUE mgsl_narray;

EXTERN VALUE mDirac;

gsl_matrix_view* gsl_matrix_view_alloc();
void gsl_matrix_view_free(gsl_matrix_view * mv);
gsl_vector_view* gsl_vector_view_alloc();
void gsl_vector_view_free(gsl_vector_view * v);
gsl_vector_complex_view* gsl_vector_complex_view_alloc();
void gsl_vector_complex_view_free(gsl_vector_view * vv);
gsl_matrix_complex_view* gsl_matrix_complex_view_alloc();
void gsl_matrix_complex_view_free(gsl_matrix_view * vv);

VALUE rb_gsl_vector_new(int argc, VALUE *argv, VALUE klass);

gsl_vector* get_cvector(VALUE v);
VALUE make_rarray_from_cvector(const gsl_vector *v);
VALUE make_rarray_from_cpermutation(const gsl_permutation *v);
gsl_vector* make_cvector_from_rarrays(VALUE a);
gsl_vector* make_cvector_from_rarray(VALUE a);
void cvector_set_from_carray(gsl_vector *v, const double *a);
void cvector_set_from_rarrays(gsl_vector *v, VALUE ary);
void cvector_set_from_rarray(gsl_vector *v, VALUE ary);
void carray_set_from_cvector(double *a, const gsl_vector *v);
void carray_set_from_rarrays(double *a, VALUE ary);
void carray_set_from_rarray(double *a, VALUE ary);
int is_vector_p(VALUE obj);
void check_vector(VALUE obj);
int is_vector_complex_p(VALUE obj);
void check_vector_complex(VALUE obj);
int is_matrix_p(VALUE obj);
void check_matrix(VALUE obj);
int is_matrix_complex_p(VALUE obj);
void check_matrix_complex(VALUE obj);
gsl_complex* make_complex(double re, double im);
int is_permutation_p(VALUE obj);
void check_permutation(VALUE obj);
int is_combination_p(VALUE obj);
void check_combination(VALUE obj);
gsl_vector* get_vector(VALUE ary);
gsl_matrix* make_matrix_clone(const gsl_matrix *m);
gsl_matrix_int* make_matrix_int_clone(const gsl_matrix_int *m);
VALUE make_matrix_clone2(VALUE vm);
gsl_matrix_complex* make_matrix_complex_clone(const gsl_matrix_complex *m);
int is_matrix_complex_p(VALUE obj);
void check_matrix_complex(VALUE obj);

double* get_vector_ptr(VALUE ary, size_t *stride, size_t *n);


gsl_matrix_complex* matrix_to_complex(const gsl_matrix *m);

void gsl_matrix_complex_mul(gsl_matrix_complex *mnew, const gsl_matrix_complex *m,
          const gsl_matrix_complex *mb);
void gsl_matrix_mul(gsl_matrix *mnew, gsl_matrix *m, gsl_matrix *b);
void gsl_matrix_complex_mul_vector(gsl_vector_complex *vnew,
           const gsl_matrix_complex *m,
           const gsl_vector_complex *v);
void gsl_matrix_mul_vector(gsl_vector *vnew,
         const gsl_matrix *m, const gsl_vector *v);
gsl_vector_complex* vector_to_complex(const gsl_vector *v);

gsl_vector* make_vector_clone(const gsl_vector *v);
gsl_vector_complex* make_vector_complex_clone(const gsl_vector_complex *v);
int gsl_vector_complex_add(gsl_vector_complex *cv, const gsl_vector_complex *cv2);
int gsl_vector_complex_sub(gsl_vector_complex *cv, const gsl_vector_complex *cv2);
int gsl_vector_complex_mul(gsl_vector_complex *cv, const gsl_vector_complex *cv2);
int gsl_vector_complex_div(gsl_vector_complex *cv, const gsl_vector_complex *cv2);
int gsl_vector_complex_add_constant(gsl_vector_complex *cv, gsl_complex b);
int gsl_vector_complex_scale(gsl_vector_complex *cv, gsl_complex b);
gsl_vector_view* rb_gsl_make_vector_view(double *data, size_t size, size_t stride);
gsl_vector_int_view* rb_gsl_make_vector_int_view(int *data, size_t size, size_t stride);

void Init_gsl_array_complex(VALUE module);
void Init_gsl_vector(VALUE module);
void Init_gsl_vector_complex(VALUE module);
void Init_gsl_matrix(VALUE module);
void Init_gsl_matrix_complex(VALUE module);
void Init_gsl_matrix(VALUE module);
void Init_gsl_permutation(VALUE module);
void Init_gsl_combination(VALUE module);
void Init_gsl_matrix_int(VALUE module);

VALUE rb_gsl_range2ary(VALUE obj);

void Init_gsl_vector_int();
gsl_vector_int_view* rb_gsl_vector_int_view_alloc(size_t n);
void rb_gsl_vector_int_view_free(gsl_vector_int_view *v);
gsl_vector_int* make_vector_int_clone(const gsl_vector_int *v);
gsl_matrix_int_view* rb_gsl_matrix_int_view_alloc();
void rb_gsl_matrix_int_view_free(gsl_matrix_int_view *v);
VALUE rb_gsl_matrix_to_i(VALUE obj);
VALUE rb_gsl_matrix_int_to_f(VALUE obj);
void gsl_matrix_int_mul_vector(gsl_vector_int *vnew,
             const gsl_matrix_int *m, const gsl_vector_int *v);
VALUE rb_gsl_vector_to_i(VALUE obj);
VALUE make_rarray_from_cvector_int(const gsl_vector_int *v);
VALUE rb_gsl_vector_int_to_f(VALUE obj);
VALUE rb_gsl_vector_uminus(VALUE obj);
VALUE rb_gsl_vector_print(VALUE obj);
void gsl_vector_print(const gsl_vector *v, VALUE klass);
int rbgsl_vector_equal(const gsl_vector *v1, const gsl_vector *v2, double eps);

#ifndef GSL_1_2_LATER
int gsl_matrix_complex_add(gsl_matrix_complex * a, const gsl_matrix_complex * b);
int gsl_matrix_complex_sub(gsl_matrix_complex * a, const gsl_matrix_complex * b);
int gsl_matrix_complex_mul_elements(gsl_matrix_complex * a, const gsl_matrix_complex * b);
int gsl_matrix_complex_div_elements(gsl_matrix_complex * a, const gsl_matrix_complex * b);
int gsl_matrix_complex_scale(gsl_matrix_complex * a, const gsl_complex x);
int gsl_matrix_complex_add_constant(gsl_matrix_complex * a, const gsl_complex x);
int gsl_matrix_complex_add_diagonal(gsl_matrix_complex * a, const gsl_complex x);
#endif

void Init_gsl_vector_init(VALUE module);
void Init_gsl_vector_int_init(VALUE module);
void Init_gsl_matrix_init(VALUE module);
void Init_gsl_matrix_int_init(VALUE module);

gsl_matrix* gsl_matrix_alloc_from_array_sizes(VALUE ary,
                 VALUE nn1, VALUE nn2);
gsl_matrix* gsl_matrix_alloc_from_arrays(int argc, VALUE *argv);

gsl_matrix* gsl_matrix_alloc_from_vector_sizes(VALUE ary,
                VALUE nn1, VALUE nn2);
gsl_matrix* gsl_matrix_alloc_from_vectors(int argc, VALUE *argv);
gsl_matrix_int* gsl_matrix_int_alloc_from_array_sizes(VALUE ary,
                 VALUE nn1, VALUE nn2);
gsl_matrix_int* gsl_matrix_int_alloc_from_arrays(int argc, VALUE *argv);

gsl_matrix_int* gsl_matrix_int_alloc_from_vector_sizes(VALUE ary,
                VALUE nn1, VALUE nn2);
gsl_matrix_int* gsl_matrix_int_alloc_from_vectors(int argc, VALUE *argv);

VALUE rb_gsl_matrix_do_something(VALUE obj, void (*f)(gsl_matrix *));
VALUE rb_gsl_matrix_int_do_something(VALUE obj, void (*f)(gsl_matrix_int *));

VALUE rb_gsl_matrix_power(VALUE a, VALUE b);
VALUE rb_gsl_matrix_int_power(VALUE a, VALUE b);

void mygsl_vector_shift(gsl_vector *p, size_t n);
void mygsl_vector_int_shift(gsl_vector_int *p, size_t n);
void mygsl_vector_shift_scale2(gsl_vector *p, size_t n);
void mygsl_vector_int_shift_scale2(gsl_vector_int *p, size_t n);

VALUE rb_gsl_vector_to_s(VALUE obj);
VALUE rb_gsl_vector_int_to_s(VALUE obj);

VALUE rb_gsl_vector_add_constant(VALUE obj, VALUE x);
VALUE rb_gsl_vector_int_add_constant(VALUE obj, VALUE x);
VALUE rb_gsl_vector_scale(VALUE obj, VALUE x);
VALUE rb_gsl_vector_scale_bang(VALUE obj, VALUE x);
VALUE rb_gsl_vector_int_scale(VALUE obj, VALUE x);
VALUE rb_gsl_vector_int_scale_bang(VALUE obj, VALUE x);
gsl_vector_int* make_cvector_int_from_rarray(VALUE ary);
void cvector_int_set_from_rarray(gsl_vector_int *v, VALUE ary);
VALUE rb_gsl_range2vector(VALUE obj);
VALUE rb_gsl_range2vector_int(VALUE obj);

void Init_gsl_block_init(VALUE module);
void Init_gsl_block_int_init(VALUE module);
void Init_gsl_block_uchar_init(VALUE module);

void Init_gsl_matrix_nmf(void);

#endif
