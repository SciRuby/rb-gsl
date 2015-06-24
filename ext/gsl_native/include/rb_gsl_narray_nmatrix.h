#ifndef RB_GSL_WITH_NARRAY_H
#define RB_GSL_WITH_NARRAY_H

#include <gsl/gsl_vector.h>

#ifdef HAVE_NARRAY_H

#include "narray.h"

gsl_vector* make_cvector_from_narray(VALUE);
void cvector_set_from_narray(gsl_vector*, VALUE);
void carray_set_from_narray(double*, VALUE);

VALUE rb_gsl_na_to_gsl_vector_view_method(VALUE na);
VALUE rb_gsl_na_to_gsl_matrix(VALUE obj, VALUE nna);
gsl_vector_view* na_to_gv_view(VALUE na);
gsl_matrix_view* na_to_gm_view(VALUE nna);

gsl_vector_int_view* na_to_gv_int_view(VALUE na);
gsl_matrix_int_view* na_to_gm_int_view(VALUE nna);
gsl_vector* na_to_gv(VALUE na);
gsl_vector_int* na_to_gv_int(VALUE na);
gsl_matrix* na_to_gm(VALUE nna);
gsl_matrix_int* na_to_gm_int(VALUE nna);
extern VALUE cNVector, cNMatrix;

gsl_vector_complex* na_to_gv_complex(VALUE na);
gsl_vector_complex_view* na_to_gv_complex_view(VALUE na);

#endif // HAVE_NARRAY_H

#ifdef HAVE_NMATRIX_H

#include "nmatrix.h"

// in array.c
gsl_vector* make_cvector_from_nvector(VALUE vec);
void cvector_set_from_nvector(gsl_vector*, VALUE);
void carray_set_from_nvector(double*, VALUE);

// in gsl_nmatrix.c
VALUE rb_gsl_nm_to_gsl_matrix(VALUE obj, VALUE nm);
VALUE rb_gsl_nv_to_gsl_vector(VALUE obj, VALUE nv);
gsl_matrix* nm_to_gm(VALUE nm);
gsl_vector* nv_to_gv(VALUE nv);

VALUE rb_gsl_nm_to_gsl_matrix_complex(VALUE obj, VALUE nm);
VALUE rb_gsl_nv_to_gsl_vector_complex(VALUE obj, VALUE nv);
gsl_matrix_complex* nm_to_gm_complex(VALUE nm);
gsl_vector_complex* nv_to_gv_complex(VALUE nv);

VALUE rb_gsl_nm_to_gsl_matrix_int(VALUE obj, VALUE nm);
VALUE rb_gsl_nv_to_gsl_vector_int(VALUE obj, VALUE nv);
gsl_matrix_int* nm_to_gm_int(VALUE nm);
gsl_vector_int* nv_to_gv_int(VALUE nv);

// from NMatrix:
extern VALUE cNVector, CNMatrix;

#endif // HAVE_NMATRIX_H

#endif // RB_GSL_WITH_NARRAY_H
