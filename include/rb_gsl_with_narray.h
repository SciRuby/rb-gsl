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
#endif // RB_GSL_WITH_NARRAY_H
