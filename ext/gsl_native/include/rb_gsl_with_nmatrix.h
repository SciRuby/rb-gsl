#ifndef RB_GSL_WITH_NMATRIX_H
#define RB_GSL_WITH_NMATRIX_H

#include <gsl/gsl_vector.h>

#ifdef HAVE_NMATRIX_H

#include "nmatrix.h"
#include "include/rb_gsl_array.h"
extern VALUE cNMatrix;

// functions to convert GSL::Vectors to 1D NMatrix
static VALUE rb_gsl_vector_to_nmatrix(VALUE obj);
static VALUE rb_gsl_vector_int_to_nmatrix(VALUE obj);
static VALUE rb_gsl_vector_complex_to_nmatrix(VALUE obj);

// functions to convert GSL::Matrix to 2D NMatrix
static VALUE rb_gsl_matrix_to_nmatrix(VALUE obj);
static VALUE rb_gsl_matrix_int_to_nmatrix(VALUE obj);
static VALUE rb_gsl_matrix_complex_to_nmatrix(VALUE obj);

// functions to convert NMatrix to GSL::Vector
gsl_vector* rb_gsl_nmatrix_to_gv(VALUE nm);
gsl_vector_int* rb_gsl_nmatrix_to_gv_int(VALUE nm);
gsl_vector_complex* rb_gsl_nmatrix_to_gv_complex(VALUE nm);

// functions to convert NMatrix to GSL::Matrix
gsl_matrix* rb_gsl_nmatrix_to_gm(VALUE nm);
gsl_matrix_int* rb_gsl_nmatrix_to_gm_int(VALUE nm);
gsl_matrix_complex* rb_gsl_nmatrix_to_gm_complex(VALUE nm);

// singleton functions to convert NMatrix to GSL::Vector
VALUE rb_gsl_nmatrix_to_gsl_vector(VALUE obj, VALUE nmatrix);
VALUE rb_gsl_nmatrix_to_gsl_vector_int(VALUE obj, VALUE nmatrix);
VALUE rb_gsl_nmatrix_to_gsl_vector_complex(VALUE obj, VALUE nmatrix);

// singleton functions to convert NMatrix to GSL::Matrix
VALUE rb_gsl_nmatrix_to_gsl_matrix(VALUE obj, VALUE nmatrix);
VALUE rb_gsl_nm_to_gsl_matrix_int(VALUE obj, VALUE nmatrix);
VALUE rb_gsl_nmatrix_to_gsl_matrix_complex(VALUE obj, VALUE nmatrix);

#endif // HAVE_NMATRIX_H
#endif // RB_GSL_WITH_NMATRIX_H
