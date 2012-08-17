//#ifdef HAVE_NMATRIX_H

#include "nmatrix.h"
#include "gsl/gsl_vector.h"

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

//#endif // HAVE_NMATRIX_H
