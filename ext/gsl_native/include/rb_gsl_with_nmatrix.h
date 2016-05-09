#ifndef RB_GSL_WITH_NMATRIX_H
#define RB_GSL_WITH_NMATRIX_H

#include <gsl/gsl_vector.h>

#ifdef HAVE_NMATRIX_H

#include "nmatrix_config.h"
#include "nmatrix.h"
#include "include/rb_gsl_array.h"
extern VALUE cNMatrix;

// nmatrix external API
extern VALUE rb_nmatrix_dense_create(nm_dtype_t dtype, size_t* shape, size_t rank, void* elements, size_t length);
extern VALUE rb_nvector_dense_create(nm_dtype_t dtype, void* elements, size_t length);

// functions to convert NMatrix to GSL::Vector
gsl_vector* rb_gsl_nmatrix_to_gv(VALUE nm);
gsl_vector_int* rb_gsl_nmatrix_to_gv_int(VALUE nm);
gsl_vector_complex* rb_gsl_nmatrix_to_gv_complex(VALUE nm);

// functions to convert NMatrix to GSL::Matrix
gsl_matrix* rb_gsl_nmatrix_to_gm(VALUE nm);
gsl_matrix_int* rb_gsl_nmatrix_to_gm_int(VALUE nm);
gsl_matrix_complex* rb_gsl_nmatrix_to_gm_complex(VALUE nm);


#endif // HAVE_NMATRIX_H
#endif // RB_GSL_WITH_NMATRIX_H
