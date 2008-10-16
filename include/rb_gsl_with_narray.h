#ifdef HAVE_NARRAY_H

#include "narray.h"
#include "gsl/gsl_vector.h"

gsl_vector* make_cvector_from_narray(VALUE);
void cvector_set_from_narray(gsl_vector*, VALUE);
void carray_set_from_narray(double*, VALUE);

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

#endif
