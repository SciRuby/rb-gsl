#ifndef RB_GSL_WITH_NMATRIX_H
#define RB_GSL_WITH_NMATRIX_H

#include <gsl/gsl_vector.h>

#ifdef HAVE_NMATRIX_H

#include "nmatrix.h"
#include "include/rb_gsl_array.h"
extern VALUE cNMatrix;


#endif // HAVE_NMATRIX_H
#endif // RB_GSL_WITH_NMATRIX_H
