/*
  narray.c

  Written by Yoshiki Tsunesada       Dec/2001
  Modified by Seiya Nishizawa        14/Apr/2004
 */

#ifdef HAVE_NARRAY_H
#include "include/rb_gsl_array.h"
#include "narray.h"
#include "include/rb_gsl_with_narray.h"

static VALUE rb_gsl_na_to_gsl_matrix_method(VALUE nna);
static VALUE rb_gsl_na_to_gsl_matrix_int_method(VALUE nna);
static VALUE rb_gsl_na_to_gsl_vector_int_method(VALUE na);
static VALUE rb_gsl_vector_int_to_na(VALUE obj);
static VALUE rb_gsl_na_to_gsl_vector_int(VALUE obj, VALUE na);

/* GSL::Vector -> NArray */

static VALUE rb_gsl_vector_to_narray(VALUE obj, VALUE klass)
{
  gsl_vector *v = NULL;
  VALUE nary;
  int shape[1];
  Data_Get_Struct(obj, gsl_vector, v);
  shape[0] = v->size;
  nary = na_make_object(NA_DFLOAT, 1, shape, klass);
  if (v->stride == 1) {
    memcpy(NA_PTR_TYPE(nary,double*), v->data, shape[0]*sizeof(double));
  } else {
    int i;
    for(i=0; i < (int) v->size; i++) {
      (NA_PTR_TYPE(nary,double*))[i] = gsl_vector_get(v, i);
    }
  }
  return nary;
}

static VALUE rb_gsl_vector_complex_to_narray(VALUE obj, VALUE klass)
{
  gsl_vector_complex *v = NULL;
  VALUE nary;
  int shape[1];
  Data_Get_Struct(obj, gsl_vector_complex, v);
  shape[0] = v->size;
  nary = na_make_object(NA_DCOMPLEX, 1, shape, klass);
  if (v->stride == 1) {
    memcpy(NA_PTR_TYPE(nary,double*), v->data, shape[0]*2*sizeof(double));
  } else {
    int i;
    for(i=0; i < (int) (2*v->size); i++) {
      (NA_PTR_TYPE(nary,gsl_complex*))[i] = gsl_vector_complex_get(v, i);
    }
  }
  return nary;
}

static VALUE rb_gsl_vector_to_na(VALUE obj)
{
  VALUE na = Qnil;

  if(VECTOR_P(obj))
    na = rb_gsl_vector_to_narray(obj, cNArray);
  else if(VECTOR_COMPLEX_P(obj))
    na = rb_gsl_vector_complex_to_narray(obj, cNArray);
  else
    rb_raise(rb_eRuntimeError, "unexpected type '%s'",
        rb_obj_classname(obj));

  return na;
}

static VALUE rb_gsl_vector_complex_to_na(VALUE obj)
{
  return rb_gsl_vector_complex_to_narray(obj, cNArray);
}

static VALUE rb_gsl_vector_to_nvector(VALUE obj)
{
  return rb_gsl_vector_to_narray(obj, cNVector);
}

static VALUE rb_gsl_vector_complex_to_nvector(VALUE obj)
{
  return rb_gsl_vector_complex_to_narray(obj, cNVector);
}

/* GSL::Vector -> NArray view */

static struct NARRAY* rb_gsl_na_view_alloc(int rank, int total, int type)
{
  struct NARRAY *na;
  na = (struct NARRAY*) malloc(sizeof(struct NARRAY));
  na->rank = rank;
  na->total = total;
  na->type = type;
  // TODO Set na->ref to a one element NArray of type NAObject that contains
  // the GSL::Vector being referenced.
  na->ref = Qtrue;       /* to initialize */
  na->shape = (int *) malloc(sizeof(int)*rank);
  return na;
}

static void rb_gsl_na_view_free(struct NARRAY *na)
{
  free((int *) na->shape);
  free((struct NARRAY *) na);
}

static VALUE rb_gsl_vector_to_narray_ref(VALUE obj, VALUE klass)
{
  gsl_vector *v = NULL;
  gsl_vector_complex *vc = NULL;
  gsl_vector_int *vi = NULL;
  VALUE nary;
  struct NARRAY *na;
  if (VECTOR_P(obj)) {
    Data_Get_Struct(obj, gsl_vector, v);
    if (v->stride != 1) {
      rb_raise(rb_eRuntimeError, "Cannot make a reference obj: stride!=1");
    }
    na = rb_gsl_na_view_alloc(1, v->size, NA_DFLOAT);
    na->shape[0] = v->size;
    na->ptr = (char *) v->data;
  } else if (VECTOR_INT_P(obj)) {
    Data_Get_Struct(obj, gsl_vector_int, vi);
    if (vi->stride != 1) {
      rb_raise(rb_eRuntimeError, "Cannot make a reference obj: stride!=1");
    }
    na = rb_gsl_na_view_alloc(1, vi->size, NA_LINT);
    na->shape[0] = vi->size;
    na->ptr = (char *) vi->data;
  } else if (VECTOR_COMPLEX_P(obj)) {
    Data_Get_Struct(obj, gsl_vector_complex, vc);
    if (vc->stride != 1) {
      rb_raise(rb_eRuntimeError, "Cannot make a reference obj: stride!=1");
    }
    na = rb_gsl_na_view_alloc(1, vc->size, NA_DCOMPLEX);
    na->shape[0] = vc->size;
    na->ptr = (char *) vc->data;
  } else {
    rb_raise(rb_eRuntimeError, "cannot convert %s to NArray reference object",
        rb_obj_classname(obj));
  }
  nary = Data_Wrap_Struct(klass, 0, rb_gsl_na_view_free, na);
  return nary;
}

static VALUE rb_gsl_vector_to_na_ref(VALUE obj)
{
  return rb_gsl_vector_to_narray_ref(obj, cNArray);
}

static VALUE rb_gsl_vector_to_nvector_ref(VALUE obj)
{
  return rb_gsl_vector_to_narray_ref(obj, cNVector);
}

static VALUE rb_gsl_vector_int_to_narray(VALUE obj, VALUE klass)
{
  gsl_vector_int *v = NULL;
  VALUE nary;
  int shape[1];
  Data_Get_Struct(obj, gsl_vector_int, v);
  shape[0] = v->size;
  nary = na_make_object(NA_LINT, 1, shape, klass);
  if (v->stride == 1) {
    memcpy(NA_PTR_TYPE(nary,int*), v->data, shape[0]*sizeof(int));
  } else {
    int i;
    for(i=0; i < (int) v->size; i++) {
      (NA_PTR_TYPE(nary,int*))[i] = gsl_vector_int_get(v, i);
    }
  }
  return nary;
}

static VALUE rb_gsl_vector_int_to_na(VALUE obj)
{
  return rb_gsl_vector_int_to_narray(obj, cNArray);
}

static VALUE rb_gsl_vector_int_to_nvector(VALUE obj)
{
  return rb_gsl_vector_int_to_narray(obj, cNVector);
}

/* singleton method */
static VALUE rb_gsl_na_to_gsl_vector(VALUE obj, VALUE na)
{
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free,
                          na_to_gv(na));
}

static VALUE rb_gsl_na_to_gsl_vector_view(VALUE obj, VALUE na)
{
  return Data_Wrap_Struct(cgsl_vector_view, 0, gsl_vector_view_free,
                          na_to_gv_view(na));
}

static VALUE rb_gsl_na_to_gsl_vector_complex(VALUE obj, VALUE na)
{
  return Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free,
                          na_to_gv_complex(na));
}

static VALUE rb_gsl_na_to_gsl_vector_complex_view(VALUE obj, VALUE na)
{
  return Data_Wrap_Struct(cgsl_vector_complex_view, 0, gsl_vector_complex_view_free,
                          na_to_gv_complex_view(na));
}

static VALUE rb_gsl_na_to_gsl_vector_int(VALUE obj, VALUE na)
{
  return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free,
                          na_to_gv_int(na));
}

static VALUE rb_gsl_na_to_gsl_vector_int_view(VALUE obj, VALUE na)
{
  return Data_Wrap_Struct(cgsl_vector_int_view, 0, rb_gsl_vector_int_view_free,
                          na_to_gv_int_view(na));
}

static VALUE rb_gsl_na_to_gsl_vector_method(VALUE na)
{
  VALUE v;

  if(NA_TYPE(na) == NA_SCOMPLEX || NA_TYPE(na) == NA_DCOMPLEX)
    v = Data_Wrap_Struct(cgsl_vector_complex, 0, gsl_vector_complex_free,
        na_to_gv_complex(na));
  else
    v = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free,
        na_to_gv(na));

  return v;
}

VALUE rb_gsl_na_to_gsl_vector_view_method(VALUE na)
{
  VALUE v;

  if(NA_TYPE(na) == NA_SCOMPLEX || NA_TYPE(na) == NA_DCOMPLEX)
    v = Data_Wrap_Struct(cgsl_vector_complex_view, 0, gsl_vector_complex_view_free,
        na_to_gv_complex_view(na));
  else
    v = Data_Wrap_Struct(cgsl_vector_view, 0, gsl_vector_view_free,
        na_to_gv_view(na));

  return v;
}

static VALUE rb_gsl_na_to_gsl_vector_int_method(VALUE na)
{
  return Data_Wrap_Struct(cgsl_vector_int, 0, gsl_vector_int_free,
                          na_to_gv_int(na));
}

static VALUE rb_gsl_na_to_gsl_vector_int_view_method(VALUE na)
{
  return Data_Wrap_Struct(cgsl_vector_int_view, 0, rb_gsl_vector_int_view_free,
                          na_to_gv_int_view(na));
}

gsl_vector* na_to_gv(VALUE na)
{
  gsl_vector *v = NULL;
  VALUE nary = na;
  v = gsl_vector_alloc(NA_TOTAL(na));
  if(NA_TYPE(na) != NA_DFLOAT) {
    nary = na_change_type(na, NA_DFLOAT);
  }
  memcpy(v->data, NA_PTR_TYPE(nary,double*), v->size*sizeof(double));
  return v;
}

gsl_vector_view* na_to_gv_view(VALUE na)
{
  gsl_vector_view *v = NULL;

  // Raise exception if na's type is not NA_DFLOAT.
  if(NA_TYPE(na) != NA_DFLOAT)
    rb_raise(rb_eTypeError, "GSL::Vector::View requires NArray be DFLOAT");

  v = gsl_vector_view_alloc();
  v->vector.data = NA_PTR_TYPE(na,double*);
  v->vector.size = NA_TOTAL(na);
  v->vector.stride = 1;
  v->vector.owner = 0;
  return v;
}

gsl_vector_complex* na_to_gv_complex(VALUE na)
{
  gsl_vector_complex *v = NULL;
  VALUE nary = na;
  v = gsl_vector_complex_alloc(NA_TOTAL(na));
  if(NA_TYPE(na) != NA_DCOMPLEX) {
    nary = na_change_type(na, NA_DCOMPLEX);
  }
  memcpy(v->data, NA_PTR_TYPE(nary,gsl_complex*), v->size*sizeof(gsl_complex));
  return v;
}

gsl_vector_complex_view* na_to_gv_complex_view(VALUE na)
{
  gsl_vector_complex_view *v = NULL;

  // Raise exception if na's type is not NA_DCOMPLEX
  if(NA_TYPE(na) != NA_DCOMPLEX)
    rb_raise(rb_eTypeError, "GSL::Vector::Complex::View requires NArray be DCOMPLEX");

  v = gsl_vector_complex_view_alloc();
  v->vector.data = NA_PTR_TYPE(na,double*);
  v->vector.size = NA_TOTAL(na);
  v->vector.stride = 1;
  v->vector.owner = 0;
  return v;
}

gsl_vector_int* na_to_gv_int(VALUE na)
{
  gsl_vector_int *v = NULL;
  VALUE nary = na;
  v = gsl_vector_int_alloc(NA_TOTAL(na));
  if(NA_TYPE(na) != NA_LINT) {
    nary = na_change_type(na, NA_LINT);
  }
  memcpy(v->data, NA_PTR_TYPE(nary,int*), v->size*sizeof(int));
  return v;
}

gsl_vector_int_view* na_to_gv_int_view(VALUE na)
{
  gsl_vector_int_view *v = NULL;

  // Raise exception if na's type is not NA_LINT
  if(NA_TYPE(na) != NA_LINT)
    rb_raise(rb_eTypeError, "GSL::Vector::Int::View requires NArray be LINT");
  v = rb_gsl_vector_int_view_alloc(NA_TOTAL(na));
  v->vector.data = NA_PTR_TYPE(na,int*);
  v->vector.size = NA_TOTAL(na);
  v->vector.stride = 1;
  v->vector.owner = 0;
  return v;
}

static VALUE rb_gsl_matrix_to_narray(VALUE obj, VALUE klass)
{
  gsl_matrix *m = NULL;
  VALUE nary;
  int shape[2];
  size_t i;
  Data_Get_Struct(obj, gsl_matrix, m);
  shape[0] = m->size2;
  shape[1] = m->size1;
  nary = na_make_object(NA_DFLOAT, 2, shape, klass);
  for (i = 0; (int) i < shape[1]; i++) {
    memcpy(NA_PTR_TYPE(nary,double*)+(i*shape[0]), m->data+(i*m->tda),
           shape[0]*sizeof(double));
  }
  return nary;
}

static VALUE rb_gsl_matrix_to_na(VALUE obj, VALUE klass)
{
  return rb_gsl_matrix_to_narray(obj, cNArray);
}

static VALUE rb_gsl_matrix_to_nmatrix(VALUE obj, VALUE klass)
{
  return rb_gsl_matrix_to_narray(obj, cNMatrix);
}

static VALUE rb_gsl_matrix_int_to_narray(VALUE obj, VALUE klass)
{
  gsl_matrix_int *m = NULL;
  VALUE nary;
  int shape[2];
  size_t i;
  Data_Get_Struct(obj, gsl_matrix_int, m);
  shape[0] = m->size2;
  shape[1] = m->size1;
  nary = na_make_object(NA_LINT, 2, shape, klass);
  for (i = 0; (int) i < shape[1]; i++) {
    memcpy(NA_PTR_TYPE(nary,int*)+(i*shape[0]), m->data+(i*m->tda),
           shape[0]*sizeof(int));
  }
  return nary;
}

static VALUE rb_gsl_matrix_int_to_na(VALUE obj)
{
  return rb_gsl_matrix_int_to_narray(obj, cNArray);
}

static VALUE rb_gsl_matrix_int_to_nmatrix(VALUE obj)
{
  return rb_gsl_matrix_int_to_narray(obj, cNMatrix);
}

static VALUE rb_gsl_matrix_to_narray_ref(VALUE obj, VALUE klass)
{
  gsl_matrix *m = NULL;
  VALUE nary;
  struct NARRAY *na;
  Data_Get_Struct(obj, gsl_matrix, m);
  if (m->tda != m->size2) {
    rb_raise(rb_eRuntimeError, "Cannot make a reference obj: non-contiguous");
  }
  na = rb_gsl_na_view_alloc(2, m->size2 * m->size1, NA_DFLOAT);
  na->shape[0] = m->size2;
  na->shape[1] = m->size1;
  na->ptr = (char *) m->data;
  nary = Data_Wrap_Struct(klass, 0, rb_gsl_na_view_free, na);
  return nary;
}

static VALUE rb_gsl_matrix_to_na_ref(VALUE obj)
{
  return rb_gsl_matrix_to_narray_ref(obj, cNArray);
}

static VALUE rb_gsl_matrix_to_nmatrix_ref(VALUE obj)
{
  return rb_gsl_matrix_to_narray_ref(obj, cNMatrix);
}

static VALUE rb_gsl_matrix_int_to_narray_ref(VALUE obj, VALUE klass)
{
  gsl_matrix_int *m = NULL;
  VALUE nary;
  struct NARRAY *na;
  Data_Get_Struct(obj, gsl_matrix_int, m);
  if (m->tda != m->size2) {
    rb_raise(rb_eRuntimeError, "Cannot make a reference obj: non-contiguous");
  }
  na = rb_gsl_na_view_alloc(2, m->size2 * m->size1, NA_LINT);
  na->shape[0] = m->size2;
  na->shape[1] = m->size1;
  na->ptr = (char *) m->data;
  nary = Data_Wrap_Struct(klass, 0, rb_gsl_na_view_free, na);
  return nary;
}

static VALUE rb_gsl_matrix_int_to_na_ref(VALUE obj)
{
  return rb_gsl_matrix_int_to_narray_ref(obj, cNArray);
}

static VALUE rb_gsl_matrix_int_to_nmatrix_ref(VALUE obj)
{
  return rb_gsl_matrix_int_to_narray_ref(obj, cNMatrix);
}

/* NArray -> GSL::Matrix */
VALUE rb_gsl_na_to_gsl_matrix(VALUE obj, VALUE nna)
{
  gsl_matrix *m = NULL;
  m = na_to_gm(nna);
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, m);
}

VALUE rb_gsl_na_to_gsl_matrix_view(VALUE obj, VALUE nna)
{
  gsl_matrix_view *m = NULL;
  m = na_to_gm_view(nna);
  return Data_Wrap_Struct(cgsl_matrix_view, 0, gsl_matrix_view_free, m);
}

VALUE rb_gsl_na_to_gsl_matrix_int(VALUE obj, VALUE nna)
{
  gsl_matrix_int *m = NULL;
  m = na_to_gm_int(nna);
  return Data_Wrap_Struct(cgsl_matrix_int, 0, gsl_matrix_int_free, m);
}

VALUE rb_gsl_na_to_gsl_matrix_int_view(VALUE obj, VALUE nna)
{
  gsl_matrix_int_view *m = NULL;
  m = na_to_gm_int_view(nna);
  return Data_Wrap_Struct(cgsl_matrix_int_view, 0, rb_gsl_matrix_int_view_free, m);
}

static VALUE rb_gsl_na_to_gsl_matrix_method(VALUE nna)
{
  gsl_matrix *m = NULL;
  m = na_to_gm(nna);
  return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, m);
}

static VALUE rb_gsl_na_to_gsl_matrix_view_method(VALUE nna)
{
  gsl_matrix_view *m = NULL;
  m = na_to_gm_view(nna);
  return Data_Wrap_Struct(cgsl_matrix_view, 0, gsl_matrix_view_free, m);
}

static VALUE rb_gsl_na_to_gsl_matrix_int_method(VALUE nna)
{
  gsl_matrix_int *m = NULL;
  m = na_to_gm_int(nna);
  return Data_Wrap_Struct(cgsl_matrix_int, 0, gsl_matrix_int_free, m);
}

static VALUE rb_gsl_na_to_gsl_matrix_int_view_method(VALUE nna)
{
  gsl_matrix_int_view *m = NULL;
  m = na_to_gm_int_view(nna);
  return Data_Wrap_Struct(cgsl_matrix_int_view, 0, rb_gsl_matrix_int_view_free, m);
}

gsl_matrix* na_to_gm(VALUE nna)
{
  gsl_matrix *m = NULL;
  VALUE ary2;
  struct NARRAY *na = NULL;
  GetNArray(nna, na);
  m = gsl_matrix_alloc(na->shape[1], na->shape[0]);
  ary2 = na_change_type(nna, NA_DFLOAT);
  memcpy(m->data, NA_PTR_TYPE(ary2,double*), m->size1*m->size2*sizeof(double));
  return m;
}

gsl_matrix_view* na_to_gm_view(VALUE nna)
{
  gsl_matrix_view *m = NULL;
  VALUE ary2;
  struct NARRAY *na = NULL;

  // Raise exception if nna's type is not NA_DFLOAT
  if(NA_TYPE(nna) != NA_DFLOAT)
    rb_raise(rb_eTypeError, "GSL::Matrix::View requires NArray be DFLOAT");
  GetNArray(nna, na);
  m = gsl_matrix_view_alloc();
  ary2 = na_change_type(nna, NA_DFLOAT);
  m->matrix.data = NA_PTR_TYPE(ary2,double*);
  m->matrix.size1 = na->shape[1];
  m->matrix.size2 = na->shape[0];
  m->matrix.tda = m->matrix.size2;
  m->matrix.owner = 0;
  return m;
}

gsl_matrix_int* na_to_gm_int(VALUE nna)
{
  gsl_matrix_int *m = NULL;
  VALUE ary2;
  struct NARRAY *na = NULL;
  GetNArray(nna, na);
  m = gsl_matrix_int_alloc(na->shape[1], na->shape[0]);
  ary2 = na_change_type(nna, NA_LINT);
  memcpy(m->data, NA_PTR_TYPE(ary2,int*), m->size1*m->size2*sizeof(int));
  return m;
}

gsl_matrix_int_view* na_to_gm_int_view(VALUE nna)
{
  gsl_matrix_int_view *m = NULL;
  VALUE ary2;
  struct NARRAY *na = NULL;

  // Raise exception if nna's type is not NA_LINT
  if(NA_TYPE(nna) != NA_LINT)
    rb_raise(rb_eTypeError, "GSL::Matrix::Int::View requires NArray be LINT");
  GetNArray(nna, na);
  m = rb_gsl_matrix_int_view_alloc(na->shape[1], na->shape[0]);
  ary2 = na_change_type(nna, NA_LINT);
  m->matrix.data = NA_PTR_TYPE(ary2,int*);
  m->matrix.size1 = na->shape[1];
  m->matrix.size2 = na->shape[0];
  m->matrix.tda = m->matrix.size2;
  m->matrix.owner = 0;
  return m;
}

#ifdef HAVE_NARRAY_H
#include "narray.h"
#include <gsl/gsl_histogram.h>
EXTERN VALUE cgsl_histogram;
static VALUE rb_gsl_narray_histogram(int argc, VALUE *argv, VALUE obj)
{
  double *ptr, *ptr_range;
  gsl_histogram *h = NULL;
  gsl_vector *ranges;
  gsl_vector_view v;
  double min, max;
  size_t i, n, size, stride;
  ptr = get_vector_ptr(obj, &stride, &size);
  v.vector.data = ptr;
  v.vector.size = size;
  v.vector.size = stride;
  switch (argc) {
  case 1:
    if (rb_obj_is_kind_of(argv[0], rb_cRange))
      argv[0] = rb_gsl_range2ary(argv[0]);
    switch (TYPE(argv[0])) {
    case T_FIXNUM:
      n = NUM2INT(argv[0]);
      min = gsl_vector_min(&v.vector) - 4*GSL_DBL_EPSILON;
      max = gsl_vector_max(&v.vector) + 4*GSL_DBL_EPSILON;
      h = gsl_histogram_alloc(n);
      gsl_histogram_set_ranges_uniform(h, min, max);
      break;
    case T_ARRAY:
      //      n = RARRAY(argv[0])->len - 1;
      n = RARRAY_LEN(argv[0]) - 1;
      h = gsl_histogram_alloc(n);
      for (i = 0; i <= n; i++) h->range[i] = NUM2DBL(rb_ary_entry(argv[0], i));
      break;
    default:
      if (VECTOR_P(argv[0])) {
        Data_Get_Struct(argv[0], gsl_vector, ranges);
        n = ranges->size - 1;
        h = gsl_histogram_alloc(n);
        gsl_histogram_set_ranges(h, ranges->data, ranges->size);
      } else if (NA_IsNArray(argv[0])) {
        ptr_range = get_vector_ptr(argv[0], &stride, &n);
        h = gsl_histogram_alloc(n);
        gsl_histogram_set_ranges(h, ptr_range, n);
      }
      break;
    }
    break;
  case 2:
    n = NUM2INT(argv[0]);
    switch (TYPE(argv[1])) {
    case T_ARRAY:
      min = NUM2DBL(rb_ary_entry(argv[1], 0));
      max = NUM2DBL(rb_ary_entry(argv[1], 1));
      break;
    default:
      rb_raise(rb_eTypeError, "wrong argument type %s (Array expected)",
               rb_class2name(CLASS_OF(argv[1])));
      break;
    }
    h = gsl_histogram_alloc(n);
    gsl_histogram_set_ranges_uniform(h, min, max);
    break;
  case 3:
    n = NUM2INT(argv[0]);
    min = NUM2DBL(argv[1]); max = NUM2DBL(argv[2]);
    h = gsl_histogram_alloc(n);
    gsl_histogram_set_ranges_uniform(h, min, max);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments %d", argc);
    break;
  }
  for (i = 0; i < size; i++)
    gsl_histogram_increment(h, ptr[i*stride]);
  return Data_Wrap_Struct(cgsl_histogram, 0, gsl_histogram_free, h);
}
#endif

/*void rb_gsl_with_narray_define_methods()*/
void Init_gsl_narray(VALUE module)
{
  rb_define_method(cgsl_vector, "to_na", rb_gsl_vector_to_na, 0);
  rb_define_alias(cgsl_vector, "to_narray", "to_na");
  rb_define_method(cgsl_vector, "to_na2", rb_gsl_vector_to_na_ref, 0);
  rb_define_alias(cgsl_vector, "to_narray_ref", "to_na2");
  rb_define_alias(cgsl_vector, "to_na_ref", "to_na2");

  rb_define_method(cgsl_vector, "to_nv", rb_gsl_vector_to_nvector, 0);
  rb_define_alias(cgsl_vector, "to_nvector", "to_nv");
  rb_define_method(cgsl_vector, "to_nv_ref", rb_gsl_vector_to_nvector_ref, 0);
  rb_define_alias(cgsl_vector, "to_nvector_ref", "to_nv_ref");

  rb_define_singleton_method(cgsl_vector, "to_gslv", rb_gsl_na_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector, "to_gv", rb_gsl_na_to_gsl_vector, 1);

  rb_define_singleton_method(cgsl_vector, "na_to_gslv", rb_gsl_na_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector, "na_to_gv", rb_gsl_na_to_gsl_vector, 1);
  rb_define_singleton_method(cgsl_vector, "to_gslv_view", rb_gsl_na_to_gsl_vector_view, 1);

  rb_define_singleton_method(cgsl_vector, "to_gv_view", rb_gsl_na_to_gsl_vector_view, 1);

  rb_define_singleton_method(cgsl_vector, "na_to_gslv_view", rb_gsl_na_to_gsl_vector_view, 1);
  rb_define_singleton_method(cgsl_vector, "na_to_gv_view", rb_gsl_na_to_gsl_vector_view, 1);

  rb_define_method(cNArray, "to_gv", rb_gsl_na_to_gsl_vector_method, 0);
  rb_define_alias(cNArray, "to_gslv", "to_gv");
  rb_define_method(cNArray, "to_gv_view", rb_gsl_na_to_gsl_vector_view_method, 0);
  rb_define_alias(cNArray, "to_gslv_view", "to_gv_view");
  rb_define_alias(cNArray, "to_gv2", "to_gv_view");
  rb_define_alias(cNArray, "to_gv_ref", "to_gv_view");
  /*****/
  rb_define_method(cgsl_vector_complex, "to_na", rb_gsl_vector_complex_to_na, 0);
  rb_define_alias(cgsl_vector_complex, "to_narray", "to_na");

  rb_define_method(cgsl_vector_complex, "to_na2", rb_gsl_vector_to_na_ref, 0);
  rb_define_alias(cgsl_vector_complex, "to_narray_ref", "to_na2");
  rb_define_alias(cgsl_vector_complex, "to_na_ref", "to_na2");

  rb_define_method(cgsl_vector_complex, "to_nv", rb_gsl_vector_complex_to_nvector, 0);
  rb_define_alias(cgsl_vector_complex, "to_nvector", "to_nv");

  rb_define_method(cgsl_vector_complex, "to_nv2", rb_gsl_vector_to_nvector_ref, 0);
  rb_define_alias(cgsl_vector_complex, "to_nv_ref", "to_nv2");
  rb_define_alias(cgsl_vector_complex, "to_nvector_ref", "to_nv2");

  rb_define_singleton_method(cgsl_vector_complex, "to_gslv", rb_gsl_na_to_gsl_vector_complex, 1);
  rb_define_singleton_method(cgsl_vector_complex, "to_gv", rb_gsl_na_to_gsl_vector_complex, 1);
  rb_define_singleton_method(cgsl_vector_complex, "na_to_gslv", rb_gsl_na_to_gsl_vector_complex, 1);
  rb_define_singleton_method(cgsl_vector_complex, "na_to_gv", rb_gsl_na_to_gsl_vector_complex, 1);

  rb_define_singleton_method(cgsl_vector_complex, "to_gslv_view", rb_gsl_na_to_gsl_vector_complex_view, 1);
  rb_define_singleton_method(cgsl_vector_complex, "to_gv_view", rb_gsl_na_to_gsl_vector_complex_view, 1);
  rb_define_singleton_method(cgsl_vector_complex, "na_to_gslv_view", rb_gsl_na_to_gsl_vector_complex_view, 1);
  rb_define_singleton_method(cgsl_vector_complex, "na_to_gv_view", rb_gsl_na_to_gsl_vector_complex_view, 1);

  /*****/
  rb_define_method(cgsl_vector_int, "to_na", rb_gsl_vector_int_to_na, 0);
  rb_define_alias(cgsl_vector_int, "to_narray", "to_na");
  rb_define_method(cgsl_vector_int, "to_nv", rb_gsl_vector_int_to_nvector, 0);
  rb_define_alias(cgsl_vector_int, "to_nvector", "to_nv");
  rb_define_method(cgsl_vector_int, "to_na2", rb_gsl_vector_to_na_ref, 0);
  rb_define_alias(cgsl_vector_int, "to_na_ref", "to_na2");
  rb_define_alias(cgsl_vector_int, "to_narray_ref", "to_na2");
  rb_define_method(cgsl_vector_int, "to_nv_ref", rb_gsl_vector_to_nvector_ref, 0);
  rb_define_alias(cgsl_vector_int, "to_nvector_ref", "to_nv_ref");

  rb_define_singleton_method(cgsl_vector_int, "to_gslv", rb_gsl_na_to_gsl_vector_int, 1);
  rb_define_singleton_method(cgsl_vector_int, "to_gv", rb_gsl_na_to_gsl_vector_int, 1);

  rb_define_singleton_method(cgsl_vector_int, "na_to_gslv", rb_gsl_na_to_gsl_vector_int, 1);
  rb_define_singleton_method(cgsl_vector_int, "na_to_gv", rb_gsl_na_to_gsl_vector_int, 1);

  rb_define_singleton_method(cgsl_vector_int, "to_gslv_view", rb_gsl_na_to_gsl_vector_int_view, 1);
  rb_define_singleton_method(cgsl_vector_int, "to_gv_view", rb_gsl_na_to_gsl_vector_int_view, 1);

  rb_define_singleton_method(cgsl_vector_int, "na_to_gslv_view", rb_gsl_na_to_gsl_vector_int_view, 1);
  rb_define_singleton_method(cgsl_vector_int, "na_to_gv_view", rb_gsl_na_to_gsl_vector_int_view, 1);

  rb_define_method(cNArray, "to_gv_int", rb_gsl_na_to_gsl_vector_int_method, 0);
  rb_define_alias(cNArray, "to_gslv_int", "to_gv_int");
  rb_define_method(cNArray, "to_gv_int_view", rb_gsl_na_to_gsl_vector_int_view_method, 0);
  rb_define_alias(cNArray, "to_gslv_int_view", "to_gv_int_view");
  /*****/

  rb_define_method(cgsl_matrix, "to_na", rb_gsl_matrix_to_na, 0);
  rb_define_alias(cgsl_matrix, "to_narray", "to_na");
  rb_define_method(cgsl_matrix, "to_nm", rb_gsl_matrix_to_nmatrix, 0);
  rb_define_alias(cgsl_matrix, "to_nmatrix", "to_nm");
  rb_define_method(cgsl_matrix, "to_na2", rb_gsl_matrix_to_na_ref, 0);
  rb_define_alias(cgsl_matrix, "to_na_ref", "to_na2");
  rb_define_alias(cgsl_matrix, "to_narray_ref", "to_na2");
  rb_define_method(cgsl_matrix, "to_nm_ref", rb_gsl_matrix_to_nmatrix_ref, 0);
  rb_define_alias(cgsl_matrix, "to_nmatrix_ref", "to_nm_ref");

  rb_define_singleton_method(cgsl_matrix, "to_gslm", rb_gsl_na_to_gsl_matrix, 1);
  rb_define_singleton_method(cgsl_matrix, "to_gm", rb_gsl_na_to_gsl_matrix, 1);

  rb_define_singleton_method(cgsl_matrix, "na_to_gslm", rb_gsl_na_to_gsl_matrix, 1);
  rb_define_singleton_method(cgsl_matrix, "na_to_gm", rb_gsl_na_to_gsl_matrix, 1);

  rb_define_singleton_method(cgsl_matrix, "to_gslm_view", rb_gsl_na_to_gsl_matrix_view, 1);
  rb_define_singleton_method(cgsl_matrix, "to_gm_view", rb_gsl_na_to_gsl_matrix_view, 1);

  rb_define_singleton_method(cgsl_matrix, "na_to_gslm_view", rb_gsl_na_to_gsl_matrix_view, 1);
  rb_define_singleton_method(cgsl_matrix, "na_to_gm_view", rb_gsl_na_to_gsl_matrix_view, 1);

  rb_define_method(cNArray, "to_gslm", rb_gsl_na_to_gsl_matrix_method, 0);
  rb_define_alias(cNArray, "to_gm", "to_gslm");
  rb_define_method(cNArray, "to_gslm_view", rb_gsl_na_to_gsl_matrix_view_method, 0);
  rb_define_alias(cNArray, "to_gm_view", "to_gslm_view");

  /*****/
  // TODO Complex matrix

  /*****/
  rb_define_method(cgsl_matrix_int, "to_na", rb_gsl_matrix_int_to_na, 0);
  rb_define_alias(cgsl_matrix_int, "to_narray", "to_na");
  rb_define_method(cgsl_matrix_int, "to_nm", rb_gsl_matrix_int_to_nmatrix, 0);
  rb_define_alias(cgsl_matrix_int, "to_nmatrix", "to_nm");

  rb_define_method(cgsl_matrix_int, "to_na2", rb_gsl_matrix_int_to_na_ref, 0);
  rb_define_alias(cgsl_matrix_int, "to_na_ref", "to_na2");
  rb_define_alias(cgsl_matrix_int, "to_narray_ref", "to_na2");
  rb_define_method(cgsl_matrix_int, "to_nm_ref", rb_gsl_matrix_int_to_nmatrix_ref, 0);
  rb_define_alias(cgsl_matrix_int, "to_nmatrix_ref", "to_nm_ref");

  rb_define_singleton_method(cgsl_matrix_int, "to_gslm", rb_gsl_na_to_gsl_matrix_int, 1);
  rb_define_singleton_method(cgsl_matrix_int, "to_gm", rb_gsl_na_to_gsl_matrix_int, 1);
  rb_define_singleton_method(cgsl_matrix_int, "na_to_gslm", rb_gsl_na_to_gsl_matrix_int, 1);
  rb_define_singleton_method(cgsl_matrix_int, "na_to_gm", rb_gsl_na_to_gsl_matrix_int, 1);

  rb_define_singleton_method(cgsl_matrix_int, "to_gslm_view", rb_gsl_na_to_gsl_matrix_int_view, 1);
  rb_define_singleton_method(cgsl_matrix_int, "to_gm_view", rb_gsl_na_to_gsl_matrix_int_view, 1);
  rb_define_singleton_method(cgsl_matrix_int, "na_to_gslm_view", rb_gsl_na_to_gsl_matrix_int_view, 1);
  rb_define_singleton_method(cgsl_matrix_int, "na_to_gm_view", rb_gsl_na_to_gsl_matrix_int_view, 1);

  rb_define_method(cNArray, "to_gm_int", rb_gsl_na_to_gsl_matrix_int_method, 0);
  rb_define_alias(cNArray, "to_gslm_int", "to_gm_int");
  rb_define_method(cNArray, "to_gm_int_view", rb_gsl_na_to_gsl_matrix_int_view_method, 0);
  rb_define_alias(cNArray, "to_gslm_int_view", "to_gm_int_view");

  rb_define_method(cNArray, "histogram", rb_gsl_narray_histogram, -1);
}

#endif
