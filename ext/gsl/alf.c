/*
 * Ruby interface for ALF, a library to calculate associated Legendre functions by Patrick Alken
 * Based on ALF-0.1
 */
#ifdef HAVE_ALF_ALF_H
#include "include/rb_gsl.h"

static VALUE cWspace;
static VALUE rb_alf_alloc(VALUE klass, VALUE lmax)
{
  alf_workspace *w = NULL;
  w = alf_alloc(FIX2INT(lmax));
  return Data_Wrap_Struct(cWspace, 0, alf_free, w);
}

static VALUE rb_alf_params(VALUE obj, VALUE csphase, VALUE cnorm, VALUE norm)
{
  alf_workspace *w;
  int ret;
  Data_Get_Struct(obj, alf_workspace, w);
  ret = alf_params(FIX2INT(csphase), FIX2INT(cnorm), (alf_norm_t) FIX2INT(norm), w);
  return INT2FIX(ret);
}

static void define_constants(VALUE klass)
{
  rb_define_const(klass, "NORM_NONE", INT2FIX((int) ALF_NORM_NONE));
  rb_define_const(klass, "NORM_SPHARM", INT2FIX((int) ALF_NORM_SPHARM));
  rb_define_const(klass, "NORM_ORTHO", INT2FIX((int) ALF_NORM_ORTHO));
  rb_define_const(klass, "NORM_SCHMIDT", INT2FIX((int) ALF_NORM_SCHMIDT));
}

/**
 * arguments:
 *   - Plm_array(x) : lmax = w->lmax, A new vector is created
 *   - Plm_array(x, result) : lmax = w->lmax, the given vector is used
 *   - Plm_array(lmax, x) : A new vector is created
 *   - Plm_array(lmax, x, result) : Same as C Plm_array()
 *   - Plm_array(x, result, deriv) : lmax = w->lmax, calcurate Plm_deriv_array(lmax, x, result, deriv)
 *   - Plm_array(lmax, x, result, deriv) : Same as C alf_Plm_deriv_array
 */
static VALUE rb_alf_Plm_array(int argc, VALUE *argv, VALUE obj)
{
  alf_workspace *w = NULL;
  gsl_vector *res = NULL, *deriv = NULL;
  int lmax;
  double x;
  VALUE ret;
  Data_Get_Struct(obj, alf_workspace, w);
  switch (argc) {
  case 1:
    x = NUM2DBL(argv[0]);
    lmax = w->lmax;
    res = gsl_vector_alloc(alf_array_size(lmax));
    ret = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, res);
    break;
  case 2: // Plm_array(x, result) or Plm_array(lmax, x)
    if (VECTOR_P(argv[1])) {
      x = NUM2DBL(argv[0]);
      Data_Get_Struct(argv[1], gsl_vector, res);
      lmax = w->lmax;
      if (res->size < alf_array_size(lmax)) {
        rb_raise(rb_eRuntimeError, "Vector length is too small. (%d for >= %d\n", (int) res->size,
                 (int) alf_array_size(lmax));
      }
      ret = argv[1];
    } else {
      lmax = FIX2INT(argv[0]);
      x = NUM2DBL(argv[1]);
      res = gsl_vector_alloc(alf_array_size(lmax));
      ret = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, res);
    }
    break;
  case 3: // Plm_array(lmax, x, result) or Plm_array(x, result, deriv)
    if (VECTOR_P(argv[1])) {
      CHECK_VECTOR(argv[2]);
      lmax = w->lmax;
      x = NUM2DBL(argv[0]);
      Data_Get_Struct(argv[1], gsl_vector, res);
      Data_Get_Struct(argv[2], gsl_vector, deriv);
      ret = argv[1];
    } else {
      lmax = FIX2INT(argv[0]);
      x = NUM2DBL(argv[1]);
      CHECK_VECTOR(argv[2]);
      Data_Get_Struct(argv[2], gsl_vector, res);
      if (res->size < alf_array_size(lmax)) {
        rb_raise(rb_eRuntimeError, "Vector length is too small. (%d for >= %d\n", (int) res->size,
                 (int) alf_array_size(lmax));
      }
      ret = argv[2];
    }
    break;
  case 4:
    CHECK_VECTOR(argv[2]); CHECK_VECTOR(argv[3])
    lmax = FIX2INT(argv[0]);
    x = NUM2DBL(argv[1]);
    Data_Get_Struct(argv[2], gsl_vector, res);
    Data_Get_Struct(argv[3], gsl_vector, deriv);
    ret = argv[2];
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of argumentso (%d for 1-3)\n", argc);
  }
  if (argc == 4 && deriv != NULL) alf_Plm_deriv_array(lmax, x, res->data, deriv->data, w);
  else alf_Plm_array(lmax, x, res->data, w);
  return ret;
}

/**
 * arguments:
 *   - Plm_deriv_array(x) : lmax = w->lmax, two new vectors are created and returned as an array
 *   - Plm_array(lmax, x) : Two new vectors are created and returned as an array
 *   - Plm_array(lmax, x, result, deriv) : Same as C alf_Plm_deriv_array()
 */
static VALUE rb_alf_Plm_deriv_array(int argc, VALUE *argv, VALUE obj)
{
  alf_workspace *w = NULL;
  gsl_vector *res = NULL, *deriv = NULL;
  int lmax;
  double x;
  VALUE ret1, ret2, ary;
  Data_Get_Struct(obj, alf_workspace, w);
  switch (argc) {
  case 1:
    x = NUM2DBL(argv[0]);
    lmax = w->lmax;
    res = gsl_vector_alloc(alf_array_size(lmax));
    deriv = gsl_vector_alloc(alf_array_size(lmax));
    ret1 = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, res);
    ret2 = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, deriv);
    break;
  case 2:
    lmax = FIX2INT(argv[1]);
    x = NUM2DBL(argv[1]);
    res = gsl_vector_alloc(alf_array_size(lmax));
    deriv = gsl_vector_alloc(alf_array_size(lmax));
    ret1 = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, res);
    ret2 = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, deriv);
    break;
  case 3:
    CHECK_VECTOR(argv[1]);
    CHECK_VECTOR(argv[2]);
    ret1 = argv[1];
    ret2 = argv[2];
    lmax = w->lmax;
    x = NUM2DBL(argv[0]);
    break;
  case 4:
    lmax = FIX2INT(argv[0]);
    x = NUM2DBL(argv[1]);
    CHECK_VECTOR(argv[2]);
    CHECK_VECTOR(argv[3]);
    Data_Get_Struct(argv[2], gsl_vector, res);
    Data_Get_Struct(argv[3], gsl_vector, deriv);
    if (res->size < alf_array_size(lmax)) {
      rb_raise(rb_eRuntimeError, "Vector length is too small. (%d for >= %d\n", (int) res->size,
               (int) alf_array_size(lmax));
    }
    if (deriv->size < alf_array_size(lmax)) {
      rb_raise(rb_eRuntimeError, "Vector length is too small. (%d for >= %d\n", (int) res->size,
               (int) alf_array_size(lmax));
    }
    ret1 = argv[2];
    ret2 = argv[3];
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of argumentso (%d for 1-4)\n", argc);
  }
  alf_Plm_deriv_array(lmax, x, res->data, deriv->data, w);
  ary = rb_ary_new2(2);
  rb_ary_store(ary, 0, ret1);
  rb_ary_store(ary, 1, ret2);
  return ary;
}

static VALUE rb_alf_array_size(VALUE module, VALUE lmax)
{
  return INT2FIX(alf_array_size(FIX2INT(lmax)));

}
static VALUE rb_alf_array_index(VALUE module, VALUE l, VALUE m)
{
  return INT2FIX(alf_array_index(FIX2INT(l), FIX2INT(m)));
}

void Init_alf(VALUE module)
{
  VALUE mALF;
  mALF = rb_define_module_under(module, "ALF");
  cWspace = rb_define_class_under(mALF, "Workspace", cGSL_Object);
  rb_define_singleton_method(cWspace, "alloc", rb_alf_alloc, 1);
  rb_define_singleton_method(mALF, "alloc", rb_alf_alloc, 1);
  rb_define_module_function(module, "alf_alloc", rb_alf_alloc, 1);

  rb_define_method(cWspace, "params", rb_alf_params, 3);

  rb_define_method(cWspace, "Plm_array", rb_alf_Plm_array, -1);
  rb_define_method(cWspace, "Plm_deriv_array", rb_alf_Plm_deriv_array, -1);

  rb_define_module_function(mALF, "array_size", rb_alf_array_size, 1);
  rb_define_module_function(mALF, "array_index", rb_alf_array_index, 2);

  define_constants(mALF);
}
#endif
