#include "include/rb_gsl.h"

static VALUE cWorkspace;

static VALUE rb_gsl_sf_mathieu_alloc(VALUE klass, VALUE n, VALUE q)
{
  gsl_sf_mathieu_workspace *w;
  w = gsl_sf_mathieu_alloc((size_t) FIX2INT(n), NUM2DBL(q));
  return Data_Wrap_Struct(klass, 0, gsl_sf_mathieu_free, w);
}

static VALUE sf_mathieu_eval(VALUE order, VALUE qq,
                             int (*f)(int, double, gsl_sf_result*))
{
  gsl_sf_result r;
  (*f)(FIX2INT(order), NUM2DBL(qq), &r);
  return rb_float_new(r.val);
}

static VALUE sf_mathieu_eval2(VALUE n1, VALUE n2, VALUE q, VALUE x,
                              int (*f)(int, int, double, double, gsl_sf_result*))
{
  gsl_sf_result r;
  (*f)(FIX2INT(n1),FIX2INT(n2),  NUM2DBL(q), NUM2DBL(x), &r);
  return rb_float_new(r.val);
}

static VALUE sf_mathieu_array_eval(int argc, VALUE *argv,
                                   int (*f)(int, int, double, gsl_sf_mathieu_workspace*, double[]))
{
  gsl_sf_mathieu_workspace *w;
  gsl_vector *v;
  int n1, n2;
  double q;
  switch (argc) {
  case 4:
    if (!rb_obj_is_kind_of(argv[3], cWorkspace)) {
      rb_raise(rb_eTypeError, "Wrong argument type 3 (%s detected, %s expected)",
               rb_class2name(CLASS_OF(argv[3])), rb_class2name(cWorkspace));
    }
    n1 = FIX2INT(argv[0]);
    n2 = FIX2INT(argv[1]);
    q = NUM2DBL(argv[2]);
    Data_Get_Struct(argv[3], gsl_sf_mathieu_workspace, w);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments. (%d for 4)", argc);
  }
  v = gsl_vector_alloc(n2 - n1 + 1);
  (*f)(n1, n2, q, w, v->data);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}

static VALUE sf_mathieu_array_eval2(int argc, VALUE *argv,
                                    int (*f)(int, int,  double, double, gsl_sf_mathieu_workspace*, double[]))
{
  gsl_sf_mathieu_workspace *w;
  gsl_vector *v;
  int n1, n2;
  double q, x;
  switch (argc) {
  case 5:
    if (!rb_obj_is_kind_of(argv[4], cWorkspace)) {
      rb_raise(rb_eTypeError, "Wrong argument type 4 (%s detected, %s expected)",
               rb_class2name(CLASS_OF(argv[4])), rb_class2name(cWorkspace));
    }
    n1 = FIX2INT(argv[0]);
    n2 = FIX2INT(argv[1]);
    q = NUM2DBL(argv[2]);
    x = NUM2DBL(argv[3]);
    Data_Get_Struct(argv[4], gsl_sf_mathieu_workspace, w);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments. (%d for 5)", argc);
  }
  v = gsl_vector_alloc(n2 - n1 + 1);
  (*f)(n1, n2, q, x, w, v->data);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}
static VALUE sf_mathieu_array_eval3(int argc, VALUE *argv,
                                    int (*f)(int, int, int, double, double, gsl_sf_mathieu_workspace*, double[]))
{
  gsl_sf_mathieu_workspace *w;
  gsl_vector *v;
  int n1, n2, n3;
  double q, x;
  switch (argc) {
  case 6:
    if (!rb_obj_is_kind_of(argv[5], cWorkspace)) {
      rb_raise(rb_eTypeError, "Wrong argument type 5 (%s detected, %s expected)",
               rb_class2name(CLASS_OF(argv[5])), rb_class2name(cWorkspace));
    }
    n1 = FIX2INT(argv[0]);
    n2 = FIX2INT(argv[1]);
    n3 = FIX2INT(argv[2]);
    q = NUM2DBL(argv[3]);
    x = NUM2DBL(argv[4]);
    Data_Get_Struct(argv[5], gsl_sf_mathieu_workspace, w);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments. (%d for 6)", argc);
  }
  v = gsl_vector_alloc(n3 - n2 + 1);
  (*f)(n1, n2, n3, q, x, w, v->data);
  return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, v);
}
static VALUE sf_mathieu_eval_int_double2(VALUE order, VALUE qq, VALUE zz,
                                         int (*f)(int, double, double, gsl_sf_result*))
{
  gsl_sf_result r;
  (*f)(FIX2INT(order), NUM2DBL(qq), NUM2DBL(zz), &r);
  return rb_float_new(r.val);
}
static VALUE sf_mathieu_eval_e_int_double2(VALUE order, VALUE qq, VALUE zz,
                                           int (*f)(int, double, double, gsl_sf_result*))
{
  gsl_sf_result *r;
  VALUE val;
  val = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, r);
  (*f)(FIX2INT(order), NUM2DBL(qq), NUM2DBL(zz), r);
  return val;
}

static VALUE sf_mathieu_eval_e_int2_double2(VALUE n1, VALUE n2, VALUE qq, VALUE zz,
                                            int (*f)(int, int, double, double, gsl_sf_result*))
{
  gsl_sf_result *r;
  VALUE val;
  val = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, r);
  (*f)(FIX2INT(n1), FIX2INT(n2), NUM2DBL(qq), NUM2DBL(zz), r);
  return val;
}
/**********/
static VALUE rb_gsl_sf_mathieu_a_e(VALUE module, VALUE order, VALUE qq)
{
#ifdef HAVE_GSL_SF_MATHIEU_A_E
  return rb_gsl_sf_eval_e_int_double(gsl_sf_mathieu_a_e, order, qq);
#else
  return rb_gsl_sf_eval_e_int_double(gsl_sf_mathieu_a, order, qq);
#endif
}
static VALUE rb_gsl_sf_mathieu_a(VALUE module, VALUE order, VALUE qq)
{
#ifdef HAVE_GSL_SF_MATHIEU_A_E
  return sf_mathieu_eval(order, qq, gsl_sf_mathieu_a_e);
#else
  return sf_mathieu_eval(order, qq, gsl_sf_mathieu_a);
#endif
}
static VALUE rb_gsl_sf_mathieu_a_array(VALUE module, int argc, VALUE *argv)
{
  return sf_mathieu_array_eval(argc, argv, gsl_sf_mathieu_a_array);
}
static VALUE rb_gsl_sf_mathieu_b_e(VALUE module, VALUE order, VALUE qq)
{
#ifdef HAVE_GSL_SF_MATHIEU_B_E
  return rb_gsl_sf_eval_e_int_double(gsl_sf_mathieu_b_e, order, qq);
#else
  return rb_gsl_sf_eval_e_int_double(gsl_sf_mathieu_b, order, qq);
#endif
}
static VALUE rb_gsl_sf_mathieu_b(VALUE module, VALUE order, VALUE qq)
{
#ifdef HAVE_GSL_SF_MATHIEU_B_E
  return sf_mathieu_eval(order, qq, gsl_sf_mathieu_b_e);
#else
  return sf_mathieu_eval(order, qq, gsl_sf_mathieu_b);
#endif
}
static VALUE rb_gsl_sf_mathieu_b_array(VALUE module, int argc, VALUE *argv)
{
  return sf_mathieu_array_eval(argc, argv, gsl_sf_mathieu_b_array);
}
static VALUE rb_gsl_sf_mathieu_ce_e(VALUE module, VALUE order, VALUE qq, VALUE zz)
{
#ifdef HAVE_GSL_SF_MATHIEU_CE_E
  return sf_mathieu_eval_e_int_double2(order, qq, zz, gsl_sf_mathieu_ce_e);
#else
  return sf_mathieu_eval_e_int_double2(order, qq, zz, gsl_sf_mathieu_ce);
#endif
}
static VALUE rb_gsl_sf_mathieu_ce(VALUE module, VALUE order, VALUE qq, VALUE zz)
{
#ifdef HAVE_GSL_SF_MATHIEU_CE_E
  return sf_mathieu_eval_int_double2(order, qq, zz, gsl_sf_mathieu_ce_e);
#else
  return sf_mathieu_eval_int_double2(order, qq, zz, gsl_sf_mathieu_ce);
#endif
}
static VALUE rb_gsl_sf_mathieu_ce_array(VALUE module, int argc, VALUE *argv)
{
  return sf_mathieu_array_eval2(argc, argv, gsl_sf_mathieu_ce_array);
}
static VALUE rb_gsl_sf_mathieu_se_e(VALUE module, VALUE order, VALUE qq, VALUE zz)
{
#ifdef HAVE_GSL_SF_MATHIEU_SE_E
  return sf_mathieu_eval_e_int_double2(order, qq, zz, gsl_sf_mathieu_se_e);
#else
  return sf_mathieu_eval_e_int_double2(order, qq, zz, gsl_sf_mathieu_se);
#endif
}
static VALUE rb_gsl_sf_mathieu_se(VALUE module, VALUE order, VALUE qq, VALUE zz)
{
#ifdef HAVE_GSL_SF_MATHIEU_SE_E
  return sf_mathieu_eval_int_double2(order, qq, zz, gsl_sf_mathieu_se_e);
#else
  return sf_mathieu_eval_int_double2(order, qq, zz, gsl_sf_mathieu_se);
#endif
}
static VALUE rb_gsl_sf_mathieu_se_array(VALUE module, int argc, VALUE *argv)
{
  return sf_mathieu_array_eval2(argc, argv, gsl_sf_mathieu_se_array);
}

/*****/
static VALUE rb_gsl_sf_mathieu_Mc_e(VALUE module, VALUE n1, VALUE n2, VALUE q, VALUE x)
{
#ifdef HAVE_GSL_SF_MATHIEU_MC_E
  return sf_mathieu_eval_e_int2_double2(n1, n2, q, x, gsl_sf_mathieu_Mc_e);
#else
  return sf_mathieu_eval_e_int2_double2(n1, n2, q, x, gsl_sf_mathieu_Mc);
#endif
}
static VALUE rb_gsl_sf_mathieu_Mc(VALUE module, VALUE n1, VALUE n2, VALUE q, VALUE x)
{
#ifdef HAVE_GSL_SF_MATHIEU_MC_E
  return sf_mathieu_eval2(n1, n2, q, x, gsl_sf_mathieu_Mc_e);
#else
  return sf_mathieu_eval2(n1, n2, q, x, gsl_sf_mathieu_Mc);
#endif
}
static VALUE rb_gsl_sf_mathieu_Mc_array(VALUE module, int argc, VALUE *argv)
{
  return sf_mathieu_array_eval3(argc, argv, gsl_sf_mathieu_Mc_array);
}
static VALUE rb_gsl_sf_mathieu_Ms_e(VALUE module, VALUE n1, VALUE n2, VALUE q, VALUE x)
{
#ifdef HAVE_GSL_SF_MATHIEU_MS_E
  return sf_mathieu_eval_e_int2_double2(n1, n2, q, x, gsl_sf_mathieu_Ms_e);
#else
  return sf_mathieu_eval_e_int2_double2(n1, n2, q, x, gsl_sf_mathieu_Ms);
#endif
}
static VALUE rb_gsl_sf_mathieu_Ms(VALUE module, VALUE n1, VALUE n2, VALUE q, VALUE x)
{
#ifdef HAVE_GSL_SF_MATHIEU_MS_E
  return sf_mathieu_eval2(n1, n2, q, x, gsl_sf_mathieu_Ms_e);
#else
  return sf_mathieu_eval2(n1, n2, q, x, gsl_sf_mathieu_Ms);
#endif
}
static VALUE rb_gsl_sf_mathieu_Ms_array(VALUE module, int argc, VALUE *argv)
{
  return sf_mathieu_array_eval3(argc, argv, gsl_sf_mathieu_Ms_array);
}
/*****/
void Init_sf_mathieu(VALUE module)
{
  VALUE mMathieu;

  mMathieu = rb_define_module_under(module, "Mathieu");
  cWorkspace = rb_define_class_under(mMathieu, "Workspace", cGSL_Object);
  rb_define_singleton_method(cWorkspace, "alloc", rb_gsl_sf_mathieu_alloc, 2);

  rb_define_module_function(module, "mathieu_a", rb_gsl_sf_mathieu_a, 2);
  rb_define_module_function(module, "mathieu_a_e", rb_gsl_sf_mathieu_a_e, 2);
  rb_define_module_function(module, "mathieu_a_array", rb_gsl_sf_mathieu_a_array, -1);
  rb_define_module_function(module, "mathieu_b", rb_gsl_sf_mathieu_b, 2);
  rb_define_module_function(module, "mathieu_b_e", rb_gsl_sf_mathieu_b_e, 2);
  rb_define_module_function(module, "mathieu_b_array", rb_gsl_sf_mathieu_b_array, -1);
  rb_define_module_function(module, "mathieu_ce", rb_gsl_sf_mathieu_ce, 3);
  rb_define_module_function(module, "mathieu_ce_e", rb_gsl_sf_mathieu_ce_e, 3);
  rb_define_module_function(module, "mathieu_ce_array", rb_gsl_sf_mathieu_ce_array, -1);
  rb_define_module_function(module, "mathieu_se", rb_gsl_sf_mathieu_se, 3);
  rb_define_module_function(module, "mathieu_se_e", rb_gsl_sf_mathieu_se_e, 3);
  rb_define_module_function(module, "mathieu_se_array", rb_gsl_sf_mathieu_se_array, -1);
  rb_define_module_function(module, "mathieu_Mc", rb_gsl_sf_mathieu_Mc, 4);
  rb_define_module_function(module, "mathieu_Mc_e", rb_gsl_sf_mathieu_Mc_e, 4);
  rb_define_module_function(module, "mathieu_Mc_array", rb_gsl_sf_mathieu_Mc_array, -1);
  rb_define_module_function(module, "mathieu_Ms", rb_gsl_sf_mathieu_Ms, 4);
  rb_define_module_function(module, "mathieu_Ms_e", rb_gsl_sf_mathieu_Ms_e, 4);
  rb_define_module_function(module, "mathieu_Ms_array", rb_gsl_sf_mathieu_Ms_array, -1);
}
