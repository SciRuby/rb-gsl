#ifdef HAVE_GSL_GSL_CQP_H

#include "include/rb_gsl.h"
#include "gsl/gsl_cqp.h"

static const gsl_cqpminimizer_type* type_by_string(VALUE t);
static const gsl_cqpminimizer_type* get_type(VALUE t)
{

  switch (TYPE(t)) {
  case T_STRING:
    return type_by_string(t);
    break;
  default:
    rb_raise(rb_eTypeError, "Wrong argument type %s.", rb_class2name(CLASS_OF(t)));
  }
}

static const gsl_cqpminimizer_type* type_by_string(VALUE t)
{
  char *name;
  name = STR2CSTR(t);
  if (strcmp(name, "mg_pdip") == 0) {
    return gsl_cqpminimizer_mg_pdip;
  } else {
    rb_raise(rb_eRuntimeError, "Unknown minimizer type %s.", name);
  }
  return NULL; /* never reach here */
}

static VALUE rb_cqpminimizer_alloc(VALUE klass,  VALUE t, VALUE n, VALUE me, VALUE mi)
{
  gsl_cqpminimizer *m;
  m = gsl_cqpminimizer_alloc(get_type(t), (size_t) FIX2INT(n), (size_t) FIX2INT(me), (size_t) FIX2INT(mi));
  return Data_Wrap_Struct(klass, 0, gsl_cqpminimizer_free, m);
}

static VALUE rb_cqpminimizer_set(VALUE obj, VALUE data)
{
  gsl_cqpminimizer *m;
  gsl_cqp_data *d;
  Data_Get_Struct(obj, gsl_cqpminimizer, m);
  Data_Get_Struct(data, gsl_cqp_data, d);
  gsl_cqpminimizer_set(m, d);
  return Qtrue;
}

static VALUE rb_cqpminimizer_name(VALUE obj)
{
  gsl_cqpminimizer *m;
  Data_Get_Struct(obj, gsl_cqpminimizer, m);
  return rb_str_new2(gsl_cqpminimizer_name(m));
}

static VALUE rb_cqpminimizer_iterate(VALUE obj)
{
  gsl_cqpminimizer *m;
  Data_Get_Struct(obj, gsl_cqpminimizer, m);
  return INT2FIX(gsl_cqpminimizer_iterate(m));
}

static VALUE rb_cqpminimizer_x(VALUE obj)
{
  gsl_cqpminimizer *m;
  Data_Get_Struct(obj, gsl_cqpminimizer, m);
  return Data_Wrap_Struct(cgsl_vector_view, 0, NULL, gsl_cqpminimizer_x(m));
}

static VALUE rb_cqpminimizer_lm_eq(VALUE obj)
{
  gsl_cqpminimizer *m;
  Data_Get_Struct(obj, gsl_cqpminimizer, m);
  return Data_Wrap_Struct(cgsl_vector_view, 0, NULL, gsl_cqpminimizer_lm_eq(m));
}
static VALUE rb_cqpminimizer_lm_ineq(VALUE obj)
{
  gsl_cqpminimizer *m;
  Data_Get_Struct(obj, gsl_cqpminimizer, m);
  return Data_Wrap_Struct(cgsl_vector_view, 0, NULL, gsl_cqpminimizer_lm_ineq(m));
}
static VALUE rb_cqpminimizer_f(VALUE obj)
{
  gsl_cqpminimizer *m;
  Data_Get_Struct(obj, gsl_cqpminimizer, m);
  return rb_float_new(gsl_cqpminimizer_f(m));
}
static VALUE rb_cqpminimizer_gap(VALUE obj)
{
  gsl_cqpminimizer *m;
  Data_Get_Struct(obj, gsl_cqpminimizer, m);
  return rb_float_new(gsl_cqpminimizer_gap(m));
}
static VALUE rb_cqpminimizer_residuals_norm(VALUE obj)
{
  gsl_cqpminimizer *m;
  Data_Get_Struct(obj, gsl_cqpminimizer, m);
  return rb_float_new(gsl_cqpminimizer_residuals_norm(m));
}
/*
static VALUE rb_cqpminimizer_minimum(VALUE obj)
{
  gsl_cqpminimizer *m;
  Data_Get_Struct(obj, gsl_cqpminimizer, m);
  return rb_float_new(gsl_cqpminimizer_minimum(m));
}
*/
static VALUE rb_cqpminimizer_test_convergence(VALUE obj, VALUE g, VALUE r)
{
  gsl_cqpminimizer *m;
  Data_Get_Struct(obj, gsl_cqpminimizer, m);
  return INT2FIX(gsl_cqpminimizer_test_convergence(m, NUM2DBL(g), NUM2DBL(r)));
}
static VALUE rb_cqpminimizer_test_infeasibility(VALUE obj, VALUE e)
{
  gsl_cqpminimizer *m;
  Data_Get_Struct(obj, gsl_cqpminimizer, m);
  return INT2FIX(gsl_cqp_minimizer_test_infeasibility(m, NUM2DBL(e)));
}

static VALUE rb_cqp_data_alloc(VALUE klass)
{
  gsl_cqp_data *d;
  d = (gsl_cqp_data*) malloc(sizeof(gsl_cqp_data));
  return Data_Wrap_Struct(klass, 0, free, d);
}

static VALUE rb_cqp_data_Q(VALUE obj)
{
  gsl_cqp_data *d;
  Data_Get_Struct(obj, gsl_cqp_data, d);
  return Data_Wrap_Struct(cgsl_matrix_view, 0, NULL,   d->Q);
}

static VALUE rb_cqp_data_q(VALUE obj)
{
  gsl_cqp_data *d;
  Data_Get_Struct(obj, gsl_cqp_data, d);
  return Data_Wrap_Struct(cgsl_vector_view, 0, NULL,   d->q);
}

static VALUE rb_cqp_data_A(VALUE obj)
{
  gsl_cqp_data *d;
  Data_Get_Struct(obj, gsl_cqp_data, d);
  return Data_Wrap_Struct(cgsl_matrix_view, 0, NULL,   d->A);
}

static VALUE rb_cqp_data_b(VALUE obj)
{
  gsl_cqp_data *d;
  Data_Get_Struct(obj, gsl_cqp_data, d);
  return Data_Wrap_Struct(cgsl_vector_view, 0, NULL,   d->b);
}

static VALUE rb_cqp_data_C(VALUE obj)
{
  gsl_cqp_data *d;
  Data_Get_Struct(obj, gsl_cqp_data, d);
  return Data_Wrap_Struct(cgsl_matrix_view, 0, NULL,   d->C);
}

static VALUE rb_cqp_data_d(VALUE obj)
{
  gsl_cqp_data *d;
  Data_Get_Struct(obj, gsl_cqp_data, d);
  return Data_Wrap_Struct(cgsl_vector_view, 0, NULL,   d->d);
}

static VALUE rb_cqp_data_set_Q(VALUE obj, VALUE mm)
{
  gsl_cqp_data *d;
  gsl_matrix *m;
  Data_Get_Struct(obj, gsl_cqp_data, d);
  CHECK_MATRIX(mm);
  Data_Get_Struct(mm, gsl_matrix, m);
  d->Q = m;
  return Qtrue;
}

static VALUE rb_cqp_data_set_q(VALUE obj, VALUE vv)
{
  gsl_cqp_data *d;
  gsl_vector *v;
  Data_Get_Struct(obj, gsl_cqp_data, d);
  CHECK_VECTOR(vv);
  Data_Get_Struct(vv, gsl_vector, v);
  d->q = v;
  return Qtrue;
}

static VALUE rb_cqp_data_set_A(VALUE obj, VALUE mm)
{
  gsl_cqp_data *d;
  gsl_matrix *m;
  Data_Get_Struct(obj, gsl_cqp_data, d);
  CHECK_MATRIX(mm);
  Data_Get_Struct(mm, gsl_matrix, m);
  d->A = m;
  return Qtrue;
}

static VALUE rb_cqp_data_set_b(VALUE obj, VALUE vv)
{
  gsl_cqp_data *d;
  gsl_vector *v;
  Data_Get_Struct(obj, gsl_cqp_data, d);
  CHECK_VECTOR(vv);
  Data_Get_Struct(vv, gsl_vector, v);
  d->b = v;
  return Qtrue;
}

static VALUE rb_cqp_data_set_C(VALUE obj, VALUE mm)
{
  gsl_cqp_data *d;
  gsl_matrix *m;
  Data_Get_Struct(obj, gsl_cqp_data, d);
  CHECK_MATRIX(mm);
  Data_Get_Struct(mm, gsl_matrix, m);
  d->C = m;
  return Qtrue;
}

static VALUE rb_cqp_data_set_d(VALUE obj, VALUE vv)
{
  gsl_cqp_data *d;
  gsl_vector *v;
  Data_Get_Struct(obj, gsl_cqp_data, d);
  CHECK_VECTOR(vv);
  Data_Get_Struct(vv, gsl_vector, v);
  d->d = v;
  return Qtrue;
}

void Init_cqp(VALUE module)
{
  VALUE mCQP, cMinimizer, cData;

  mCQP = rb_define_module_under(module, "CQP");
  cMinimizer = rb_define_class_under(mCQP, "Minimizer", cGSL_Object);
  cData = rb_define_class_under(mCQP, "Data", cGSL_Object);

  rb_define_singleton_method(cMinimizer, "alloc", rb_cqpminimizer_alloc, 4);

  rb_define_method(cMinimizer, "set", rb_cqpminimizer_set, 1);
  rb_define_method(cMinimizer, "name", rb_cqpminimizer_name, 0);
  rb_define_method(cMinimizer, "iterate", rb_cqpminimizer_iterate, 0);
  rb_define_method(cMinimizer, "x", rb_cqpminimizer_x, 0);
  rb_define_method(cMinimizer, "lm_eq", rb_cqpminimizer_lm_eq, 0);
  rb_define_method(cMinimizer, "lm_ineq", rb_cqpminimizer_lm_ineq, 0);
  rb_define_method(cMinimizer, "f", rb_cqpminimizer_f, 0);
  rb_define_method(cMinimizer, "gap", rb_cqpminimizer_gap, 0);
  rb_define_method(cMinimizer, "residuals_norm", rb_cqpminimizer_residuals_norm, 0);
/*  rb_define_method(cMinimizer, "minimum", rb_cqpminimizer_minimum, 0);  */
  rb_define_method(cMinimizer, "test_convergence", rb_cqpminimizer_test_convergence, 2);
  rb_define_method(cMinimizer, "test_infeasibility", rb_cqpminimizer_test_infeasibility, 1);

  /*****/
  rb_define_singleton_method(cData, "alloc", rb_cqp_data_alloc, 0);

  rb_define_method(cData, "Q", rb_cqp_data_Q, 0);
  rb_define_method(cData, "q", rb_cqp_data_q, 0);
  rb_define_method(cData, "A", rb_cqp_data_A, 0);
  rb_define_method(cData, "b", rb_cqp_data_b, 0);
  rb_define_method(cData, "C", rb_cqp_data_C, 0);
  rb_define_method(cData, "d", rb_cqp_data_d, 0);

  rb_define_method(cData, "set_Q", rb_cqp_data_set_Q, 1);
  rb_define_method(cData, "set_q", rb_cqp_data_set_q, 1);
  rb_define_method(cData, "set_A", rb_cqp_data_set_A, 1);
  rb_define_method(cData, "set_b", rb_cqp_data_set_b, 1);
  rb_define_method(cData, "set_C", rb_cqp_data_set_C, 1);
  rb_define_method(cData, "set_d", rb_cqp_data_set_d, 1);
  rb_define_alias(cData, "Q=", "set_Q");
  rb_define_alias(cData, "q=", "set_q");
  rb_define_alias(cData, "A=", "set_A");
  rb_define_alias(cData, "b=", "set_b");
  rb_define_alias(cData, "C=", "set_C");
  rb_define_alias(cData, "d=", "set_d");
}

#endif

