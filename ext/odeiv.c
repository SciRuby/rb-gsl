/*
  odeiv.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/
#include "rb_gsl_config.h"
#include "rb_gsl_odeiv.h"
#include "rb_gsl_array.h"
#include "rb_gsl_function.h"

#ifndef CHECK_SYSTEM
#define CHECK_SYSTEM(x) if(CLASS_OF(x)!=cgsl_odeiv_system)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (GSL::Odeiv::System expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

#ifndef CHECK_STEP
#define CHECK_STEP(x) if(CLASS_OF(x)!=cgsl_odeiv_step)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (GSL::Odeiv::Step expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

#ifndef CHECK_CONTROL
#define CHECK_CONTROL(x) if(CLASS_OF(x)!=cgsl_odeiv_control)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (GSL::Odeiv::Control expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

#ifndef CHECK_EVOLVE
#define CHECK_EVOLVE(x) if(CLASS_OF(x)!=cgsl_odeiv_evolve)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (GSL::Odeiv::Evolve expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

#ifndef CHECK_SOLVER
#define CHECK_SOLVER(x) if(CLASS_OF(x)!=cgsl_odeiv_solver)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (GSL::Odeiv::Solver expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

static VALUE cgsl_odeiv_system;
static VALUE cgsl_odeiv_step;
static VALUE cgsl_odeiv_control;
static VALUE cgsl_odeiv_evolve;
static VALUE cgsl_odeiv_solver;

enum {
  GSL_ODEIV_STEP_RK2,
  GSL_ODEIV_STEP_RK4,
  GSL_ODEIV_STEP_RKF45,
  GSL_ODEIV_STEP_RKCK,
  GSL_ODEIV_STEP_RK8PD,
  GSL_ODEIV_STEP_RK2IMP,
  GSL_ODEIV_STEP_RK4IMP,
  GSL_ODEIV_STEP_BSIMP,
  GSL_ODEIV_STEP_GEAR1,
  GSL_ODEIV_STEP_GEAR2,
  GSL_ODEIV_STEP_RK2SIMP,
};

typedef struct __rb_gsl_odeiv_evolve {
  gsl_odeiv_evolve *e;
  gsl_odeiv_control *c;
  gsl_odeiv_step *s;
  gsl_odeiv_system *sys;
} gsl_odeiv_solver;

static int calc_func(double t, const double y[], double dydt[], void *data);
static int calc_jac(double t, const double y[], double *dfdy, double dfdt[], void *data);

static int calc_func(double t, const double y[], double dydt[], void *data)
{
  VALUE ary, params, proc;
  // local variable "result" declared and set, but never used
  //VALUE result;
  VALUE vy, vdydt;
  gsl_vector_view ytmp, dydttmp;
  size_t dim;

  ary = (VALUE) data;
  proc = rb_ary_entry(ary, 0);
  dim = FIX2INT(rb_ary_entry(ary, 2));
  params = rb_ary_entry(ary, 3);

  ytmp.vector.data = (double *) y;
  ytmp.vector.stride = 1;
  ytmp.vector.size = dim;
  dydttmp.vector.data = dydt;
  dydttmp.vector.stride = 1;
  dydttmp.vector.size = dim;
  vy = Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, &ytmp);
  vdydt = Data_Wrap_Struct(cgsl_vector_view, 0, NULL, &dydttmp);

  if (NIL_P(params)) /*result =*/ rb_funcall((VALUE) proc, RBGSL_ID_call, 3, rb_float_new(t), 
					 vy, vdydt);
  else /*result =*/ rb_funcall((VALUE) proc, RBGSL_ID_call, 4, rb_float_new(t), vy, vdydt, params);
 
  return GSL_SUCCESS;
}

static int calc_jac(double t, const double y[], double *dfdy, double dfdt[], void *data)
{
  VALUE params, proc, ary;
  VALUE vdfdt;
  // local variable "result" declared and set, but never used
  //VALUE result;
  VALUE vy, vmjac;
  gsl_vector_view ytmp, dfdttmp;
  gsl_matrix_view mv;
  size_t dim;
  
  ary = (VALUE) data;
  proc = rb_ary_entry(ary, 1);
  if (NIL_P(proc)) rb_raise(rb_eRuntimeError, "df function not given");

  dim = FIX2INT(rb_ary_entry(ary, 2));
  params = rb_ary_entry(ary, 3);

  ytmp.vector.data = (double *) y;
  ytmp.vector.size = dim;
  ytmp.vector.stride = 1;
  dfdttmp.vector.data = dfdt;
  dfdttmp.vector.size = dim;
  dfdttmp.vector.stride = 1;
  mv = gsl_matrix_view_array(dfdy, dim, dim);
  vy = Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, &ytmp);
  vmjac = Data_Wrap_Struct(cgsl_matrix_view, 0, NULL, &mv);
  vdfdt = Data_Wrap_Struct(cgsl_vector_view, 0, NULL, &dfdttmp);
  if (NIL_P(params)) /*result =*/ rb_funcall((VALUE) proc, RBGSL_ID_call, 4, rb_float_new(t),
					 vy, vmjac, vdfdt);
  else /*result =*/ rb_funcall((VALUE) proc, RBGSL_ID_call, 5, rb_float_new(t), 
			   vy, vmjac, vdfdt, params);
  return GSL_SUCCESS;
}

static void gsl_odeiv_system_mark(gsl_odeiv_system *sys);
static void gsl_odeiv_system_mark(gsl_odeiv_system *sys)
{
  rb_gc_mark((VALUE) sys->params);
}

static gsl_odeiv_system* make_sys(int argc, VALUE *argv);
static void set_sys(int argc, VALUE *argv, gsl_odeiv_system *sys);
static VALUE rb_gsl_odeiv_system_new(int argc, VALUE *argv, VALUE klass)
{
  VALUE obj;
  gsl_odeiv_system *sys = NULL;
  sys = make_sys(argc, argv);
  obj = Data_Wrap_Struct(klass, gsl_odeiv_system_mark, free, sys);
  return obj;
}

static VALUE rb_gsl_odeiv_system_set(int argc, VALUE *argv, VALUE obj)
{
  gsl_odeiv_system *sys = NULL;
  Data_Get_Struct(obj, gsl_odeiv_system, sys);
  set_sys(argc, argv, sys);
  return obj;
}

static gsl_odeiv_system* make_sys(int argc, VALUE *argv)
{
  gsl_odeiv_system *sys = NULL;
  sys = ALLOC(gsl_odeiv_system);
  sys->function = &calc_func;
  sys->jacobian = &calc_jac;
  sys->params = NULL;
  set_sys(argc, argv, sys);
  return sys;
}

static void set_sys(int argc, VALUE *argv, gsl_odeiv_system *sys)
{
  size_t dimension;
  VALUE ary, vjac, dim;
  VALUE vparams;
  int itmp;
  size_t i, j;

  if (argc < 2) rb_raise(rb_eArgError, "too few arguments");
  CHECK_PROC(argv[0]);

  if (sys == NULL) {
    sys = ALLOC(gsl_odeiv_system);
    sys->function = &calc_func;
    sys->jacobian = &calc_jac;
  }

  if (sys->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) sys->params = ary;*/
    sys->params = (void *) ary;
  } else {
    ary = (VALUE) sys->params;
  }
  rb_ary_store(ary, 1, Qnil);   /* function to calc J */
  rb_ary_store(ary, 3, Qnil);   /* parameters */

  itmp = 1;
  if (rb_obj_is_kind_of(argv[1], rb_cProc)) {
    vjac = argv[1];
    itmp = 2;
  } else {
    vjac = Qnil;
  }
  if ((dim =argv[itmp++]) == Qnil) dim = argv[itmp++];
  switch (argc - itmp) {
  case 0:
    vparams = Qnil;
    break;
  case 1:
    vparams = argv[itmp];
    break;
  default:
    vparams = rb_ary_new2(argc-itmp);
    for (i = itmp, j = 0; (int) i < argc; i++, j++) rb_ary_store(vparams, j, argv[i]);
  }
  dimension = FIX2INT(dim);
  sys->dimension = dimension;
  rb_ary_store(ary, 0, argv[0]); 
  rb_ary_store(ary, 1, vjac);
  rb_ary_store(ary, 2, dim);  
  rb_ary_store(ary, 3, vparams); 
}

static VALUE rb_gsl_odeiv_system_set_params(int argc, VALUE *argv, VALUE obj);
static VALUE rb_gsl_odeiv_system_set_params(int argc, VALUE *argv, VALUE obj)
{
  VALUE vparams, ary;
  gsl_odeiv_system *sys = NULL;
  size_t i;
  Data_Get_Struct(obj, gsl_odeiv_system, sys);

  ary = (VALUE) sys->params;
  switch (argc) {
  case 0:
    vparams = Qnil;
    break;
  case 1:
    vparams = argv[0];
    break;
  default:
    vparams = rb_ary_new2(argc);
    for (i = 0; (int) i < argc; i++) rb_ary_store(vparams, i, argv[i]);
  }
  //  rb_ary_delete_at(ary, 3);
  rb_ary_store(ary, 3, vparams); 
  return obj;
}

static VALUE rb_gsl_odeiv_system_params(VALUE obj)
{
  VALUE ary;
  gsl_odeiv_system *sys = NULL;
  Data_Get_Struct(obj, gsl_odeiv_system, sys);
  ary = (VALUE) sys->params;
  return rb_ary_entry(ary, 3);
}

static VALUE rb_gsl_odeiv_system_function(VALUE obj)
{
  VALUE ary;
  gsl_odeiv_system *sys = NULL;
  Data_Get_Struct(obj, gsl_odeiv_system, sys);
  ary = (VALUE) sys->params;
  return rb_ary_entry(ary, 0);
}

static VALUE rb_gsl_odeiv_system_jacobian(VALUE obj)
{
  VALUE ary;
  gsl_odeiv_system *sys = NULL;
  Data_Get_Struct(obj, gsl_odeiv_system, sys);
  ary = (VALUE) sys->params;
  return rb_ary_entry(ary, 1);
}

static VALUE rb_gsl_odeiv_system_dimension(VALUE obj)
{
  gsl_odeiv_system *sys = NULL;
  Data_Get_Struct(obj, gsl_odeiv_system, sys);
  return INT2FIX(sys->dimension);
}

static const gsl_odeiv_step_type* rb_gsl_odeiv_step_type_get(VALUE tt);

static gsl_odeiv_step* make_step(VALUE tt, VALUE dim);
static VALUE rb_gsl_odeiv_step_new(int argc, VALUE *argv, VALUE klass)
{
  VALUE obj;
  gsl_odeiv_step *s = NULL;
  switch (argc) {
  case 1:
    CHECK_FIXNUM(argv[0]);
    s = make_step(INT2FIX(GSL_ODEIV_STEP_RKF45), argv[0]);
    break;
  case 2:
    CHECK_FIXNUM(argv[1]);
    s = make_step(argv[0], argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 1 or 2)", argc);
  }
  obj = Data_Wrap_Struct(klass, 0, gsl_odeiv_step_free, s);
  return obj;
}

static gsl_odeiv_step* make_step(VALUE tt, VALUE dim)
{
  const gsl_odeiv_step_type *T;
  T = rb_gsl_odeiv_step_type_get(tt);
  return gsl_odeiv_step_alloc(T, FIX2INT(dim));
}

static const gsl_odeiv_step_type* rb_gsl_odeiv_step_type_get(VALUE tt)
{
  const gsl_odeiv_step_type *T;
  int type;
  char name[64];
  switch (TYPE(tt)) {
  case T_FIXNUM:
    type = FIX2INT(tt);
    switch (type) {
    case GSL_ODEIV_STEP_RK2: T = gsl_odeiv_step_rk2; break;
    case GSL_ODEIV_STEP_RK4: T = gsl_odeiv_step_rk4; break;
    case GSL_ODEIV_STEP_RKF45: T = gsl_odeiv_step_rkf45; break;
    case GSL_ODEIV_STEP_RKCK: T = gsl_odeiv_step_rkck; break;
    case GSL_ODEIV_STEP_RK8PD: T = gsl_odeiv_step_rk8pd; break;
    case GSL_ODEIV_STEP_RK2IMP: T = gsl_odeiv_step_rk2imp; break;
    case GSL_ODEIV_STEP_RK4IMP: T = gsl_odeiv_step_rk4imp; break;
    case GSL_ODEIV_STEP_BSIMP: T = gsl_odeiv_step_bsimp; break;
    case GSL_ODEIV_STEP_GEAR1: T = gsl_odeiv_step_gear1; break;
    case GSL_ODEIV_STEP_GEAR2: T = gsl_odeiv_step_gear2; break;
#ifdef GSL_1_6_LATER
    case GSL_ODEIV_STEP_RK2SIMP: T = gsl_odeiv_step_rk2simp; break;
#endif
    default:
      rb_raise(rb_eArgError, "wrong argument type (Fixnum expected)");
      break;
    }
    break;
  case T_STRING:
    strcpy(name, STR2CSTR(tt));
    if (str_tail_grep(name, "rk2") == 0) T = gsl_odeiv_step_rk2;
    else if (str_tail_grep(name, "rk4") == 0) T = gsl_odeiv_step_rk4;
    else if (str_tail_grep(name, "rkf45") == 0) T = gsl_odeiv_step_rkf45;
    else if (str_tail_grep(name, "rkck") == 0) T = gsl_odeiv_step_rkck;
    else if (str_tail_grep(name, "rk8pd") == 0) T = gsl_odeiv_step_rk8pd;
    else if (str_tail_grep(name, "rk2imp") == 0) T = gsl_odeiv_step_rk2imp;
    else if (str_tail_grep(name, "rk4imp") == 0) T = gsl_odeiv_step_rk4imp;
    else if (str_tail_grep(name, "bsimp") == 0) T = gsl_odeiv_step_bsimp;
    else if (str_tail_grep(name, "gear1") == 0) T = gsl_odeiv_step_gear1;
    else if (str_tail_grep(name, "gear2") == 0) T = gsl_odeiv_step_gear2;
#ifdef GSL_1_6_LATER
    else if (str_tail_grep(name, "rk2simp") == 0) T = gsl_odeiv_step_rk2simp;
#endif
    else {
      rb_raise(rb_eArgError, "wrong argument type %s", name);
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong argument type %s (String or Fixnum expected)",
	     rb_class2name(CLASS_OF(tt)));
    break;
  }
  return T;
}

static VALUE rb_gsl_odeiv_step_reset(VALUE obj)
{
  gsl_odeiv_step *s = NULL;
  Data_Get_Struct(obj, gsl_odeiv_step, s);
  return INT2FIX(gsl_odeiv_step_reset(s));
}

static VALUE rb_gsl_odeiv_step_name(VALUE obj)
{
  gsl_odeiv_step *s = NULL;
  Data_Get_Struct(obj, gsl_odeiv_step, s);
  return rb_str_new2(gsl_odeiv_step_name(s));
}

static VALUE rb_gsl_odeiv_step_order(VALUE obj)
{
  gsl_odeiv_step *s = NULL;
  Data_Get_Struct(obj, gsl_odeiv_step, s);
  return INT2FIX(gsl_odeiv_step_order(s));
}

static VALUE rb_gsl_odeiv_step_dimension(VALUE obj)
{
  gsl_odeiv_step *s = NULL;
  Data_Get_Struct(obj, gsl_odeiv_step, s);
  return INT2FIX(s->dimension);
}

static VALUE rb_gsl_odeiv_step_apply(int argc, VALUE *argv, VALUE obj)
{
  gsl_odeiv_step *s = NULL;
  gsl_odeiv_system *sys = NULL;
  gsl_vector *y = NULL, *yerr = NULL;
  gsl_vector *vtmp1 = NULL, *vtmp2 = NULL;
  double *dydt_in = NULL, *dydt_out = NULL;
  double t, h;
  switch (argc) {
  case 5:
    break;
  case 7:
    if (VECTOR_P(argv[5])) {
      Data_Get_Struct(argv[5], gsl_vector, vtmp2);
      if (vtmp2) dydt_out = vtmp2->data;
    }
    /* no break */
  case 6:
    if (VECTOR_P(argv[4])) {
      Data_Get_Struct(argv[4], gsl_vector, vtmp1);
      if (vtmp1) dydt_in = vtmp1->data;
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 5, 6 or 7)", argc);
    break;
  }
  Need_Float(argv[0]); Need_Float(argv[1]);
  CHECK_VECTOR(argv[2]); CHECK_VECTOR(argv[3]);
  CHECK_SYSTEM(argv[argc-1]);
  Data_Get_Struct(obj, gsl_odeiv_step, s);
  t = NUM2DBL(argv[0]);
  h = NUM2DBL(argv[1]);
  Data_Get_Struct(argv[2], gsl_vector, y);
  Data_Get_Struct(argv[3], gsl_vector, yerr);
  Data_Get_Struct(argv[argc-1], gsl_odeiv_system, sys);
  return INT2FIX(gsl_odeiv_step_apply(s, t, h, y->data, yerr->data, 
				      dydt_in, dydt_out, sys));
}

static VALUE rb_gsl_odeiv_step_info(VALUE obj)
{
  gsl_odeiv_step *s;
  char buf[256];
  Data_Get_Struct(obj, gsl_odeiv_step, s);
  sprintf(buf, "Class:      %s\n", rb_class2name(CLASS_OF(obj)));
  sprintf(buf, "%sSuperClass: %s\n", buf, rb_class2name(RCLASS_SUPER(CLASS_OF(obj))));
  sprintf(buf, "%sType:       %s\n", buf, gsl_odeiv_step_name(s));
  sprintf(buf, "%sDimension:  %d\n", buf, (int) s->dimension);
  return rb_str_new2(buf);
}

static gsl_odeiv_control* make_control_standard(VALUE epsabs, 
					       VALUE epsrel,
						VALUE ay, VALUE adydt);
static gsl_odeiv_control* make_control_y(VALUE epsabs, VALUE epsrel);
static VALUE rb_gsl_odeiv_control_standard_new(VALUE klass, VALUE epsabs, 
					       VALUE epsrel,
					       VALUE ay, VALUE adydt)
{
  gsl_odeiv_control *c = NULL;
  c = make_control_standard(epsabs, epsrel, ay, adydt); 
  return Data_Wrap_Struct(klass, 0, gsl_odeiv_control_free, c);
}

static gsl_odeiv_control* make_control_standard(VALUE epsabs, 
					       VALUE epsrel,
					       VALUE ay, VALUE adydt)
{
  Need_Float(epsabs); Need_Float(epsrel);
  Need_Float(ay);   Need_Float(adydt);
  return gsl_odeiv_control_standard_new(NUM2DBL(epsabs), NUM2DBL(epsrel), 
					NUM2DBL(ay), NUM2DBL(adydt));
}

static gsl_odeiv_control* make_control_y(VALUE epsabs, VALUE epsrel)
{
  Need_Float(epsabs); Need_Float(epsrel);
  return gsl_odeiv_control_y_new(NUM2DBL(epsabs), NUM2DBL(epsrel));
}

static VALUE rb_gsl_odeiv_control_y_new(VALUE klass, VALUE epsabs, 
					VALUE epsrel)
{
  gsl_odeiv_control *c = NULL;
  c = make_control_y(epsabs, epsrel);
  return Data_Wrap_Struct(klass, 0, gsl_odeiv_control_free, c);
}

static VALUE rb_gsl_odeiv_control_yp_new(VALUE klass, VALUE epsabs, 
					 VALUE epsrel)
{
  gsl_odeiv_control *c = NULL;
  Need_Float(epsabs); Need_Float(epsrel);
  c = gsl_odeiv_control_yp_new(NUM2DBL(epsabs), NUM2DBL(epsrel));
  return Data_Wrap_Struct(klass, 0, gsl_odeiv_control_free, c);
}

#ifdef GSL_1_2_LATER
static VALUE rb_gsl_odeiv_control_scaled_new(VALUE klass, VALUE epsabs, 
					     VALUE epsrel,
					     VALUE ay, VALUE adydt,
					     VALUE sc, VALUE dd)
{
  gsl_odeiv_control *c = NULL;
  gsl_vector *v = NULL;
  Need_Float(epsabs); Need_Float(epsrel);
  Need_Float(ay);   Need_Float(adydt);
  CHECK_FIXNUM(dd);
  CHECK_VECTOR(sc);
  Data_Get_Struct(sc, gsl_vector, v);
  c = gsl_odeiv_control_scaled_new(NUM2DBL(epsabs), NUM2DBL(epsrel), 
				   NUM2DBL(ay), NUM2DBL(adydt), v->data,
				   FIX2INT(dd));
  return Data_Wrap_Struct(klass, 0, gsl_odeiv_control_free, c);
}
#endif

static VALUE rb_gsl_odeiv_control_init(VALUE obj, VALUE epsabs, 
				       VALUE epsrel,
				       VALUE ay, VALUE adydt)
{
  gsl_odeiv_control *c = NULL;
  Need_Float(epsabs); Need_Float(epsrel);
  Need_Float(ay);   Need_Float(adydt);
  Data_Get_Struct(obj, gsl_odeiv_control, c);
  gsl_odeiv_control_init(c, NUM2DBL(epsabs), NUM2DBL(epsrel), 
			 NUM2DBL(ay), NUM2DBL(adydt));
  return obj;
}

static VALUE rb_gsl_odeiv_control_name(VALUE obj)
{
  gsl_odeiv_control *c = NULL;
  Data_Get_Struct(obj, gsl_odeiv_control, c);
  return rb_str_new2(gsl_odeiv_control_name(c));
}

static VALUE rb_gsl_odeiv_control_hadjust(VALUE obj, VALUE ss, VALUE yy0,
					  VALUE yyerr, VALUE ddydt, VALUE hh)
{
  gsl_odeiv_control *c = NULL;
  gsl_odeiv_step *s = NULL;
  gsl_vector *y0 = NULL, *yerr = NULL, *dydt = NULL;
  double h;
  int status;
  CHECK_VECTOR(yy0);
  CHECK_VECTOR(yyerr);
  CHECK_VECTOR(ddydt);
  Data_Get_Struct(obj, gsl_odeiv_control, c);
  Data_Get_Struct(ss, gsl_odeiv_step, s);
  Data_Get_Struct(yy0, gsl_vector, y0);
  Data_Get_Struct(yyerr, gsl_vector, yerr);
  Data_Get_Struct(ddydt, gsl_vector, dydt);
  h = NUM2DBL(hh);
  status = gsl_odeiv_control_hadjust(c, s, y0->data, yerr->data, 
				     dydt->data, &h);
  return rb_ary_new3(2, rb_float_new(h), INT2FIX(status));
}

static gsl_odeiv_evolve* make_evolve(VALUE dim);
static VALUE rb_gsl_odeiv_evolve_new(VALUE klass, VALUE dim) 
{
  gsl_odeiv_evolve *e = NULL;
  e = make_evolve(dim);
  return Data_Wrap_Struct(klass, 0, gsl_odeiv_evolve_free, e);
}

static gsl_odeiv_evolve* make_evolve(VALUE dim)
{
  return gsl_odeiv_evolve_alloc(FIX2INT(dim));
}

static VALUE rb_gsl_odeiv_evolve_reset(VALUE obj)
{
  gsl_odeiv_evolve *e = NULL;
  Data_Get_Struct(obj, gsl_odeiv_evolve, e);
  return INT2FIX(gsl_odeiv_evolve_reset(e));
}

static VALUE rb_gsl_odeiv_evolve_count(VALUE obj)
{
  gsl_odeiv_evolve *e = NULL;
  Data_Get_Struct(obj, gsl_odeiv_evolve, e);
  return INT2FIX(e->count);
}

static VALUE rb_gsl_odeiv_evolve_dimension(VALUE obj)
{
  gsl_odeiv_evolve *e = NULL;
  Data_Get_Struct(obj, gsl_odeiv_evolve, e);
  return INT2FIX(e->dimension);
}

static VALUE rb_gsl_odeiv_evolve_failed_steps(VALUE obj)
{
  gsl_odeiv_evolve *e = NULL;
  Data_Get_Struct(obj, gsl_odeiv_evolve, e);
  return INT2FIX(e->failed_steps);
}

static VALUE rb_gsl_odeiv_evolve_last_step(VALUE obj)
{
  gsl_odeiv_evolve *e = NULL;
  Data_Get_Struct(obj, gsl_odeiv_evolve, e);
  return rb_float_new(e->last_step);
}

static VALUE rb_gsl_odeiv_evolve_y0(VALUE obj)
{
  gsl_odeiv_evolve *e = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_odeiv_evolve, e);
  v = gsl_vector_view_alloc();
  v->vector.data = e->y0;
  v->vector.size = e->dimension;
  v->vector.stride = 1;
  v->vector.owner = 0;
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_odeiv_evolve_yerr(VALUE obj)
{
  gsl_odeiv_evolve *e = NULL;
  gsl_vector_view *v = NULL;
  Data_Get_Struct(obj, gsl_odeiv_evolve, e);
  v = gsl_vector_view_alloc();
  v->vector.data = e->yerr;
  v->vector.size = e->dimension;
  v->vector.stride = 1;
  v->vector.owner = 0;
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, gsl_vector_view_free, v);
}

static VALUE rb_gsl_odeiv_evolve_apply(VALUE obj, VALUE cc, VALUE ss, VALUE sss,
				       VALUE tt, VALUE tt1, VALUE hh, VALUE yy)
{
  gsl_odeiv_evolve *e = NULL;
  gsl_odeiv_control *c = NULL;
  gsl_odeiv_step *s = NULL;
  gsl_odeiv_system *sys = NULL;
  gsl_vector *y = NULL;
  double t, h;
  int status;
  CHECK_STEP(ss); CHECK_SYSTEM(sss);
  CHECK_VECTOR(yy);
  Data_Get_Struct(obj, gsl_odeiv_evolve, e);
  if (NIL_P(cc)) {
    c = NULL;
  } else {
    CHECK_CONTROL(cc); 
    Data_Get_Struct(cc, gsl_odeiv_control, c);
  }
  Data_Get_Struct(ss, gsl_odeiv_step, s);
  Data_Get_Struct(sss, gsl_odeiv_system, sys);
  Data_Get_Struct(yy, gsl_vector, y);
  /*  if (TYPE(tt) != T_FLOAT) rb_raise(rb_eTypeError, "argument 4 Float expected");
      if (TYPE(hh) != T_FLOAT) rb_raise(rb_eTypeError, "argument 6 Float expected");*/
  t = NUM2DBL(tt);
  h = NUM2DBL(hh);
  status = gsl_odeiv_evolve_apply(e, c, s, sys, &t, NUM2DBL(tt1), &h, y->data);
  /*  RFLOAT(tt)->value = t;
      RFLOAT(hh)->value = h;*/
  return rb_ary_new3(3, rb_float_new(t), rb_float_new(h), INT2FIX(status));
}

static void rb_gsl_odeiv_solver_free(gsl_odeiv_solver *gde);
static void gsl_odeiv_solver_mark(gsl_odeiv_solver *gos);
static VALUE rb_gsl_odeiv_solver_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_odeiv_solver *gos = NULL;
  VALUE epsabs, epsrel, ay, adydt;
  VALUE dim;
  if (argc < 4) rb_raise(rb_eArgError, "too few arguments");
  Check_Type(argv[1], T_ARRAY);
  CHECK_PROC(argv[2]);
  if (rb_obj_is_kind_of(argv[3], rb_cProc) || NIL_P(argv[3])) {
    dim = argv[4];
  } else {
    dim = argv[3];
  }
  gos = ALLOC(gsl_odeiv_solver);
  gos->s = make_step(argv[0], dim);
  //  switch (RARRAY(argv[1])->len) {
  switch (RARRAY_LEN(argv[1])) {
  case 2:
    epsabs = rb_ary_entry(argv[1], 0);
    epsrel = rb_ary_entry(argv[1], 1);
    gos->c = make_control_y(epsabs, epsrel);
    break;
  case 4:
    epsabs = rb_ary_entry(argv[1], 0);
    epsrel = rb_ary_entry(argv[1], 1);
    ay = rb_ary_entry(argv[1], 2);
    adydt = rb_ary_entry(argv[1], 3);
    gos->c = make_control_standard(epsabs, epsrel, ay, adydt);
    break;
  default:
    rb_raise(rb_eArgError, "size of the argument 1 must be 2 or 4");
    break;
  }
  gos->sys = make_sys(argc - 2, argv + 2);
  gos->e = make_evolve(dim);
  return Data_Wrap_Struct(klass,  gsl_odeiv_solver_mark, rb_gsl_odeiv_solver_free, gos);
  //  return Data_Wrap_Struct(klass,  0, rb_gsl_odeiv_solver_free, gos);
}

static void gsl_odeiv_solver_mark(gsl_odeiv_solver *gos)
{
  rb_gc_mark((VALUE) gos->sys->params);
}

static VALUE rb_gsl_odeiv_solver_evolve(VALUE obj)
{
  gsl_odeiv_solver *gos = NULL;
  Data_Get_Struct(obj, gsl_odeiv_solver, gos);
  return Data_Wrap_Struct(cgsl_odeiv_evolve, 0, NULL, gos->e);
}

static VALUE rb_gsl_odeiv_solver_set_evolve(VALUE obj, VALUE ee)
{
  gsl_odeiv_solver *gos = NULL;
  gsl_odeiv_evolve *e = NULL;
  CHECK_EVOLVE(ee);
  Data_Get_Struct(obj, gsl_odeiv_solver, gos);
  Data_Get_Struct(ee, gsl_odeiv_evolve, e);
  gos->e = e;
  return obj;
}

static VALUE rb_gsl_odeiv_solver_step(VALUE obj)
{
  gsl_odeiv_solver *gos = NULL;
  Data_Get_Struct(obj, gsl_odeiv_solver, gos);
  return Data_Wrap_Struct(cgsl_odeiv_step, 0, NULL,gos->s);
}

static VALUE rb_gsl_odeiv_solver_set_step(VALUE obj, VALUE ss)
{
  gsl_odeiv_solver *gos = NULL;
  gsl_odeiv_step *s = NULL;
  CHECK_STEP(ss);
  Data_Get_Struct(obj, gsl_odeiv_solver, gos);
  Data_Get_Struct(ss, gsl_odeiv_step, s);
  gos->s = s;
  return obj;
}

static VALUE rb_gsl_odeiv_solver_control(VALUE obj)
{
  gsl_odeiv_solver *gos = NULL;
  Data_Get_Struct(obj, gsl_odeiv_solver, gos);
  return Data_Wrap_Struct(cgsl_odeiv_control, 0, NULL, gos->c);
}

static VALUE rb_gsl_odeiv_solver_set_control(VALUE obj, VALUE cc)
{
  gsl_odeiv_solver *gos = NULL;
  gsl_odeiv_control *c = NULL;
  CHECK_CONTROL(cc);
  Data_Get_Struct(obj, gsl_odeiv_solver, gos);
  Data_Get_Struct(cc, gsl_odeiv_control, c);
  gos->c = c;
  return obj;
}

static VALUE rb_gsl_odeiv_solver_sys(VALUE obj)
{
  gsl_odeiv_solver *gos = NULL;
  Data_Get_Struct(obj, gsl_odeiv_solver, gos);
  return Data_Wrap_Struct(cgsl_odeiv_system, 0, NULL, gos->sys);
}

static VALUE rb_gsl_odeiv_solver_set_sys(VALUE obj, VALUE ss)
{
  gsl_odeiv_solver *gos = NULL;
  gsl_odeiv_system *sys = NULL;
  CHECK_SYSTEM(ss);
  Data_Get_Struct(obj, gsl_odeiv_solver, gos);
  Data_Get_Struct(ss, gsl_odeiv_system, sys);
  gos->sys = sys;
  return obj;
}

static VALUE rb_gsl_odeiv_solver_apply(VALUE obj, VALUE tt, VALUE tt1, VALUE hh, 
				       VALUE yy)
{
  gsl_odeiv_solver *gos = NULL;
  gsl_vector *y = NULL;
  double t, h;
  int status;
  CHECK_VECTOR(yy);
  Need_Float(tt1);
  Data_Get_Struct(obj, gsl_odeiv_solver, gos);
  Data_Get_Struct(yy, gsl_vector, y);
  /*  if (TYPE(tt) != T_FLOAT) rb_raise(rb_eTypeError, "argument 0 Float expected");
      if (TYPE(hh) != T_FLOAT) rb_raise(rb_eTypeError, "argument 2 Float expected");*/
  t = NUM2DBL(tt);
  h = NUM2DBL(hh);
  status = gsl_odeiv_evolve_apply(gos->e, gos->c, gos->s, 
				  gos->sys, &t, NUM2DBL(tt1), &h, y->data);
  /*  RFLOAT(tt)->value = t;
      RFLOAT(hh)->value = h;
      return INT2FIX(status);*/
  return rb_ary_new3(3, rb_float_new(t), rb_float_new(h), INT2FIX(status));
}

static void rb_gsl_odeiv_solver_free(gsl_odeiv_solver *gos)
{
  free((gsl_odeiv_solver *) gos);
}

static VALUE rb_gsl_odeiv_solver_reset(VALUE obj)
{
  gsl_odeiv_solver *gos = NULL;
  Data_Get_Struct(obj, gsl_odeiv_solver, gos);
  gsl_odeiv_step_reset(gos->s);
  gsl_odeiv_evolve_reset(gos->e);
  return obj;
}

static VALUE rb_gsl_odeiv_solver_dim(VALUE obj)
{
  gsl_odeiv_solver *gos = NULL;
  Data_Get_Struct(obj, gsl_odeiv_solver, gos);
  return INT2FIX(gos->sys->dimension);
}

static VALUE rb_gsl_odeiv_solver_set_params(int argc, VALUE *argv, VALUE obj)
{
  gsl_odeiv_solver *gos = NULL;
  Data_Get_Struct(obj, gsl_odeiv_solver, gos);
  rb_gsl_odeiv_system_set_params(argc, argv, 
				 Data_Wrap_Struct(cgsl_odeiv_system, 0, NULL, gos->sys));
  return obj;
}

static VALUE rb_gsl_odeiv_solver_params(VALUE obj)
{
  VALUE ary;
  gsl_odeiv_solver *solver = NULL;
  Data_Get_Struct(obj, gsl_odeiv_solver, solver);
  ary = (VALUE) solver->sys->params;
  return rb_ary_entry(ary, 3);
}

void Init_gsl_odeiv(VALUE module)
{
  VALUE mgsl_odeiv;
  mgsl_odeiv = rb_define_module_under(module, "Odeiv");
  rb_define_const(mgsl_odeiv, "HADJ_DEC", INT2FIX(GSL_ODEIV_HADJ_DEC));
  rb_define_const(mgsl_odeiv, "HADJ_INC", INT2FIX(GSL_ODEIV_HADJ_INC));
  rb_define_const(mgsl_odeiv, "HADJ_NIL", INT2FIX(GSL_ODEIV_HADJ_NIL));
  cgsl_odeiv_step = rb_define_class_under(mgsl_odeiv, "Step", 
					    cGSL_Object);
  rb_define_singleton_method(cgsl_odeiv_step, "alloc", rb_gsl_odeiv_step_new, -1);

  /*****/
  rb_define_const(cgsl_odeiv_step, "RK2", INT2FIX(GSL_ODEIV_STEP_RK2));
  rb_define_const(cgsl_odeiv_step, "RK4", INT2FIX(GSL_ODEIV_STEP_RK4));
  rb_define_const(cgsl_odeiv_step, "RKF45", INT2FIX(GSL_ODEIV_STEP_RKF45));
  rb_define_const(cgsl_odeiv_step, "RKCK", INT2FIX(GSL_ODEIV_STEP_RKCK));
  rb_define_const(cgsl_odeiv_step, "RK8PD", INT2FIX(GSL_ODEIV_STEP_RK8PD));
  rb_define_const(cgsl_odeiv_step, "RK2IMP", INT2FIX(GSL_ODEIV_STEP_RK2IMP));
  rb_define_const(cgsl_odeiv_step, "RK4IMP", INT2FIX(GSL_ODEIV_STEP_RK4IMP));
  rb_define_const(cgsl_odeiv_step, "BSIMP", INT2FIX(GSL_ODEIV_STEP_BSIMP));
  rb_define_const(cgsl_odeiv_step, "GEAR1", INT2FIX(GSL_ODEIV_STEP_GEAR1));
  rb_define_const(cgsl_odeiv_step, "GEAR2", INT2FIX(GSL_ODEIV_STEP_GEAR2));
  rb_define_const(cgsl_odeiv_step, "RK2SIMP", INT2FIX(GSL_ODEIV_STEP_RK2SIMP));
  /*****/

  rb_define_method(cgsl_odeiv_step, "reset", rb_gsl_odeiv_step_reset, 0);
  rb_define_method(cgsl_odeiv_step, "name", rb_gsl_odeiv_step_name, 0);
  rb_define_method(cgsl_odeiv_step, "order", rb_gsl_odeiv_step_order, 0);
  rb_define_method(cgsl_odeiv_step, "dimension", rb_gsl_odeiv_step_dimension, 0);
  rb_define_alias(cgsl_odeiv_step, "dim", "dimension");
  rb_define_method(cgsl_odeiv_step, "apply", rb_gsl_odeiv_step_apply, -1);
  rb_define_method(cgsl_odeiv_step, "info", rb_gsl_odeiv_step_info, 0);

  /****/

  cgsl_odeiv_control = rb_define_class_under(mgsl_odeiv, "Control", 
					    cGSL_Object);
  rb_define_singleton_method(cgsl_odeiv_control, "alloc", rb_gsl_odeiv_control_standard_new, 4);
  rb_define_singleton_method(cgsl_odeiv_control, "standard_alloc", rb_gsl_odeiv_control_standard_new, 4);
 rb_define_singleton_method(cgsl_odeiv_control, "y_new", rb_gsl_odeiv_control_y_new, 2);
 rb_define_singleton_method(cgsl_odeiv_control, "yp_new", rb_gsl_odeiv_control_yp_new, 2);
#ifdef GSL_1_2_LATER
 rb_define_singleton_method(cgsl_odeiv_control, "scaled_alloc", rb_gsl_odeiv_control_scaled_new, 5);
#endif

  rb_define_method(cgsl_odeiv_control, "init", rb_gsl_odeiv_control_init, 4);
  rb_define_method(cgsl_odeiv_control, "name", rb_gsl_odeiv_control_name, 0);
  rb_define_method(cgsl_odeiv_control, "hadjust", rb_gsl_odeiv_control_hadjust, 5);

 /****/
  cgsl_odeiv_evolve = rb_define_class_under(mgsl_odeiv, "Evolve", 
					     cGSL_Object);
  rb_define_singleton_method(cgsl_odeiv_evolve, "alloc", rb_gsl_odeiv_evolve_new, 1);

  rb_define_method(cgsl_odeiv_evolve, "reset", rb_gsl_odeiv_evolve_reset, 0);
  rb_define_method(cgsl_odeiv_evolve, "apply", rb_gsl_odeiv_evolve_apply, 7);
  rb_define_method(cgsl_odeiv_evolve, "count", rb_gsl_odeiv_evolve_count, 0);
  rb_define_method(cgsl_odeiv_evolve, "dimension", rb_gsl_odeiv_evolve_dimension, 0);
  rb_define_method(cgsl_odeiv_evolve, "failed_steps", rb_gsl_odeiv_evolve_failed_steps, 0);
  rb_define_method(cgsl_odeiv_evolve, "last_step", rb_gsl_odeiv_evolve_last_step, 0);
  rb_define_method(cgsl_odeiv_evolve, "y0", rb_gsl_odeiv_evolve_y0, 0);
  rb_define_method(cgsl_odeiv_evolve, "yerr", rb_gsl_odeiv_evolve_yerr, 0);
  /*****/

  cgsl_odeiv_system = rb_define_class_under(mgsl_odeiv, "System", 
					    cGSL_Object);
  rb_define_singleton_method(cgsl_odeiv_system, "alloc", rb_gsl_odeiv_system_new, -1);
  rb_define_method(cgsl_odeiv_system, "set", rb_gsl_odeiv_system_set, -1);
  rb_define_method(cgsl_odeiv_system, "set_params", 
			     rb_gsl_odeiv_system_set_params, -1);
  rb_define_method(cgsl_odeiv_system, "params", 
			     rb_gsl_odeiv_system_params, 0);
  rb_define_method(cgsl_odeiv_system, "function", 
			     rb_gsl_odeiv_system_function, 0);
  rb_define_alias(cgsl_odeiv_system, "func", "function");
  rb_define_method(cgsl_odeiv_system, "jacobian", 
			     rb_gsl_odeiv_system_jacobian, 0);
  rb_define_alias(cgsl_odeiv_system, "jac", "jacobian");
  rb_define_method(cgsl_odeiv_system, "dimension", rb_gsl_odeiv_system_dimension, 0);
  rb_define_alias(cgsl_odeiv_system, "dim", "dimension");

  /*****/
  cgsl_odeiv_solver = rb_define_class_under(mgsl_odeiv, "Solver", cGSL_Object);
  rb_define_singleton_method(cgsl_odeiv_solver, "alloc", rb_gsl_odeiv_solver_new, -1);
  rb_define_method(cgsl_odeiv_solver, "step", rb_gsl_odeiv_solver_step, 0);
  rb_define_method(cgsl_odeiv_solver, "control", rb_gsl_odeiv_solver_control, 0);
  rb_define_method(cgsl_odeiv_solver, "evolve", rb_gsl_odeiv_solver_evolve, 0);
  rb_define_method(cgsl_odeiv_solver, "sys", rb_gsl_odeiv_solver_sys, 0);
  rb_define_method(cgsl_odeiv_solver, "apply", rb_gsl_odeiv_solver_apply, 4);

  rb_define_method(cgsl_odeiv_solver, "set_evolve", rb_gsl_odeiv_solver_set_evolve, 1);
  rb_define_method(cgsl_odeiv_solver, "set_step", rb_gsl_odeiv_solver_set_step, 1);
  rb_define_method(cgsl_odeiv_solver, "set_control", rb_gsl_odeiv_solver_set_control, 1);
  rb_define_method(cgsl_odeiv_solver, "set_system", rb_gsl_odeiv_solver_set_sys, 1);
  rb_define_method(cgsl_odeiv_solver, "reset", rb_gsl_odeiv_solver_reset, 0);
  rb_define_method(cgsl_odeiv_solver, "dim", rb_gsl_odeiv_solver_dim, 0);
  rb_define_alias(cgsl_odeiv_solver, "dimension", "dim");
  rb_define_method(cgsl_odeiv_solver, "set_params", rb_gsl_odeiv_solver_set_params, -1);
  rb_define_method(cgsl_odeiv_solver, "params", rb_gsl_odeiv_solver_params, 0);
}

#undef CHECK_SOLVER
#undef CHECK_EVOLVE
#undef CHECK_CONTROL
#undef CHECK_STEP
#undef CHECK_SYSTEM
