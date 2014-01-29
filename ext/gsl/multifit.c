/*
  multifit.c
  Ruby/GSL: Ruby extension library for GSL (GNU Scientific Library)
    (C) Copyright 2004 by Yoshiki Tsunesada

  Ruby/GSL is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License.
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY.
*/

#include "include/rb_gsl_fit.h"
#include "include/rb_gsl_function.h"
#include "include/rb_gsl_common.h"

#ifdef HAVE_NDLINEAR_GSL_MULTIFIT_NDLINEAR_H
#include <ndlinear/gsl_multifit_ndlinear.h>
void Init_ndlinear(VALUE module);
#endif

#ifndef CHECK_WORKSPACE
#define CHECK_WORKSPACE(x) if(CLASS_OF(x)!=cgsl_multifit_workspace)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (GSL::MultiFit::Workspace expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

#ifndef CHECK_MULTIFIT_FUNCTION_FDF
#define CHECK_MULTIFIT_FUNCTION_FDF(x) if(CLASS_OF(x)!=cgsl_multifit_function_fdf)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (GSL::MultiFit::Workspace expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

static VALUE cgsl_multifit_workspace;
static VALUE cgsl_multifit_function_fdf;

enum MultiFit_Solver_Type{
  GSL_MULTIFIT_FDFSOLVER_LMSDER,
  GSL_MULTIFIT_FDFSOLVER_LMDER,
};

static VALUE rb_gsl_multifit_workspace_new(VALUE klass, VALUE n, VALUE p)
{
  gsl_multifit_linear_workspace *w = NULL;
  CHECK_FIXNUM(n); CHECK_FIXNUM(p);
  w = gsl_multifit_linear_alloc(FIX2INT(n), FIX2INT(p));
  return Data_Wrap_Struct(cgsl_multifit_workspace, 0, gsl_multifit_linear_free, w);
}

static VALUE rb_gsl_multifit_linear(int argc, VALUE *argv, VALUE obj)
{
  gsl_multifit_linear_workspace *space = NULL;
  gsl_matrix *x = NULL, *cov = NULL;
  gsl_vector *y = NULL, *c = NULL;
  double chisq;
  int status, flag = 0;
  VALUE vc, vcov;
  if (argc != 2 && argc != 3) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
  Data_Get_Matrix(argv[0], x);
  Data_Get_Vector(argv[1], y);
  if (argc == 3) {
    CHECK_WORKSPACE(argv[2]);
    Data_Get_Struct(argv[2], gsl_multifit_linear_workspace, space);
  } else {
    space = gsl_multifit_linear_alloc(x->size1, x->size2);
    flag = 1;
  }
  cov = gsl_matrix_alloc(x->size2, x->size2);
  c = gsl_vector_alloc(x->size2);
  status = gsl_multifit_linear(x, y, c, cov, &chisq, space);
  if (flag == 1) gsl_multifit_linear_free(space);
  vc = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, c);
  vcov = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, cov);
  return rb_ary_new3(4, vc, vcov, rb_float_new(chisq), INT2FIX(status));
}

static VALUE rb_gsl_multifit_wlinear(int argc, VALUE *argv, VALUE obj)
{
  gsl_multifit_linear_workspace *space = NULL;
  gsl_matrix *x = NULL, *cov = NULL;
  gsl_vector *w = NULL, *y = NULL, *c = NULL;
  double chisq;
  int status, flag = 0;
  VALUE vc, vcov;
  if (argc != 3 && argc != 4) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 3 or 4)", argc);
  Data_Get_Matrix(argv[0], x);
  Data_Get_Vector(argv[1], w);
  Data_Get_Vector(argv[2], y);
  if (argc == 4) {
    CHECK_WORKSPACE(argv[3]);
    Data_Get_Struct(argv[3], gsl_multifit_linear_workspace, space);
  } else {
    space = gsl_multifit_linear_alloc(x->size1, x->size2);
    flag = 1;
  }
  cov = gsl_matrix_alloc(x->size2, x->size2);
  c = gsl_vector_alloc(x->size2);
  status = gsl_multifit_wlinear(x, w, y, c, cov, &chisq, space);
  if (flag == 1) gsl_multifit_linear_free(space);
  vc = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, c);
  vcov = Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, cov);
  return rb_ary_new3(4, vc, vcov, rb_float_new(chisq), INT2FIX(status));
}

static void calc_X_power(gsl_matrix *X, gsl_vector *x, size_t order)
{
  size_t i, j;
  double val;
  for (i = 0; i < x->size; i++) {
    val = 1.0;
    gsl_matrix_set(X, i, 0, val);
    for (j = 1; j <= order; j++) {
      val *= gsl_vector_get(x, i);
      gsl_matrix_set(X, i, j, val);
    }
  }
}

static void calc_X_legendre(gsl_matrix *X, gsl_vector *v, size_t order)
{
  size_t i, n;
  double a, b, c, x;
  for (i = 0; i < v->size; i++) {
    a = 1.0;
    b = gsl_vector_get(v, i);
    gsl_matrix_set(X, i, 0, a);
    gsl_matrix_set(X, i, 1, b);
    for (n = 2; n <= order; n++) {
      x =  gsl_vector_get(v, i);
      c = ((2*n-1)*x*b - (n-1)*a)/n;
      gsl_matrix_set(X, i, n, c);
      a = b;
      b = c;
    }
  }
}

static VALUE rb_gsl_multifit_XXXfit(int argc, VALUE *argv, VALUE obj,
				    void (*fn)(gsl_matrix*, gsl_vector*,size_t))
{
  gsl_multifit_linear_workspace *space = NULL;
  gsl_matrix *X = NULL, *cov = NULL;
  gsl_vector *x, *y = NULL, *w, *c = NULL, *err;
  gsl_vector_view xx, yy, ww;
  size_t order, i;
  double chisq;
  int status, flag = 0, flagw = 0;
  VALUE vc, verr;
  if (argc != 3 && argc != 4 && argc != 5) 
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 3 or 4)", argc);
  x = &xx.vector;
  w = &ww.vector;
  y = &yy.vector;
  Data_Get_Vector(argv[0], x);
  if (argc >= 3 && VECTOR_P(argv[2])) {
    Data_Get_Vector(argv[1], w);
    Data_Get_Vector(argv[2], y);
    order = NUM2INT(argv[3]);
    flagw = 1;
  } else {
    Data_Get_Vector(argv[1], y);
    order = NUM2INT(argv[2]);
    flagw = 0;
  }
  if (rb_obj_is_kind_of(argv[argc-1], cgsl_multifit_workspace)) {
    Data_Get_Struct(argv[argc-1], gsl_multifit_linear_workspace, space);
  } else {
    space = gsl_multifit_linear_alloc(x->size, order + 1);
    flag = 1;
  }
  cov = gsl_matrix_alloc(order + 1, order + 1);
  c = gsl_vector_alloc(order + 1);
  X = gsl_matrix_alloc(x->size, order + 1);
  (*fn)(X, x, order);
  if (flagw == 0) status = gsl_multifit_linear(X, y, c, cov, &chisq, space);
  else   status = gsl_multifit_wlinear(X, w, y, c, cov, &chisq, space);
  if (flag == 1) gsl_multifit_linear_free(space);
  err = gsl_vector_alloc(order + 1);
  vc = Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, c);
  verr = Data_Wrap_Struct(cgsl_poly, 0, gsl_vector_free, err);
  for (i = 0; i < err->size; i++) 
    gsl_vector_set(err, i, sqrt(chisq/((double)x->size-err->size)*gsl_matrix_get(cov, i, i)));
  
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  return rb_ary_new3(4, vc, verr, rb_float_new(chisq), INT2FIX(status));
}

static VALUE rb_gsl_multifit_polyfit(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_multifit_XXXfit(argc, argv, obj, calc_X_power);
}

static VALUE rb_gsl_multifit_legfit(int argc, VALUE *argv, VALUE obj)
{
  return rb_gsl_multifit_XXXfit(argc, argv, obj, calc_X_legendre);
}

/**********/

static VALUE rb_gsl_multifit_fdfsolver_new(int argc, VALUE *argv, VALUE klass)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *solver = NULL;
  size_t n, p;
  char name[64];
  switch (argc) {
  case 3:
    switch (TYPE(argv[0])) {
    case T_STRING:
      strcpy(name, STR2CSTR(argv[0]));
      if (str_tail_grep(name, "lmsder") == 0) {
	T = gsl_multifit_fdfsolver_lmsder;
      } else if (str_tail_grep(name, "lmder") == 0) {
	T = gsl_multifit_fdfsolver_lmder;
      } else {
	rb_raise(rb_eTypeError, "unknown solver type %s (lmsder of lmder)",
		 name);
      }
      break;
    case T_FIXNUM:
      switch (FIX2INT(argv[0])) {
      case GSL_MULTIFIT_FDFSOLVER_LMSDER:
	T = gsl_multifit_fdfsolver_lmsder;
	break;
      case GSL_MULTIFIT_FDFSOLVER_LMDER:
	T = gsl_multifit_fdfsolver_lmder;
	break;
      default:
	rb_raise(rb_eTypeError, 
		 "unknown solver type (GSL::MultiFit::FdfSolver::LMSDER or LMDER expected)");
	break;
      }
      break;
    default:
      rb_raise(rb_eTypeError, "wrong argument type %s (Fixnum or String)",
	       rb_class2name(CLASS_OF(argv[0])));
      break;
    }
    CHECK_FIXNUM(argv[1]); CHECK_FIXNUM(argv[2]);
    n = FIX2INT(argv[1]);
    p = FIX2INT(argv[2]);
    break;
  case 2:
    CHECK_FIXNUM(argv[0]); CHECK_FIXNUM(argv[1]);
    T = gsl_multifit_fdfsolver_lmsder;
    n = FIX2INT(argv[0]);
    p = FIX2INT(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  solver = gsl_multifit_fdfsolver_alloc(T, n, p);
  return Data_Wrap_Struct(klass, 0, gsl_multifit_fdfsolver_free, solver);
}

static VALUE rb_gsl_multifit_fdfsolver_set(VALUE obj, VALUE ff, VALUE xx)
{
  gsl_multifit_fdfsolver *solver = NULL;
  gsl_multifit_function_fdf *f = NULL;
  gsl_vector *x = NULL;
  int status;
  CHECK_MULTIFIT_FUNCTION_FDF(ff);
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  Data_Get_Struct(ff, gsl_multifit_function_fdf, f);
  Data_Get_Vector(xx, x);
  status = gsl_multifit_fdfsolver_set(solver, f, x);
  return INT2FIX(status);
}

static VALUE rb_gsl_multifit_fdfsolver_name(VALUE obj)
{
  gsl_multifit_fdfsolver *solver = NULL;
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  return rb_str_new2(gsl_multifit_fdfsolver_name(solver));
}

static VALUE rb_gsl_multifit_fdfsolver_iterate(VALUE obj)
{
  gsl_multifit_fdfsolver *solver = NULL;
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  return INT2FIX(gsl_multifit_fdfsolver_iterate(solver));
}

static VALUE rb_gsl_multifit_fdfsolver_position(VALUE obj)
{
  gsl_multifit_fdfsolver *solver = NULL;
  gsl_vector *x = NULL;
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  x = gsl_multifit_fdfsolver_position(solver);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, x);
}

static VALUE rb_gsl_multifit_fdfsolver_print_state(VALUE obj, VALUE i)
{
  gsl_multifit_fdfsolver *solver = NULL;
  CHECK_FIXNUM(i);
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  printf("iter: %d x = %15.8f %15.8f %15.8f |f(x)| = %g\n",
	 (int) FIX2INT(i), gsl_vector_get(solver->x, 0), gsl_vector_get(solver->x, 1),
	 gsl_vector_get(solver->x, 2), gsl_blas_dnrm2(solver->f));
  return Qtrue;
}

static VALUE rb_gsl_multifit_fdfsolver_fdf(VALUE obj)
{
  gsl_multifit_fdfsolver *solver = NULL;
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  return Data_Wrap_Struct(cgsl_multifit_function_fdf, 0, NULL, solver->fdf);
}

static VALUE rb_gsl_multifit_fdfsolver_test_delta(VALUE obj, VALUE r, VALUE a)
{
  gsl_multifit_fdfsolver *solver = NULL;
  Need_Float(r); Need_Float(a);
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  return INT2FIX(gsl_multifit_test_delta(solver->dx, solver->x, NUM2DBL(r), NUM2DBL(a)));
}

static VALUE rb_gsl_multifit_fdfsolver_test_gradient(int argc, VALUE *argv, VALUE obj)
{
  gsl_multifit_fdfsolver *solver = NULL;
  gsl_vector *g = NULL;
  int status;
  double epsabs;
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  switch (argc) {
  case 1:
    Need_Float(argv[0]);
    g = gsl_vector_alloc(solver->x->size);
    gsl_multifit_gradient(solver->J, solver->f, g);
    epsabs = NUM2DBL(argv[0]);
    status = gsl_multifit_test_gradient(g, epsabs);
    gsl_vector_free(g);
    break;
  case 2:
    Need_Float(argv[1]);
    Data_Get_Vector(argv[0], g);
    epsabs = NUM2DBL(argv[1]);
    status = gsl_multifit_test_gradient(g, epsabs);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
  return INT2FIX(status);
}

static VALUE rb_gsl_multifit_fdfsolver_gradient(int argc, VALUE *argv, VALUE obj)
{
  gsl_multifit_fdfsolver *solver = NULL;
  gsl_vector *g = NULL;
  // local variable "status" declared and set, but never used
  //int status;
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  if (argc == 1) {
    Data_Get_Vector(argv[0], g);
    return INT2FIX(gsl_multifit_gradient(solver->J, solver->f, g));
  } else {
    g = gsl_vector_alloc(solver->x->size);
    /*status =*/ gsl_multifit_gradient(solver->J, solver->f, g);
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, g);
  }
}

static VALUE rb_gsl_multifit_fdfsolver_covar(int argc, VALUE *argv, VALUE obj)
{
  gsl_multifit_fdfsolver *solver = NULL;
  gsl_matrix *covar = NULL;
  double epsrel;
  // local variable "status" declared and set, but never used
  //int status;
  if (argc < 1) rb_raise(rb_eArgError, "too few arguments");
  Need_Float(argv[0]);
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  epsrel = NUM2DBL(argv[0]);
  switch (argc) {
  case 1:
    covar = gsl_matrix_alloc(solver->x->size, solver->x->size);
    /*status =*/ gsl_multifit_covar(solver->J, epsrel, covar);
    return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, covar);
    break;
  case 2:
    Data_Get_Matrix(argv[1], covar);
    return INT2FIX(gsl_multifit_covar(solver->J, epsrel, covar));
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
}

static VALUE rb_gsl_multifit_fdfsolver_x(VALUE obj)
{
  gsl_multifit_fdfsolver *solver = NULL;
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, solver->x);
}

static VALUE rb_gsl_multifit_fdfsolver_dx(VALUE obj)
{
  gsl_multifit_fdfsolver *solver = NULL;
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, solver->dx);
}

static VALUE rb_gsl_multifit_fdfsolver_f(VALUE obj)
{
  gsl_multifit_fdfsolver *solver = NULL;
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, solver->f);
}

static VALUE rb_gsl_multifit_fdfsolver_J(VALUE obj)
{
  gsl_multifit_fdfsolver *solver = NULL;
  Data_Get_Struct(obj, gsl_multifit_fdfsolver, solver);
  return Data_Wrap_Struct(cgsl_matrix_view_ro, 0, NULL, solver->J);
}

/* singleton */
static VALUE rb_gsl_multifit_test_delta(VALUE obj, VALUE ddx, VALUE xx, VALUE a, VALUE r)
{
  gsl_vector *dx = NULL, *x = NULL;
  Need_Float(a); Need_Float(r);
  Data_Get_Vector(ddx, dx);
  Data_Get_Vector(xx, x);
  return INT2FIX(gsl_multifit_test_delta(dx, x, NUM2DBL(a), NUM2DBL(r)));
}

static VALUE rb_gsl_multifit_test_gradient(VALUE obj, VALUE gg, VALUE a)
{
  gsl_vector *g = NULL;
  Need_Float(a);
  Data_Get_Vector(gg, g);
  return INT2FIX(gsl_multifit_test_gradient(g, NUM2DBL(a)));
}

static VALUE rb_gsl_multifit_gradient(int argc, VALUE *argv, VALUE obj)
{
  gsl_vector *g = NULL, *f = NULL;
  gsl_matrix *J = NULL;
  // local variable "status" declared and set, but never used
  //int status;
  switch (argc) {
  case 2:
    Data_Get_Matrix(argv[0], J);
    Data_Get_Vector(argv[1], f);  
    g = gsl_vector_alloc(f->size);
    /*status =*/ gsl_multifit_gradient(J, f, g);
    return Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, g);
    break;  
  case 3:
    Data_Get_Matrix(argv[0], J);
    Data_Get_Vector(argv[1], f);  
    Data_Get_Vector(argv[2], g);  
    return INT2FIX(gsl_multifit_gradient(J, f, g));
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
}

static VALUE rb_gsl_multifit_covar(int argc, VALUE *argv, VALUE obj)
{
  gsl_matrix *J = NULL, *covar = NULL;
  double epsrel;
  int status;
  switch (argc) {
  case 2:
    Need_Float(argv[1]);
    Data_Get_Matrix(argv[0], J);
    epsrel = NUM2DBL(argv[1]);
    covar = gsl_matrix_alloc(J->size2, J->size2);
    /*status =*/ gsl_multifit_covar(J, epsrel, covar);
    return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, covar);
    break;
  case 3:
    Need_Float(argv[1]);
    Data_Get_Matrix(argv[0], J);
    epsrel = NUM2DBL(argv[1]);
    Data_Get_Matrix(argv[2], covar);
    status = gsl_multifit_covar(J, epsrel, covar);
    return INT2FIX(status);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }
}

static void rb_gsl_multifit_define_const(VALUE klass)
{
  rb_define_const(klass, "LMSDER", 
		  INT2FIX(GSL_MULTIFIT_FDFSOLVER_LMSDER));
  rb_define_const(klass, "LMDER", 
		  INT2FIX(GSL_MULTIFIT_FDFSOLVER_LMDER));

  rb_define_const(klass, "Lmsder", 
		  INT2FIX(GSL_MULTIFIT_FDFSOLVER_LMSDER));
  rb_define_const(klass, "Lmder", 
		  INT2FIX(GSL_MULTIFIT_FDFSOLVER_LMDER));
}

static void gsl_multifit_function_fdf_free(gsl_multifit_function_fdf *f);
static void gsl_multifit_function_fdf_free(gsl_multifit_function_fdf *f)
{
  free((gsl_multifit_function_fdf *) f);
}

static void gsl_multifit_function_fdf_mark(gsl_multifit_function_fdf *F);
static void gsl_multifit_function_fdf_mark(gsl_multifit_function_fdf *F)
{
  rb_gc_mark((VALUE) F->params);
}

static int gsl_multifit_function_fdf_f(const gsl_vector *x, void *params,
				       gsl_vector *f);
static int gsl_multifit_function_fdf_df(const gsl_vector *x, void *params,
					gsl_matrix *J);
static int gsl_multifit_function_fdf_fdf(const gsl_vector *x, void *params,
					 gsl_vector *f, gsl_matrix *J);

static VALUE rb_gsl_multifit_function_fdf_set_procs(int argc, VALUE *argv, VALUE obj);

static VALUE rb_gsl_multifit_function_fdf_new(int argc, VALUE *argv, VALUE klass)
{
  gsl_multifit_function_fdf *func = NULL;
  VALUE obj;
  func = ALLOC(gsl_multifit_function_fdf);
  func->f = &gsl_multifit_function_fdf_f;
  func->df = &gsl_multifit_function_fdf_df;
  func->fdf = &gsl_multifit_function_fdf_fdf;
  func->params = NULL;
  obj = Data_Wrap_Struct(klass, gsl_multifit_function_fdf_mark, gsl_multifit_function_fdf_free, func);
  switch (argc) {
  case 0:
  case 1:
    break;
  case 2:
  case 3:
    rb_gsl_multifit_function_fdf_set_procs(argc, argv, obj);
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 0-3)", argc);
    break;
  }
  return obj;
}			    

static VALUE rb_gsl_multifit_function_fdf_set_procs(int argc, VALUE *argv, VALUE obj)
{
  gsl_multifit_function_fdf *func = NULL;
  VALUE ary;
  Data_Get_Struct(obj, gsl_multifit_function_fdf, func);
  if (func->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) func->params = ary;*/
    func->params = (void *) ary;
  } else {
    ary = (VALUE) func->params;
  }
  rb_ary_store(ary, 0, argv[0]);
  rb_ary_store(ary, 1, argv[1]);
  switch (argc) {
	case 2:
		break;
  case 3:
    if (TYPE(argv[2]) == T_FIXNUM) {
      func->p = FIX2INT(argv[2]);
      rb_ary_store(ary, 2, Qnil);
    } else rb_ary_store(ary, 2, argv[2]);
    break;
  case 4:
    if (TYPE(argv[2]) == T_FIXNUM) {
      func->p = FIX2INT(argv[2]);
      rb_ary_store(ary, 2, argv[3]);
    } else {
      func->p = FIX2INT(argv[3]);
      rb_ary_store(ary, 2, argv[2]);
    }
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 3 or 4)", argc);
    break;
  }
  return obj;
}

static VALUE rb_gsl_multifit_function_fdf_set_data(int argc, VALUE *argv, VALUE obj)
{
  VALUE ary, ary2;
  gsl_multifit_function_fdf *func = NULL;
  Data_Get_Struct(obj, gsl_multifit_function_fdf, func);
  if (func->params == NULL) {
    ary = rb_ary_new2(4);
    /*    (VALUE) func->params = ary;*/
    func->params = (void *) ary;
  } else {
    ary = (VALUE) func->params;
  }
  switch (argc) {
  case 2:
    ary2 = rb_ary_new3(2, argv[0], argv[1]);  /* t, y */
    break;
  case 3:
    ary2 = rb_ary_new3(3, argv[0], argv[1], argv[2]);  /* t, y, sigma */
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments (%d for 2 or 3)", argc);
    break;
  }
  func->n = NUM2INT(rb_funcall(argv[0], rb_intern("size"), 0));
  rb_ary_store(ary, 3, ary2);
  return obj;
}

static int gsl_multifit_function_fdf_f(const gsl_vector *x, void *params,
				       gsl_vector *f)
{
  VALUE vt_y_sigma, vt, vy, vsigma, vf, vx, proc, ary;
  ary = (VALUE) params;
  vt_y_sigma = rb_ary_entry(ary, 3);
  proc = rb_ary_entry(ary, 0);
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  vf = Data_Wrap_Struct(cgsl_vector, 0, NULL, f);
  switch (RARRAY_LEN(vt_y_sigma)) {
  case 2:
    vt = rb_ary_entry(vt_y_sigma, 0);
    vy = rb_ary_entry(vt_y_sigma, 1);
    rb_funcall(proc, RBGSL_ID_call, 4, vx, vt, vy, vf);
    break;
  case 3:
    vt = rb_ary_entry(vt_y_sigma, 0);
    vy = rb_ary_entry(vt_y_sigma, 1);
    vsigma = rb_ary_entry(vt_y_sigma,2);
    rb_funcall(proc, RBGSL_ID_call, 5, vx, vt, vy, vsigma, vf);
    break;
  default:
    rb_raise(rb_eArgError, "bad argument");
    break;    
  }
  return GSL_SUCCESS;
}

static int gsl_multifit_function_fdf_df(const gsl_vector *x, void *params,
					gsl_matrix *J)
{
  VALUE vt_y_sigma, vt, vy, vsigma, vJ, vx, proc, ary;
  ary = (VALUE) params;
  vt_y_sigma = rb_ary_entry(ary, 3);
  proc = rb_ary_entry(ary, 1);
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  vJ = Data_Wrap_Struct(cgsl_matrix, 0, NULL, J);
  switch (RARRAY_LEN(vt_y_sigma)) {
  case 2:
    vt = rb_ary_entry(vt_y_sigma, 0);
    vy = rb_ary_entry(vt_y_sigma, 1);
    rb_funcall(proc, RBGSL_ID_call, 4, vx, vt, vy, vJ);
    break;
  case 3:
    vt = rb_ary_entry(vt_y_sigma, 0);
    vy = rb_ary_entry(vt_y_sigma, 1);
    vsigma = rb_ary_entry(vt_y_sigma,2);
    rb_funcall(proc, RBGSL_ID_call, 5, vx, vt, vy, vsigma, vJ);
    break;
  default:
    rb_raise(rb_eArgError, "bad argument");
    break;    
  }
  return GSL_SUCCESS;
}

static int gsl_multifit_function_fdf_fdf(const gsl_vector *x, void *params,
					 gsl_vector *f, gsl_matrix *J)
{
  VALUE vt_y_sigma, vt, vy, vsigma, vf, vJ, vx, proc_f, proc_df, proc_fdf;
  VALUE ary;
  ary = (VALUE) params;
  vt_y_sigma = rb_ary_entry(ary, 3);
  proc_f = rb_ary_entry(ary, 0);
  proc_df = rb_ary_entry(ary, 1);
  proc_fdf = rb_ary_entry(ary, 2);
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  vf = Data_Wrap_Struct(cgsl_vector, 0, NULL, f);
  vJ = Data_Wrap_Struct(cgsl_matrix, 0, NULL, J);
  switch (RARRAY_LEN(vt_y_sigma)) {
  case 2:
    vt = rb_ary_entry(vt_y_sigma, 0);
    vy = rb_ary_entry(vt_y_sigma, 1);
    if (NIL_P(proc_fdf)) {
      rb_funcall(proc_f, RBGSL_ID_call, 4, vx, vt, vy, vf);
      rb_funcall(proc_df, RBGSL_ID_call, 4, vx, vt, vy, vJ);
    } else 
      rb_funcall(proc_fdf, RBGSL_ID_call, 5, vx, vt, vy, vf, vJ);
    break;
  case 3:
    vt = rb_ary_entry(vt_y_sigma, 0);
    vy = rb_ary_entry(vt_y_sigma, 1);
    vsigma = rb_ary_entry(vt_y_sigma,2);
    if (NIL_P(proc_fdf)) {
      rb_funcall(proc_f, RBGSL_ID_call, 5, vx, vt, vy, vsigma, vf);
      rb_funcall(proc_df, RBGSL_ID_call, 5, vx, vt, vy, vsigma, vJ);
    } else 
      rb_funcall(proc_fdf, RBGSL_ID_call, 6, vx, vt, vy, vsigma, vf, vJ);
    break;
  default:
    rb_raise(rb_eArgError, "bad argument");
    break;    
  }
  return GSL_SUCCESS;
}

static VALUE rb_gsl_multifit_function_fdf_params(VALUE obj)
{
  gsl_multifit_function_fdf *f = NULL;
  Data_Get_Struct(obj, gsl_multifit_function_fdf, f);
  return (VALUE) f->params;
}

static VALUE rb_gsl_multifit_function_fdf_n(VALUE obj)
{
  gsl_multifit_function_fdf *f = NULL;
  Data_Get_Struct(obj, gsl_multifit_function_fdf, f);
  return INT2FIX(f->n);
}

static VALUE rb_gsl_multifit_function_fdf_set_n(VALUE obj, VALUE n)
{
  gsl_multifit_function_fdf *f = NULL;
  Data_Get_Struct(obj, gsl_multifit_function_fdf, f);
  f->n = FIX2INT(n);
  return obj;
}

static VALUE rb_gsl_multifit_function_fdf_p(VALUE obj)
{
  gsl_multifit_function_fdf *f = NULL;
  Data_Get_Struct(obj, gsl_multifit_function_fdf, f);
  return INT2FIX(f->p);
}

/*****/
struct fitting_xydata {
  gsl_vector *x, *y, *w;
};

/* Gaussian fit 
   y = y0 + A exp(-(x-x0/sigma)^2)
       v[0] = y0
       v[1] = A
       v[2] = x0
       v[3] = var = sigma^2
*/
static int Gaussian_f(const gsl_vector *v, void *data, gsl_vector *f)
{
  gsl_vector *x, *y, *w;
  double y0, A, x0, var, xi, yi, wi;
  size_t i;
  struct fitting_xydata *xydata;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  y = xydata->y;
  w = xydata->w;
  var = gsl_vector_get(v, 3);
  x0 = gsl_vector_get(v, 2);
  A = gsl_vector_get(v, 1);
  y0 = gsl_vector_get(v, 0);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    yi = gsl_vector_get(y, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    gsl_vector_set(f, i, (A*exp(-(xi - x0)*(xi - x0)/var/2.0) + y0 - yi)*wi);
  }
  return GSL_SUCCESS;
}

static int Gaussian_df(const gsl_vector *v, void *data, gsl_matrix *J)
{
  double A, x0, var, xi, yy, wi;
  size_t i;
  struct fitting_xydata *xydata;
  gsl_vector *x, *w;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  w = xydata->w;
  var = gsl_vector_get(v, 3);
  x0 = gsl_vector_get(v, 2);
  A = gsl_vector_get(v, 1);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yy = exp(-(xi - x0)*(xi - x0)/var/2.0);
    gsl_matrix_set(J, i, 3, A*yy*(xi - x0)*(xi - x0)/2/var/var*wi);
    gsl_matrix_set(J, i, 2, A*yy*(xi - x0)/var*wi);
    gsl_matrix_set(J, i, 1, yy*wi);
    gsl_matrix_set(J, i, 0, 1.0*wi);
  }
  return GSL_SUCCESS;
}

static int Gaussian_fdf(const gsl_vector *v, void *data,
			gsl_vector *f, gsl_matrix *J)
{
  Gaussian_f(v, data, f);
  Gaussian_df(v, data, J);
  return GSL_SUCCESS;
}

/* 2 Gaussian fit 
   y = y0 + A1 exp(-(x-x01/sigma1)^2) + A2 exp(-(x-x02/sigma2)^2)
       v[0] = y0
       v[1] = A1
       v[2] = x01
       v[3] = var1 = sigma1^2
       v[4] = A2
       v[5] = x02
       v[6] = var2 = sigma2^2
*/
static int Gaussian_2peaks_f(const gsl_vector *v, void *data, gsl_vector *f)
{
  gsl_vector *x, *y, *w;
  double y0, A1, x01, var1, A2, x02, var2, xi, yi, wi;
  size_t i;
  struct fitting_xydata *xydata;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  y = xydata->y;
  w = xydata->w;
  y0 = gsl_vector_get(v, 0);
  A1 = gsl_vector_get(v, 1);
  x01 = gsl_vector_get(v, 2);
  var1 = gsl_vector_get(v, 3);
  A2 = gsl_vector_get(v, 4);
  x02 = gsl_vector_get(v, 5);
  var2 = gsl_vector_get(v, 6);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    yi = gsl_vector_get(y, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    gsl_vector_set(f, i, (A1*exp(-(xi - x01)*(xi - x01)/var1/2.0) + A2*exp(-(xi - x02)*(xi - x02)/var2/2.0) + y0 - yi)*wi);
  }
  return GSL_SUCCESS;
}

static int Gaussian_2peaks_df(const gsl_vector *v, void *data, gsl_matrix *J)
{
  double A1, x01, var1, A2, x02, var2, xi, yy, wi;
  size_t i;
  struct fitting_xydata *xydata;
  gsl_vector *x, *w;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  w = xydata->w;
  A1 = gsl_vector_get(v, 1);
  x01 = gsl_vector_get(v, 2);
  var1 = gsl_vector_get(v, 3);
  A2 = gsl_vector_get(v, 4);
  x02 = gsl_vector_get(v, 5);
  var2 = gsl_vector_get(v, 6);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yy = exp(-(xi - x01)*(xi - x01)/var1/2.0);
    gsl_matrix_set(J, i, 0, 1.0*wi);
    gsl_matrix_set(J, i, 1, yy*wi);
    gsl_matrix_set(J, i, 2, A1*yy*(xi - x01)/var1*wi);
    gsl_matrix_set(J, i, 3, A1*yy*(xi - x01)*(xi - x01)/2/var1/var1*wi);

    yy = exp(-(xi - x02)*(xi - x02)/var2/2.0);
    gsl_matrix_set(J, i, 4, yy*wi);
    gsl_matrix_set(J, i, 5, A2*yy*(xi - x02)/var2*wi);
    gsl_matrix_set(J, i, 6, A2*yy*(xi - x02)*(xi - x02)/2/var2/var2*wi);
  }
  return GSL_SUCCESS;
}

static int Gaussian_2peaks_fdf(const gsl_vector *v, void *data,
			gsl_vector *f, gsl_matrix *J)
{
  Gaussian_2peaks_f(v, data, f);
  Gaussian_2peaks_df(v, data, J);
  return GSL_SUCCESS;
}

/* Exponential fit 
   y = y0 + A exp(-bx)
       v[0] = y0
       v[1] = A
       v[2] = b
*/
static int Exponential_f(const gsl_vector *v, void *data, gsl_vector *f)
{
  gsl_vector *x, *y, *w;
  double y0, A, b, xi, yi, wi;
  size_t i;
  struct fitting_xydata *xydata;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  y = xydata->y;
  w = xydata->w;
  y0 = gsl_vector_get(v, 0);
  A = gsl_vector_get(v, 1);
  b = gsl_vector_get(v, 2);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yi = gsl_vector_get(y, i);
    gsl_vector_set(f, i, (A*exp(-b*xi) + y0 - yi)*wi);
  }
  return GSL_SUCCESS;
}

static int Exponential_df(const gsl_vector *v, void *data, gsl_matrix *J)
{
  double A, b, xi, yy, wi;
  size_t i;
  struct fitting_xydata *xydata;
  gsl_vector *x, *w;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  w = xydata->w;
  A = gsl_vector_get(v, 1);
  b = gsl_vector_get(v, 2);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yy = exp(-xi*b);
    gsl_matrix_set(J, i, 0, 1.0*wi);
    gsl_matrix_set(J, i, 1, yy*wi); 
    gsl_matrix_set(J, i, 2, -A*yy*xi*wi);
  }
  return GSL_SUCCESS;
}

static int Exponential_fdf(const gsl_vector *v, void *data,
			gsl_vector *f, gsl_matrix *J)
{
  Exponential_f(v, data, f);
  Exponential_df(v, data, J);
  return GSL_SUCCESS;
}


/* Double Exponential fit 
   y = y0 + A1 exp(-b1 x) + A2 exp(-b2 x)
       v[0] = y0
       v[1] = A1
       v[2] = b1
       v[3] = A2
       v[4] = b2
*/
static int DblExponential_f(const gsl_vector *v, void *data, gsl_vector *f)
{
  gsl_vector *x, *y, *w;
  double y0, A1, b1, A2, b2, xi, yi, wi;
  size_t i;
  struct fitting_xydata *xydata;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  y = xydata->y;
  w = xydata->w;
  y0 = gsl_vector_get(v, 0);
  A1 = gsl_vector_get(v, 1);
  b1 = gsl_vector_get(v, 2);
  A2 = gsl_vector_get(v, 3);
  b2 = gsl_vector_get(v, 4);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yi = gsl_vector_get(y, i);
    gsl_vector_set(f, i, (A1*exp(-b1*xi) + A2*exp(-b2*xi) + y0 - yi)*wi);
  }
  return GSL_SUCCESS;
}

static int DblExponential_df(const gsl_vector *v, void *data, gsl_matrix *J)
{
  double A1, b1, A2, b2, xi, yy1, yy2, wi;
  size_t i;
  struct fitting_xydata *xydata;
  gsl_vector *x, *w;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  w = xydata->w;
  A1 = gsl_vector_get(v, 1);
  b1 = gsl_vector_get(v, 2);
  A2 = gsl_vector_get(v, 3);
  b2 = gsl_vector_get(v, 4);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yy1 = exp(-xi*b1);
    yy2 = exp(-xi*b2);
    gsl_matrix_set(J, i, 0, 1.0*wi);
    gsl_matrix_set(J, i, 1, yy1*wi); 
    gsl_matrix_set(J, i, 2, -A1*yy1*xi*wi);
    gsl_matrix_set(J, i, 3, yy2*wi); 
    gsl_matrix_set(J, i, 4, -A2*yy2*xi*wi);
  }
  return GSL_SUCCESS;
}

static int DblExponential_fdf(const gsl_vector *v, void *data,
			gsl_vector *f, gsl_matrix *J)
{
  DblExponential_f(v, data, f);
  DblExponential_df(v, data, J);
  return GSL_SUCCESS;
}

/* Lorentzian fit 
   y = y0 + A/((x-x0)^2 + B)
       v[0] = y0
       v[1] = A
       v[2] = x0
       v[3] = B
*/
static int Lorentzian_f(const gsl_vector *v, void *data, gsl_vector *f)
{
  gsl_vector *x, *y, *w;
  double y0, A, x0, B, xi, yi, wi;
  size_t i;
  struct fitting_xydata *xydata;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  y = xydata->y;
  w = xydata->w;
  y0 = gsl_vector_get(v, 0);
  A = gsl_vector_get(v, 1);
  x0 = gsl_vector_get(v, 2);
  B = gsl_vector_get(v, 3);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yi = gsl_vector_get(y, i);
    gsl_vector_set(f, i, (A/(gsl_pow_2(xi-x0)+B) + y0 - yi)*wi);
  }
  return GSL_SUCCESS;
}

static int Lorentzian_df(const gsl_vector *v, void *data, gsl_matrix *J)
{
  double A, B, x0, xi, yy, wi;
  size_t i;
  struct fitting_xydata *xydata;
  gsl_vector *x, *w;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  w = xydata->w;
  A = gsl_vector_get(v, 1);
  x0 = gsl_vector_get(v, 2);
  B = gsl_vector_get(v, 3);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yy = gsl_pow_2(xi-x0)+B;
    gsl_matrix_set(J, i, 0, 1.0*wi);
    gsl_matrix_set(J, i, 1, 1.0/yy*wi); 
    gsl_matrix_set(J, i, 2, 2.0*A*(xi-x0)/yy/yy*wi);
    gsl_matrix_set(J, i, 3, -A/yy/yy*wi);
  }
  return GSL_SUCCESS;
}

static int Lorentzian_fdf(const gsl_vector *v, void *data,
			gsl_vector *f, gsl_matrix *J)
{
  Lorentzian_f(v, data, f);
  Lorentzian_df(v, data, J);
  return GSL_SUCCESS;
}

/* Sinusoldal fit
   y = y0 + A sin(fc x + phi)
       v[0] = y0
       v[1] = A
       v[2] = fc
       v[3] = phi
*/
static int Sin_f(const gsl_vector *v, void *data, gsl_vector *f)
{
  gsl_vector *x, *y, *w;
  double y0, A, fc, phi, xi, yi, wi;
  size_t i;
  struct fitting_xydata *xydata;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  y = xydata->y;
  w = xydata->w;
  y0 = gsl_vector_get(v, 0);
  A = gsl_vector_get(v, 1);
  fc = gsl_vector_get(v, 2);
  phi = gsl_vector_get(v, 3);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yi = gsl_vector_get(y, i);
    gsl_vector_set(f, i, (A*sin(fc*xi+phi) + y0 - yi)*wi);
  }
  return GSL_SUCCESS;
}

static int Sin_df(const gsl_vector *v, void *data, gsl_matrix *J)
{
  double A, fc, phi, xi, ys, yc, wi;
  size_t i;
  struct fitting_xydata *xydata;
  gsl_vector *x, *w;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  w = xydata->w;
  A = gsl_vector_get(v, 1);
  fc = gsl_vector_get(v, 2);
  phi = gsl_vector_get(v, 3);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    ys = sin(fc*xi + phi);
    yc = cos(fc*xi + phi);
    gsl_matrix_set(J, i, 0, 1.0*wi);
    gsl_matrix_set(J, i, 1, ys*wi); 
    gsl_matrix_set(J, i, 2, A*yc*xi*wi);
    gsl_matrix_set(J, i, 3, A*yc*wi);
  }
  return GSL_SUCCESS;
}

static int Sin_fdf(const gsl_vector *v, void *data,
			gsl_vector *f, gsl_matrix *J)
{
  Sin_f(v, data, f);
  Sin_df(v, data, J);
  return GSL_SUCCESS;
}

/* Hill's equation fit
   y = y0 + (m - y0)/(1 + (xhalf/x)^r)
       v[0] = y0
       v[1] = m
       v[2] = xhalf
       v[3] = r
*/
static int Hill_f(const gsl_vector *v, void *data, gsl_vector *f)
{
  gsl_vector *x, *y, *w;
  double y0, m, xhalf, r, xi, yi, wi;
  size_t i;
  struct fitting_xydata *xydata;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  y = xydata->y;
  w = xydata->w;
  y0 = gsl_vector_get(v, 0);
  m = gsl_vector_get(v, 1);
  xhalf = gsl_vector_get(v, 2);
  r = gsl_vector_get(v, 3);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yi = gsl_vector_get(y, i);
    gsl_vector_set(f, i, ((m-y0)/(1.0+pow(xhalf/xi, r)) + y0 - yi)*wi);
  }
  return GSL_SUCCESS;
}

static int Hill_df(const gsl_vector *v, void *data, gsl_matrix *J)
{
  double y0, m, xhalf, r, yy, xi, wi, a;
  size_t i;
  struct fitting_xydata *xydata;
  gsl_vector *x, *w;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  w = xydata->w;
  y0 = gsl_vector_get(v, 0);
  m = gsl_vector_get(v, 1);
  xhalf = gsl_vector_get(v, 2);
  r = gsl_vector_get(v, 3);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    a = pow(xhalf/xi, r);
    yy = (1.0 + a);
    gsl_matrix_set(J, i, 0, (1.0 - 1.0/yy)*wi);
    gsl_matrix_set(J, i, 1, 1.0/yy*wi); 
    gsl_matrix_set(J, i, 2, -(m-y0)*r/xhalf*a/yy/yy*wi);
    gsl_matrix_set(J, i, 3, -(m-y0)/yy/yy*a*log(xhalf/xi)*wi);
  }
  return GSL_SUCCESS;
}

static int Hill_fdf(const gsl_vector *v, void *data,
			gsl_vector *f, gsl_matrix *J)
{
  Hill_f(v, data, f);
  Hill_df(v, data, J);
  return GSL_SUCCESS;
}

/* Sigmoidal fit
   y = y0 + m/(1 + exp((x0-x)/r))
       v[0] = y0
       v[1] = m
       v[2] = x0
       v[3] = r
*/
static int Sigmoid_f(const gsl_vector *v, void *data, gsl_vector *f)
{
  gsl_vector *x, *y, *w;
  double y0, m, x0, r, xi, yi, wi;
  size_t i;
  struct fitting_xydata *xydata;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  y = xydata->y;
  w = xydata->w;
  y0 = gsl_vector_get(v, 0);
  m = gsl_vector_get(v, 1);
  x0 = gsl_vector_get(v, 2);
  r = gsl_vector_get(v, 3);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yi = gsl_vector_get(y, i);
    gsl_vector_set(f, i, (m/(1.0+exp((x0-xi)/r)) + y0 - yi)*wi);
  }
  return GSL_SUCCESS;
}

static int Sigmoid_df(const gsl_vector *v, void *data, gsl_matrix *J)
{
  double m, x0, r, xi, wi, a, yy;
  size_t i;
  struct fitting_xydata *xydata;
  gsl_vector *x, *w;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  w = xydata->w;
  m = gsl_vector_get(v, 1);
  x0 = gsl_vector_get(v, 2);
  r = gsl_vector_get(v, 3);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    a = exp((x0 - xi)/r);
    yy = 1.0 + a;
    gsl_matrix_set(J, i, 0, wi);
    gsl_matrix_set(J, i, 1, 1.0/yy*wi); 
    gsl_matrix_set(J, i, 2, -m*a/r/yy/yy*wi);
    gsl_matrix_set(J, i, 3, m*a*(x0-xi)/r/r/yy/yy*wi);
  }
  return GSL_SUCCESS;
}

static int Sigmoid_fdf(const gsl_vector *v, void *data,
			gsl_vector *f, gsl_matrix *J)
{
  Sigmoid_f(v, data, f);
  Sigmoid_df(v, data, J);
  return GSL_SUCCESS;
}


/* Power-law fit
   y = y0 + A x^r
       v[0] = y0
       v[1] = A
       v[2] = r
*/
static int Power_f(const gsl_vector *v, void *data, gsl_vector *f)
{
  gsl_vector *x, *y, *w;
  double y0, A, r, xi, yi, wi;
  size_t i;
  struct fitting_xydata *xydata;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  y = xydata->y;
  w = xydata->w;
  y0 = gsl_vector_get(v, 0);
  A = gsl_vector_get(v, 1);
  r = gsl_vector_get(v, 2);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yi = gsl_vector_get(y, i);
    gsl_vector_set(f, i, (A*pow(xi, r) + y0 - yi)*wi);
  }
  return GSL_SUCCESS;
}

static int Power_df(const gsl_vector *v, void *data, gsl_matrix *J)
{
  double A, r, xi, wi, a;
  size_t i;
  struct fitting_xydata *xydata;
  gsl_vector *x, *w;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  w = xydata->w;
  A = gsl_vector_get(v, 1);
  r = gsl_vector_get(v, 2);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    a = pow(xi, r);
    gsl_matrix_set(J, i, 0, wi);
    gsl_matrix_set(J, i, 1, a*wi); 
    gsl_matrix_set(J, i, 2, A*a*log(xi)*wi);
    
  }
  return GSL_SUCCESS;
}

static int Power_fdf(const gsl_vector *v, void *data,
			gsl_vector *f, gsl_matrix *J)
{
  Power_f(v, data, f);
  Power_df(v, data, J);
  return GSL_SUCCESS;
}


/* Lognormal fit
   y = y0 + A exp[ -(log(x/x0)/width)^2 ]
       v[0] = y0
       v[1] = A
       v[2] = x0
       v[3] = width
*/
static int Lognormal_f(const gsl_vector *v, void *data, gsl_vector *f)
{
  gsl_vector *x, *y, *w;
  double y0, A, x0, width, xi, yi, wi;
  size_t i;
  struct fitting_xydata *xydata;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  y = xydata->y;
  w = xydata->w;
  y0 = gsl_vector_get(v, 0);
  A = gsl_vector_get(v, 1);
  x0 = gsl_vector_get(v, 2);
  width = gsl_vector_get(v, 3);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yi = gsl_vector_get(y, i);
    gsl_vector_set(f, i, (A*exp(-gsl_pow_2(log(xi/x0)/width)) + y0 - yi)*wi);
  }
  return GSL_SUCCESS;
}

static int Lognormal_df(const gsl_vector *v, void *data, gsl_matrix *J)
{
  double A, x0, width, xi, wi, a, b;
  size_t i;
  struct fitting_xydata *xydata;
  gsl_vector *x, *w;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  w = xydata->w;
  A = gsl_vector_get(v, 1);
  x0 = gsl_vector_get(v, 2);
  width = gsl_vector_get(v, 3);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    a = log(xi/x0)/width;
    b = exp(-a*a);
    gsl_matrix_set(J, i, 0, wi);
    gsl_matrix_set(J, i, 1, b*wi); 
    gsl_matrix_set(J, i, 2, 2.0*A*b*a*a*a/width/x0*wi);
    gsl_matrix_set(J, i, 3, 2.0*A*b*a*a*a*a/width*wi);
  }
  return GSL_SUCCESS;
}

static int Lognormal_fdf(const gsl_vector *v, void *data,
			gsl_vector *f, gsl_matrix *J)
{
  Lognormal_f(v, data, f);
  Lognormal_df(v, data, J);
  return GSL_SUCCESS;
}

/* Rayleigh fit 
   y = A exp(-x*x/2/var)
       v[0] = A
       v[1] = var = sigma^2
*/
static int Rayleigh_f(const gsl_vector *v, void *data, gsl_vector *f)
{
  gsl_vector *x, *y, *w;
  double A, var, xi, yi, wi;
  size_t i;
  struct fitting_xydata *xydata;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  y = xydata->y;
  w = xydata->w;
  var = gsl_vector_get(v, 1);
  A = gsl_vector_get(v, 0);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    yi = gsl_vector_get(y, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    gsl_vector_set(f, i, (A*xi*exp(-xi*xi/var/2.0) - yi)*wi);
  }
  return GSL_SUCCESS;
}

static int Rayleigh_df(const gsl_vector *v, void *data, gsl_matrix *J)
{
  double A, var, xi, yy, wi;
  size_t i;
  struct fitting_xydata *xydata;
  gsl_vector *x, *w;
  xydata = (struct fitting_xydata*) data;
  x = xydata->x;
  w = xydata->w;
  var = gsl_vector_get(v, 1);
  A = gsl_vector_get(v, 0);
  for (i = 0; i < x->size; i++) {
    xi = gsl_vector_get(x, i);
    if (w) wi = gsl_vector_get(w, i);
    else wi = 1.0;
    yy = xi*exp(-xi*xi/var/2.0);
    gsl_matrix_set(J, i, 1, A*yy*xi*xi/2/var/var*wi);
    gsl_matrix_set(J, i, 0, yy*wi);
  }
  return GSL_SUCCESS;
}

static int Rayleigh_fdf(const gsl_vector *v, void *data,
			gsl_vector *f, gsl_matrix *J)
{
  Rayleigh_f(v, data, f);
  Rayleigh_df(v, data, J);
  return GSL_SUCCESS;
}

static void set_fittype(gsl_multifit_function_fdf *f, const char *fittype, size_t *p, gsl_vector **v, int *flag)
{
  if (str_tail_grep(fittype, "aussian_2peaks") == 0) {
    f->f = Gaussian_2peaks_f;
    f->df = Gaussian_2peaks_df;
    f->fdf = Gaussian_2peaks_fdf;
    *p = 7;
    if (*v == NULL) {
      *v = gsl_vector_alloc(*p);
      gsl_vector_set(*v, 0, 0);   /* y0 = 0 */
      gsl_vector_set(*v, 1, 1);   /* A = 1 */
      gsl_vector_set(*v, 2, 0.0); /* x0 = 0 */
      gsl_vector_set(*v, 3, 1);   /* initial values, var = 1 */
      gsl_vector_set(*v, 4, 1);   /* A = 1 */
      gsl_vector_set(*v, 5, 1.0); /* x0 = 1 */
      gsl_vector_set(*v, 6, 1);   /* initial values, var = 1 */
      *flag = 1;
    }
  } else if (str_head_grep(fittype, "gaus") == 0) {
    f->f = Gaussian_f;
    f->df = Gaussian_df;
    f->fdf = Gaussian_fdf;
    *p = 4;
    if (*v == NULL) {
      *v = gsl_vector_alloc(*p);
      gsl_vector_set(*v, 3, 1);   /* initial values, var = 1 */
      gsl_vector_set(*v, 2, 0.0); /* x0 = 0 */
      gsl_vector_set(*v, 1, 1);   /* A = 1 */
      gsl_vector_set(*v, 0, 0);   /* y0 = 0 */
      *flag = 1;
    }
  } else if (str_head_grep(fittype, "exp") == 0) {
    f->f = Exponential_f;
    f->df = Exponential_df;
    f->fdf = Exponential_fdf;
    *p = 3;
    if (*v == NULL) {
      *v = gsl_vector_alloc(*p);
      gsl_vector_set(*v, 0, 0);   /* y0 = 0 */
      gsl_vector_set(*v, 1, 1);   /* A = 1 */
      gsl_vector_set(*v, 2, 1);   
      *flag = 1;
    }
 } else if (str_head_grep(fittype, "rayleigh") == 0) {
    f->f = Rayleigh_f;
    f->df = Rayleigh_df;
    f->fdf = Rayleigh_fdf;
    *p = 2;
    if (*v == NULL) {
      *v = gsl_vector_alloc(*p);
      gsl_vector_set(*v, 0, 1);   /* A = 1 */
      gsl_vector_set(*v, 1, 1);   /* sigma = 1 */
      *flag = 1;
    }
  } else if (str_head_grep(fittype, "dblexp") == 0) {
    f->f = DblExponential_f;
    f->df = DblExponential_df;
    f->fdf = DblExponential_fdf;
    *p = 5;
    if (*v == NULL) {
      *v = gsl_vector_alloc(*p);
      gsl_vector_set(*v, 0, 0);   /* y0 = 0 */
      gsl_vector_set(*v, 1, 1);   /* A = 1 */
      gsl_vector_set(*v, 2, 1);   
      gsl_vector_set(*v, 3, 1);   /* A = 1 */
      gsl_vector_set(*v, 4, 1);   
      *flag = 1;
    }
  } else if (str_head_grep(fittype, "lor") == 0) {
    f->f = Lorentzian_f;
    f->df = Lorentzian_df;
    f->fdf = Lorentzian_fdf;
    *p = 4;
    if (*v == NULL) {
      *v = gsl_vector_alloc(*p);
      gsl_vector_set(*v, 0, 0); 
      gsl_vector_set(*v, 1, 1);
      gsl_vector_set(*v, 2, 0);  
      gsl_vector_set(*v, 3, 1);  
      *flag = 1;
    }
  } else if (str_head_grep(fittype, "sin") == 0) {
    f->f = Sin_f;
    f->df = Sin_df;
    f->fdf = Sin_fdf;
    *p = 4;
    if (*v == NULL) {
      *v = gsl_vector_alloc(*p);
      gsl_vector_set(*v, 0, 0); 
      gsl_vector_set(*v, 1, 1);
      gsl_vector_set(*v, 2, 1);  
      gsl_vector_set(*v, 3, 0);  
      *flag = 1;
    }
  } else if (str_head_grep(fittype, "hill") == 0) {
    f->f = Hill_f;
    f->df = Hill_df;
    f->fdf = Hill_fdf;
    *p = 4;
    if (*v == NULL) {
      *v = gsl_vector_alloc(*p);
      gsl_vector_set(*v, 0, 0); 
      gsl_vector_set(*v, 1, 1);
      gsl_vector_set(*v, 2, 1);  
      gsl_vector_set(*v, 3, 1);  
      *flag = 1;
    }
  } else if (str_head_grep(fittype, "sigmoid") == 0 || str_head_grep(fittype, "fermi") == 0) {
    f->f = Sigmoid_f;
    f->df = Sigmoid_df;
    f->fdf = Sigmoid_fdf;
    *p = 4;
    if (*v == NULL) {
      *v = gsl_vector_alloc(*p);
      gsl_vector_set(*v, 0, 0); 
      gsl_vector_set(*v, 1, 1);
      gsl_vector_set(*v, 2, 0);  
      gsl_vector_set(*v, 3, 1);  
      *flag = 1;
    }
  } else if (str_head_grep(fittype, "power") == 0) {
    f->f = Power_f;
    f->df = Power_df;
    f->fdf = Power_fdf;
    *p = 3;
    if (*v == NULL) {
      *v = gsl_vector_alloc(*p);
      gsl_vector_set(*v, 0, 0); 
      gsl_vector_set(*v, 1, 1);
      gsl_vector_set(*v, 2, -1);  
      *flag = 1;
    }
  } else if (str_head_grep(fittype, "lognormal") == 0) {
    f->f = Lognormal_f;
    f->df = Lognormal_df;
    f->fdf = Lognormal_fdf;
    *p = 4;
    if (*v == NULL) {
      *v = gsl_vector_alloc(*p);
      gsl_vector_set(*v, 0, 0); 
      gsl_vector_set(*v, 1, 1);
      gsl_vector_set(*v, 2, 1);  
      gsl_vector_set(*v, 3, 1);  
      *flag = 1;
    }
  } else {
    rb_raise(rb_eRuntimeError, "Unknown fit type (gaussian expected)");
  }

}

/* Singleton method */
static VALUE rb_gsl_multifit_fit(int argc, VALUE *argv, VALUE klass)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *solver;
  int status;
  size_t iter = 0, i;
  size_t n, dof;      /* # of data points */
  size_t p;           /* # of fitting parameters */
  gsl_multifit_function_fdf f;
  gsl_matrix *covar = NULL;
  gsl_vector *v = NULL;
  gsl_vector *x, *y, *w = NULL;
  gsl_vector_view xx, yy, ww;
  gsl_vector *vout, *verr;
  int flag = 0;
  double chi2;
  char fittype[256];
  struct fitting_xydata xydata;
  if (argc < 3) rb_raise(rb_eArgError, "too few arguments");
  switch (TYPE(argv[argc-1])) {
  case T_ARRAY:
    v = get_vector(argv[argc-1]);
    flag = 1;
    argc--;
    break;
  case T_STRING:
    /* do nothing */
    break;
  default:
    Data_Get_Vector(argv[argc-1], v);
    flag = 0;
    argc--;
    break;
  }
  x = &xx.vector;
  y = &yy.vector;
  w = &ww.vector;
  switch (argc) {
  case 3:
    Data_Get_Vector(argv[0], x);
    Data_Get_Vector(argv[1], y);
    w = NULL;
    strcpy(fittype, STR2CSTR(argv[2]));
    break;
  case 4:
    Data_Get_Vector(argv[0], x);
    Data_Get_Vector(argv[1], w);
    Data_Get_Vector(argv[2], y);
    strcpy(fittype, STR2CSTR(argv[3]));
    break;
  default:
    rb_raise(rb_eArgError, "wrong number of arguments");
    break;
  }

  xydata.x = x;
  xydata.y = y;
  xydata.w = w;
  n = x->size;

  set_fittype(&f, fittype, &p, &v, &flag);

  f.n = n;
  f.p = p;
  f.params = &xydata;

  T = gsl_multifit_fdfsolver_lmsder;
  solver = gsl_multifit_fdfsolver_alloc(T, n, p);
  gsl_multifit_fdfsolver_set(solver, &f, v);

  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);
    if (status) break;
    status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-6, 1e-6);
  } while (status == GSL_CONTINUE);
  vout = gsl_vector_alloc(p);
  verr = gsl_vector_alloc(p);
  gsl_vector_memcpy(vout, solver->x);
  covar = gsl_matrix_alloc(p, p);
  chi2 = gsl_pow_2(gsl_blas_dnrm2(solver->f));   /* not reduced chi-square */
  dof = n - p;
  gsl_multifit_covar(solver->J, 0.0, covar);
  for (i = 0; i < p; i++)
    gsl_vector_set(verr, i, sqrt(chi2/dof*gsl_matrix_get(covar, i, i)));
  gsl_matrix_free(covar);
  if (flag == 1) gsl_vector_free(v);
  gsl_multifit_fdfsolver_free(solver);
  return rb_ary_new3(4,
		     Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, vout),
		     Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, verr),
		     rb_float_new(chi2), INT2FIX(dof));
}

#ifdef GSL_1_8_LATER
static VALUE rb_gsl_multifit_linear_est(VALUE module, VALUE xx, VALUE cc, VALUE ccov)
{
  gsl_vector *x, *c;
  gsl_matrix *cov;
  double y, y_err;

  Data_Get_Vector(xx, x);
  Data_Get_Vector(cc, c);
  Data_Get_Matrix(ccov, cov);
  gsl_multifit_linear_est(x, c, cov, &y, &y_err);
  return rb_ary_new3(2, rb_float_new(y), rb_float_new(y_err));
}
#endif
#ifdef GSL_1_11_LATER
static VALUE rb_gsl_multifit_linear_residuals(int argc, VALUE argv[], VALUE module)
{
  gsl_vector *y, *c, *r;
  gsl_matrix *X;
  VALUE ret;
  switch (argc) {
  case 3:
    Data_Get_Matrix(argv[0], X);
    Data_Get_Vector(argv[1], y);
    Data_Get_Vector(argv[2], c);
    r = gsl_vector_alloc(y->size);
    ret = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, r);
    break;
  case 4:
    Data_Get_Matrix(argv[0], X);
    Data_Get_Vector(argv[1], y);
    Data_Get_Vector(argv[2], c);
    Data_Get_Vector(argv[3], r);
    ret = argv[3];
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments %d (3 or 4).\n", argc);
  }
  gsl_multifit_linear_residuals(X, y, c, r);
  return ret;
}

#endif

void Init_gsl_multifit(VALUE module)
{
  VALUE mgsl_multifit;
  VALUE cgsl_multifit_fdfsolver;
  VALUE cgsl_multifit_solver;

  mgsl_multifit = rb_define_module_under(module, "MultiFit");

  cgsl_multifit_workspace = rb_define_class_under(mgsl_multifit, "Workspace", 
						  cGSL_Object);
  rb_define_singleton_method(cgsl_multifit_workspace, "new", 
			     rb_gsl_multifit_workspace_new, 2);
  rb_define_singleton_method(cgsl_multifit_workspace, "alloc", 
			     rb_gsl_multifit_workspace_new, 2);
  rb_define_singleton_method(mgsl_multifit, "alloc", 
			     rb_gsl_multifit_workspace_new, 2);			     
  rb_define_module_function(mgsl_multifit, "linear", rb_gsl_multifit_linear, -1);
  rb_define_module_function(mgsl_multifit, "wlinear", rb_gsl_multifit_wlinear, -1);

  cgsl_multifit_solver = rb_define_class_under(mgsl_multifit, "Solver", cGSL_Object);
  cgsl_multifit_fdfsolver = rb_define_class_under(mgsl_multifit, "FdfSolver", 
						  cgsl_multifit_solver);
  cgsl_multifit_function_fdf = rb_define_class_under(mgsl_multifit, "Function_fdf",
						     cGSL_Object);

  /*****/

  rb_define_singleton_method(cgsl_multifit_fdfsolver, "new", 
			     rb_gsl_multifit_fdfsolver_new, -1);
  rb_define_singleton_method(cgsl_multifit_fdfsolver, "alloc", 
			     rb_gsl_multifit_fdfsolver_new, -1);

  rb_define_singleton_method(cgsl_multifit_function_fdf, "new", 
			     rb_gsl_multifit_function_fdf_new, -1);
  rb_define_singleton_method(cgsl_multifit_function_fdf, "alloc", 
			     rb_gsl_multifit_function_fdf_new, -1);

  /*****/

  rb_define_method(cgsl_multifit_fdfsolver, "set", rb_gsl_multifit_fdfsolver_set, 2);
  rb_define_method(cgsl_multifit_fdfsolver, "name", rb_gsl_multifit_fdfsolver_name, 0);
  rb_define_method(cgsl_multifit_fdfsolver, "iterate", rb_gsl_multifit_fdfsolver_iterate, 0);
  rb_define_method(cgsl_multifit_fdfsolver, "position", rb_gsl_multifit_fdfsolver_position, 0);
  //  rb_define_alias(cgsl_multifit_fdfsolver, "x", "position");
  rb_define_method(cgsl_multifit_fdfsolver, "print_state", rb_gsl_multifit_fdfsolver_print_state, 1);
  rb_define_method(cgsl_multifit_fdfsolver, "fdf", rb_gsl_multifit_fdfsolver_fdf, 0);
  rb_define_method(cgsl_multifit_fdfsolver, "test_delta", rb_gsl_multifit_fdfsolver_test_delta, 2);
  rb_define_method(cgsl_multifit_fdfsolver, "test_gradient", rb_gsl_multifit_fdfsolver_test_gradient, -1);
  rb_define_method(cgsl_multifit_fdfsolver, "covar", rb_gsl_multifit_fdfsolver_covar, -1);
  rb_define_method(cgsl_multifit_fdfsolver, "gradient", rb_gsl_multifit_fdfsolver_gradient, -1);
  rb_define_method(cgsl_multifit_fdfsolver, "x", rb_gsl_multifit_fdfsolver_x, 0);
  rb_define_method(cgsl_multifit_fdfsolver, "dx", rb_gsl_multifit_fdfsolver_dx, 0);
  rb_define_method(cgsl_multifit_fdfsolver, "f", rb_gsl_multifit_fdfsolver_f, 0);
  rb_define_method(cgsl_multifit_fdfsolver, "J", rb_gsl_multifit_fdfsolver_J, 0);
  rb_define_alias(cgsl_multifit_fdfsolver, "jac", "J");
  rb_define_alias(cgsl_multifit_fdfsolver, "jacobian", "J");
  rb_define_alias(cgsl_multifit_fdfsolver, "Jacobian", "J");

  /*****/
  rb_define_module_function(mgsl_multifit, "polyfit", rb_gsl_multifit_polyfit, -1);
  rb_define_module_function(mgsl_multifit, "legfit", rb_gsl_multifit_legfit, -1);
  /*****/

  rb_define_singleton_method(mgsl_multifit, "test_delta", rb_gsl_multifit_test_delta, 4);
  rb_define_singleton_method(mgsl_multifit, "test_gradient", rb_gsl_multifit_test_gradient, 2);
  rb_define_singleton_method(mgsl_multifit, "gradient", rb_gsl_multifit_gradient, -1);
  rb_define_singleton_method(mgsl_multifit, "covar", rb_gsl_multifit_covar, -1);

  rb_gsl_multifit_define_const(cgsl_multifit_fdfsolver);

  rb_define_method(cgsl_multifit_function_fdf, "set_procs", 
		   rb_gsl_multifit_function_fdf_set_procs, -1);
  rb_define_method(cgsl_multifit_function_fdf, "set_data", 
		   rb_gsl_multifit_function_fdf_set_data, -1);
  rb_define_method(cgsl_multifit_function_fdf, "params", 
		   rb_gsl_multifit_function_fdf_params, 0);
  rb_define_alias(cgsl_multifit_function_fdf, "param", "params");

  rb_define_method(cgsl_multifit_function_fdf, "n", 
		   rb_gsl_multifit_function_fdf_n, 0);
  rb_define_method(cgsl_multifit_function_fdf, "set_n", 
		   rb_gsl_multifit_function_fdf_set_n, 1);
  rb_define_alias(cgsl_multifit_function_fdf, "n=", "set_n");
  rb_define_method(cgsl_multifit_function_fdf, "p", 
		   rb_gsl_multifit_function_fdf_p, 0);
  rb_define_alias(cgsl_multifit_function_fdf, "np", "p");

  /*****/
  rb_define_singleton_method(cgsl_multifit_fdfsolver, "fit", rb_gsl_multifit_fit, -1);

  /***/
#ifdef GSL_1_8_LATER
  rb_define_module_function(mgsl_multifit, "linear_est", rb_gsl_multifit_linear_est, 3);
  rb_define_module_function(module, "multifit_linear_est", rb_gsl_multifit_linear_est, 3);
#endif
#ifdef GSL_1_11_LATER
  rb_define_module_function(mgsl_multifit, "linear_residuals", rb_gsl_multifit_linear_residuals, -1);
  rb_define_module_function(module, "multifit_linear_residuals", rb_gsl_multifit_linear_residuals, -1);
#endif

#ifdef HAVE_NDLINEAR_GSL_MULTIFIT_NDLINEAR_H
  Init_ndlinear(mgsl_multifit);
#endif

}

#ifdef CHECK_WORKSPACE
#undef CHECK_WORKSPACE
#endif
#ifdef CHECK_MULTIFIT_FUNCTION_FDF
#undef CHECK_MULTIFIT_FUNCTION_FDF
#endif
