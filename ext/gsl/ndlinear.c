#include "rb_gsl.h"

#ifdef HAVE_NDLINEAR_GSL_MULTIFIT_NDLINEAR_H
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>
#include <ndlinear/gsl_multifit_ndlinear.h>

static VALUE cWorkspace;

enum Index_Ndlinear {
  INDEX_NDIM = 0,
  INDEX_N = 1,
  INDEX_PROCS = 2,
  INDEX_PARAMS = 3,
  INDEX_FUNCS = 4,
  INDEX_NDIM_I = 5,  
  
  NDLINEAR_ARY_SIZE = 6,
};

static void multifit_ndlinear_mark(gsl_multifit_ndlinear_workspace *w)
{
  rb_gc_mark((VALUE) w->params); 
}

typedef int (*UFUNC)(double, double[], void*);
typedef struct ufunc_struct
{
  UFUNC *fptr;
} ufunc_struct;

static VALUE cUFunc;
static ufunc_struct* ufunc_struct_alloc(size_t n_dim) {
  ufunc_struct *p;
  p = (ufunc_struct*) malloc(sizeof(ufunc_struct));
  p->fptr = malloc(sizeof(UFUNC)*n_dim);  
  return p;
}
static void ufunc_struct_free(ufunc_struct *p)
{
  free(p->fptr);
  free(p);
}

static int func_u(double x, double y[], void *data);
static VALUE rb_gsl_multifit_ndlinear_alloc(int argc, VALUE *argv, VALUE klass)
{
  gsl_multifit_ndlinear_workspace *w;
  int istart = 0;
  size_t n_dim = 0, *N, i;
  struct ufunc_struct *p;
  VALUE params, wspace, pp;
  switch (argc) {
  case 4:
    istart = 1;
    CHECK_FIXNUM(argv[0]);
    n_dim = FIX2INT(argv[0]);
    /* no break */
  case 3:  
    if (TYPE(argv[istart]) != T_ARRAY) {
      rb_raise(rb_eTypeError, "Wrong argument type %s (Array expected)",
        rb_class2name(CLASS_OF(argv[istart])));
    }
    if (TYPE(argv[istart+1]) != T_ARRAY) {
      rb_raise(rb_eTypeError, "Wrong argument type %s (Array expected)",
        rb_class2name(CLASS_OF(argv[istart+1])));
    }
    //    n_dim = RARRAY(argv[istart])->len;
    n_dim = RARRAY_LEN(argv[istart]);
    N = (size_t*) malloc(sizeof(size_t)*n_dim);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 3 or 4)", argc);
  }
  for (i = 0; i < n_dim; i++) {
    N[i] = FIX2INT(rb_ary_entry(argv[istart], i));
  }

  params = rb_ary_new2(NDLINEAR_ARY_SIZE);
  rb_ary_store(params, INDEX_NDIM, INT2FIX((int) n_dim));
  rb_ary_store(params, INDEX_N, argv[istart]);   /* N */
  rb_ary_store(params, INDEX_PROCS, argv[istart+1]); /* procs */
  rb_ary_store(params, INDEX_PARAMS, argv[istart+2]); /* params */  
  rb_ary_store(params, INDEX_NDIM_I, INT2FIX(0)); /* for the first parameter */    
  
  p = ufunc_struct_alloc(n_dim);
  for (i = 0; i < n_dim; i++) p->fptr[i] = func_u;
  pp = Data_Wrap_Struct(cUFunc, 0, ufunc_struct_free, p);  
  rb_ary_store(params, INDEX_FUNCS, pp);  

  w = gsl_multifit_ndlinear_alloc(n_dim, N, p->fptr, (void*) params);
    
  free((size_t*) N);

  wspace = Data_Wrap_Struct(cWorkspace, multifit_ndlinear_mark, gsl_multifit_ndlinear_free, w);

  return wspace;
}

static int func_u(double x, double y[], void *data)
{
  VALUE ary, vN, procs, proc, vy, params;
  gsl_vector_view ytmp;
  size_t i, n_dim;
  int rslt;
  ary = (VALUE) data;
  n_dim = FIX2INT(rb_ary_entry(ary, INDEX_NDIM));
  vN = rb_ary_entry(ary, INDEX_N);
  procs = rb_ary_entry(ary, INDEX_PROCS);
  params = rb_ary_entry(ary, INDEX_PARAMS);
  i = FIX2INT(rb_ary_entry(ary, INDEX_NDIM_I));
  proc = rb_ary_entry(procs, i);
  
  ytmp.vector.data = (double*) y;
  ytmp.vector.stride = 1;
  ytmp.vector.size = FIX2INT(rb_ary_entry(vN, i));
  vy = Data_Wrap_Struct(cgsl_vector_view, 0, NULL, &ytmp);

  rslt = rb_funcall((VALUE) proc, RBGSL_ID_call, 3, rb_float_new(x), vy, params);

  /* for the next parameter */
  i += 1;
  if (i == n_dim) i = 0;
  rb_ary_store(ary, INDEX_NDIM_I, INT2FIX(i));
  
  return GSL_SUCCESS;
}

static VALUE rb_gsl_multifit_ndlinear_design(int argc, VALUE *argv, VALUE obj)
{
  gsl_multifit_ndlinear_workspace *w;
  gsl_matrix *vars = NULL, *X = NULL;
  int argc2, flag = 0, ret;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (!rb_obj_is_kind_of(argv[argc-1], cWorkspace)) {
      rb_raise(rb_eTypeError, "Wrong argument type %s (GSL::MultiFit::Ndlinear::Workspace expected)",
        rb_class2name(CLASS_OF(argv[argc-1])));
    }
    Data_Get_Struct(argv[argc-1], gsl_multifit_ndlinear_workspace, w);
    argc2 = argc-1;
    break;
  default:
    Data_Get_Struct(obj, gsl_multifit_ndlinear_workspace, w);
    argc2 = argc;
  }
  switch (argc2) {
  case 1:
      CHECK_MATRIX(argv[0]);
      Data_Get_Struct(argv[0], gsl_matrix, vars);
      X = gsl_matrix_alloc(vars->size1, w->n_coeffs);
      flag = 1;
      break;
  case 2:
      CHECK_MATRIX(argv[0]);
      CHECK_MATRIX(argv[1]);
      Data_Get_Struct(argv[0], gsl_matrix, vars);
      Data_Get_Struct(argv[1], gsl_matrix, X);            
      break;
  default:
      rb_raise(rb_eArgError, "Wrong number of arguments.");
  }
  ret = gsl_multifit_ndlinear_design(vars, X, w);
  
  if (flag == 1) {
    return Data_Wrap_Struct(cgsl_matrix, 0, gsl_matrix_free, X);
  } else {
    return INT2FIX(ret);
  }
}

static VALUE rb_gsl_multifit_ndlinear_est(int argc, VALUE *argv, VALUE obj)
{
  gsl_multifit_ndlinear_workspace *w;
  gsl_vector *x = NULL, *c = NULL;
  gsl_matrix *cov = NULL;
  double y, yerr;
  int argc2;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (!rb_obj_is_kind_of(argv[argc-1], cWorkspace)) {
      rb_raise(rb_eTypeError, "Wrong argument type %s (GSL::MultiFit::Ndlinear::Workspace expected)",
        rb_class2name(CLASS_OF(argv[argc-1])));
    }
    Data_Get_Struct(argv[argc-1], gsl_multifit_ndlinear_workspace, w);    
    argc2 = argc-1;
    break;
  default:
    Data_Get_Struct(obj, gsl_multifit_ndlinear_workspace, w);
    argc2 = argc;  
  }
  switch (argc2) {
  case 3:
    CHECK_VECTOR(argv[0]);
    CHECK_VECTOR(argv[1]);
    CHECK_MATRIX(argv[2]);
    Data_Get_Struct(argv[0], gsl_vector, x);
    Data_Get_Struct(argv[1], gsl_vector, c);
    Data_Get_Struct(argv[2], gsl_matrix, cov);   
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments.");  
  }
  gsl_multifit_ndlinear_est(x, c, cov, &y, &yerr, w);
  return rb_ary_new3(2, rb_float_new(y), rb_float_new(yerr));
}

static VALUE rb_gsl_multifit_ndlinear_calc(int argc, VALUE *argv, VALUE obj)
{
  gsl_multifit_ndlinear_workspace *w;
  gsl_vector *x = NULL, *c = NULL;
  double val;
  int argc2;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (!rb_obj_is_kind_of(argv[argc-1], cWorkspace)) {
      rb_raise(rb_eTypeError, 
	       "Wrong argument type %s (GSL::MultiFit::Ndlinear::Workspace expected)",
        rb_class2name(CLASS_OF(argv[argc-1])));
    }
    Data_Get_Struct(argv[argc-1], gsl_multifit_ndlinear_workspace, w);    
    argc2 = argc-1;
    break;
  default:
    Data_Get_Struct(obj, gsl_multifit_ndlinear_workspace, w);
    argc2 = argc;  
  }
  switch (argc2) {
  case 2:
    CHECK_VECTOR(argv[0]);
    CHECK_VECTOR(argv[1]);
    Data_Get_Struct(argv[0], gsl_vector, x);
    Data_Get_Struct(argv[1], gsl_vector, c);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments.");  
  }
  val = gsl_multifit_ndlinear_calc(x, c, w);
  return rb_float_new(val);
}

static VALUE rb_gsl_multifit_ndlinear_n_coeffs(VALUE obj)
{
  gsl_multifit_ndlinear_workspace *w;
  Data_Get_Struct(obj, gsl_multifit_ndlinear_workspace, w);
  return INT2FIX(w->n_coeffs);
}

static VALUE rb_gsl_multifit_ndlinear_n_dim(VALUE obj)
{
  gsl_multifit_ndlinear_workspace *w;
  Data_Get_Struct(obj, gsl_multifit_ndlinear_workspace, w);
  return INT2FIX(w->n_dim);
}

static VALUE rb_gsl_multifit_ndlinear_N(VALUE obj)
{
  gsl_multifit_ndlinear_workspace *w;
  VALUE ary;
  Data_Get_Struct(obj, gsl_multifit_ndlinear_workspace, w);
  ary = (VALUE) w->params;
  return rb_ary_entry(ary, INDEX_N);
}
/*
static VALUE rb_gsl_multifit_linear_Rsq(VALUE module, VALUE vy, VALUE vchisq)
{
  gsl_vector *y;
  double chisq, Rsq;
  CHECK_VECTOR(vy);
  Data_Get_Struct(vy, gsl_vector, y);
  chisq = NUM2DBL(vchisq);
  gsl_multifit_linear_Rsq(y, chisq, &Rsq);
  return rb_float_new(Rsq);
}
*/
void Init_ndlinear(VALUE module)
{
  VALUE mNdlinear;
  mNdlinear = rb_define_module_under(module, "Ndlinear");
  cUFunc = rb_define_class_under(mNdlinear, "UFunc", rb_cObject);
  cWorkspace = rb_define_class_under(mNdlinear, "Workspace", cGSL_Object);
  
  rb_define_singleton_method(mNdlinear, "alloc", 
                            rb_gsl_multifit_ndlinear_alloc, -1);
  rb_define_singleton_method(cWorkspace, "alloc", 
                            rb_gsl_multifit_ndlinear_alloc, -1);    
                            
  rb_define_singleton_method(mNdlinear, "design", 
                            rb_gsl_multifit_ndlinear_design, -1);
  rb_define_singleton_method(cWorkspace, "design", 
                            rb_gsl_multifit_ndlinear_design, -1);
  rb_define_method(cWorkspace, "design",rb_gsl_multifit_ndlinear_est, -1);
  rb_define_singleton_method(mNdlinear, "est", 
                            rb_gsl_multifit_ndlinear_est, -1);
  rb_define_singleton_method(cWorkspace, "est", 
                            rb_gsl_multifit_ndlinear_est, -1);
  rb_define_method(cWorkspace, "est",rb_gsl_multifit_ndlinear_est, -1);  
  
  rb_define_singleton_method(mNdlinear, "calc", 
                            rb_gsl_multifit_ndlinear_calc, -1);
  rb_define_singleton_method(cWorkspace, "calc", 
                            rb_gsl_multifit_ndlinear_calc, -1);
  rb_define_method(cWorkspace, "calc",rb_gsl_multifit_ndlinear_calc, -1);  

  rb_define_method(cWorkspace, "n_coeffs",rb_gsl_multifit_ndlinear_n_coeffs, 0);    
  rb_define_method(cWorkspace, "n_dim",rb_gsl_multifit_ndlinear_n_dim, 0);      
  rb_define_method(cWorkspace, "N",rb_gsl_multifit_ndlinear_N, 0);
  
  //  rb_define_module_function(module, "linear_Rsq", rb_gsl_multifit_linear_Rsq, 2);
}

#endif

