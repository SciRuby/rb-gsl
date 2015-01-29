#ifdef HAVE_JACOBI_H
#include "include/rb_gsl.h"
#include "jacobi.h"

static VALUE jac_eval3_e(VALUE x, VALUE a, VALUE b,
                    int (*f)(double, double, double, gsl_sf_result*))
{
  gsl_sf_result *result;
  VALUE obj;
  obj = Data_Make_Struct(cgsl_sf_result, gsl_sf_result, 0, free, result);
  (*f)(NUM2DBL(x), NUM2DBL(a), NUM2DBL(b), result);
  return obj;
}

static VALUE jac_eval3(VALUE xx, VALUE aa, VALUE bb, double (*f)(double, double, double))
{
  gsl_vector *x, *y;
  double a, b;
  size_t i, len;
  VALUE ary;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif  
  a = NUM2DBL(aa);
  b = NUM2DBL(bb);
  if (VECTOR_P(xx)) {
    Data_Get_Struct(xx, gsl_vector, x);
    y = gsl_vector_alloc(x->size);
    for (i = 0; i < x->size; i++) {
      gsl_vector_set(y, i, (*f)(gsl_vector_get(x, i), a, b));
    }
    return Data_Wrap_Struct(VECTOR_ROW_COL(CLASS_OF(xx)), 0, gsl_vector_free, y);
  } else if (TYPE(xx) == T_ARRAY) {
    //    len = RARRAY(xx)->len;
    len = RARRAY_LEN(xx);
    ary = rb_ary_new2(len);
    for (i = 0; i < len; i++) {
      rb_ary_store(ary, i, rb_float_new((*f)(NUM2DBL(rb_ary_entry(xx, i)), a, b)));
    }
    return ary;
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(xx)) {
    GetNArray(xx, na);
    len = na->total;
    ptr1 = (double*) na->ptr;
    ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(xx));
    ptr2 = NA_PTR_TYPE(ary, double*);
    for (i = 0; i < len; i++) {
      ptr2[i] = (*f)(ptr1[i], a, b);
    }    
    return ary;
#endif    
  } else {
    return rb_float_new((*f)(NUM2DBL(xx), a, b));
  }  
}

static VALUE rb_jac_jacobi_eval(int argc, VALUE *argv, 
          double (*f)(double, int, double, double),
          int (*f2)(int, const double*, int, double*, double, double, double*))
{
  gsl_vector *x, *ws, *y;
  double a, b;
  VALUE ary;
  size_t len, i;
  int n, flag = 0;
#ifdef HAVE_NARRAY_H
  double *ptr1, *ptr2;
  struct NARRAY *na;
#endif
  if (argc < 4) rb_raise(rb_eArgError, "Too few arguments (%d for >= 4)", argc);
  if (VECTOR_P(argv[0])) {
    Data_Get_Struct(argv[0], gsl_vector, x);
    y = gsl_vector_alloc(x->size);
    ary = Data_Wrap_Struct(VECTOR_ROW_COL(CLASS_OF(x)), 0, gsl_vector_free, y);
    switch (argc) {
    case 4:
      ws = gsl_vector_alloc(2*x->size);
      flag = 1;
      break;
    case 5:
      CHECK_VECTOR(argv[4]);
      Data_Get_Struct(argv[4], gsl_vector, ws);
      break;
    default:  
      rb_raise(rb_eArgError, "Too many arguments (%d for 4 or 5)", argc);
    }
    (*f2)(x->size, x->data, FIX2INT(argv[1]), y->data, NUM2DBL(argv[2]), NUM2DBL(argv[3]), 
                    ws->data);
    if (flag == 1) gsl_vector_free(ws);
    return ary;
  } else if (TYPE(argv[0]) == T_ARRAY) {
    n = FIX2INT(argv[1]);
    a = NUM2DBL(argv[2]);
    b = NUM2DBL(argv[3]);
    //    len = RARRAY(argv[0])->len;
    len = RARRAY_LEN(argv[0]);
    ary = rb_ary_new2(len);
    for (i = 0; i < len; i++) {
      rb_ary_store(ary, i, rb_float_new((*f)(NUM2DBL(rb_ary_entry(argv[0], i)), n, a, b)));
    }
    return ary;
#ifdef HAVE_NARRAY_H
  } else if (NA_IsNArray(argv[0])) {
    GetNArray(argv[0], na);
    len = na->total;
    ptr1 = (double*) na->ptr;
    ary = na_make_object(NA_DFLOAT, na->rank, na->shape, CLASS_OF(argv[0]));
    ptr2 = NA_PTR_TYPE(ary, double*);
    n = FIX2INT(argv[1]);
    a = NUM2DBL(argv[2]);
    b = NUM2DBL(argv[3]);    
    ws = gsl_vector_alloc(len);
    (*f2)(len, ptr1, n, ptr2, a, b, ws->data);
    gsl_vector_free(ws);
    return ary;
#endif
  } else {
    return rb_float_new((*f)(NUM2DBL(argv[0]), FIX2INT(argv[1]), NUM2DBL(argv[2]), NUM2DBL(argv[3])));    
  }  

}
static VALUE rb_jac_jacobi_P0_e(VALUE module, VALUE x, VALUE a, VALUE b)
{  
  return jac_eval3_e(x, a, b, jac_jacobi_P0_e);
}

static VALUE rb_jac_jacobi_P0(VALUE module, VALUE x, VALUE a, VALUE b)
{  
  return jac_eval3(x, a, b, jac_jacobi_P0);
}

static VALUE rb_jac_jacobi_P1_e(VALUE module, VALUE x, VALUE a, VALUE b)
{  
  return jac_eval3_e(x, a, b, jac_jacobi_P1_e);
}

static VALUE rb_jac_jacobi_P1(VALUE module, VALUE x, VALUE a, VALUE b)
{  
  return jac_eval3(x, a, b, jac_jacobi_P1);
}


static VALUE rb_jac_jacobi(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_jacobi_eval(argc, argv, jac_jacobi, jac_jacobi_array);
}

static VALUE rb_jac_djacobi_P0_e(VALUE module, VALUE x, VALUE a, VALUE b)
{  
  return jac_eval3_e(x, a, b, jac_djacobi_P0_e);
}

static VALUE rb_jac_djacobi_P0(VALUE module, VALUE x, VALUE a, VALUE b)
{  
  return jac_eval3(x, a, b, jac_djacobi_P0);
}

static VALUE rb_jac_djacobi_P1_e(VALUE module, VALUE x, VALUE a, VALUE b)
{  
  return jac_eval3_e(x, a, b, jac_djacobi_P1_e);
}

static VALUE rb_jac_djacobi_P1(VALUE module, VALUE x, VALUE a, VALUE b)
{  
  return jac_eval3(x, a, b, jac_djacobi_P1);
}

static VALUE rb_jac_djacobi(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_jacobi_eval(argc, argv, jac_djacobi, jac_djacobi_array);  
}

static VALUE rb_jac_zeros_eval(int argc, VALUE *argv, VALUE module,
    int (*f)(double*, int, double, double))
{
  gsl_vector *x;
  int m, status;
  double a, b;
  VALUE xx;
  switch (argc) {
  case 3:
    if (FIXNUM_P(argv[0])) {
      m = FIX2INT(argv[0]);
      a = NUM2DBL(argv[1]);
      b = NUM2DBL(argv[2]);
      x = gsl_vector_alloc(m);
      xx = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, x);
    } else if (VECTOR_P(argv[0])) {
      Data_Get_Struct(argv[0], gsl_vector, x);
      m = x->size;
      a = NUM2DBL(argv[1]);
      b = NUM2DBL(argv[2]);
      xx = argv[0];      
    } else {
      rb_raise(rb_eTypeError, "Wrong argument type %s (Fixnum or GSL::Vector expected)",
              rb_class2name(CLASS_OF(argv[0])));
    }
    break;
  case 4:
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, x);    
    m = FIX2INT(argv[1]);    
    a = NUM2DBL(argv[2]);
    b = NUM2DBL(argv[3]);  
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 3 or 4)", argc);
  }
  status = (*f)(x->data, m, a, b);
  if (status != GSL_SUCCESS)
    rb_raise(rb_eRuntimeError, "Something wrong. (error code %d)", status);
  return xx;
}

static VALUE rb_jac_jacobi_zeros(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_zeros_eval(argc, argv, module, jac_jacobi_zeros);
}

static void jac_define_const(VALUE module)
{
  rb_define_const(module, "GJ", INT2FIX(JAC_GJ));
  rb_define_const(module, "GLJ", INT2FIX(JAC_GLJ));  
  rb_define_const(module, "GRJM", INT2FIX(JAC_GRJM));
  rb_define_const(module, "GRJP", INT2FIX(JAC_GRJP));
}

static VALUE rb_jac_quadrature_alloc(VALUE klass, VALUE vQ)
{
  jac_quadrature *q;
  
  q = jac_quadrature_alloc(FIX2INT(vQ));
  
  return Data_Wrap_Struct(klass, 0, jac_quadrature_free, q);
}

static VALUE rb_jac_quadrature_Q(VALUE obj)
{
  jac_quadrature *q;
  Data_Get_Struct(obj, jac_quadrature, q);
  return INT2FIX(q->Q);
}

static VALUE rb_jac_quadrature_type(VALUE obj)
{
  jac_quadrature *q;
  Data_Get_Struct(obj, jac_quadrature, q);
  return INT2FIX((int) q->type);
}

static VALUE rb_jac_quadrature_alpha(VALUE obj)
{
  jac_quadrature *q;
  Data_Get_Struct(obj, jac_quadrature, q);
  return NUM2DBL(q->alpha);
}

static VALUE rb_jac_quadrature_beta(VALUE obj)
{
  jac_quadrature *q;
  Data_Get_Struct(obj, jac_quadrature, q);
  return NUM2DBL(q->beta);
}

static VALUE rb_jac_quadrature_x(VALUE obj)
{
  jac_quadrature *q;
  gsl_vector_view *v;
  Data_Get_Struct(obj, jac_quadrature, q);
  v = gsl_vector_view_alloc();
  v->vector.data = q->x;
  v->vector.size = q->Q;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_vector_view, 0, gsl_vector_view_free, v);
}

static VALUE rb_jac_quadrature_w(VALUE obj)
{
  jac_quadrature *q;
  gsl_vector_view *v;
  Data_Get_Struct(obj, jac_quadrature, q);
  v = gsl_vector_view_alloc();
  v->vector.data = q->w;
  v->vector.size = q->Q;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_vector_view, 0, gsl_vector_view_free, v);
}

static VALUE rb_jac_quadrature_D(VALUE obj)
{
  jac_quadrature *q;
  gsl_vector_view *v;
  Data_Get_Struct(obj, jac_quadrature, q);
  v = gsl_vector_view_alloc();
  v->vector.data = q->D;
  v->vector.size = q->Q;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_vector_view, 0, gsl_vector_view_free, v);
}

static VALUE rb_jac_quadrature_xp(VALUE obj)
{
  jac_quadrature *q;
  gsl_vector_view *v;
  Data_Get_Struct(obj, jac_quadrature, q);
  v = gsl_vector_view_alloc();
  v->vector.data = q->w;
  v->vector.size = q->np;
  v->vector.stride = 1;
  return Data_Wrap_Struct(cgsl_vector_view, 0, gsl_vector_view_free, v);
}

static VALUE rb_jac_interpmat_alloc(int argc, VALUE *argv, VALUE obj)
{
  int err;
  jac_quadrature *q;
  gsl_vector *xp;
  int np;
  Data_Get_Struct(obj, jac_quadrature, q);  
  switch (argc) {
  case 1:
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, xp);  
    np = xp->size;
    break;
  case 2:
    CHECK_VECTOR(argv[1]);
    Data_Get_Struct(argv[1], gsl_vector, xp);  
    np = FIX2INT(argv[0]);  
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 1 or 2)", argc);
  }
  err = jac_interpmat_alloc(q, np, xp->data);
  return FIX2INT(err);
}

static VALUE rb_jac_interpmat_free(VALUE obj)
{
  jac_quadrature *q;  
  Data_Get_Struct(obj, jac_quadrature, q);    
  jac_interpmat_free(q);
  return Qtrue;
}

static VALUE rb_jac_quadrature_zwd(int argc, VALUE *argv, VALUE obj)
{
  jac_quadrature *q;
  gsl_vector *ws;
  int flag = 0, type, status;
  double a, b;
  Data_Get_Struct(obj, jac_quadrature, q);      
  switch (argc) {
  case 3:
    type = FIX2INT(argv[0]);
    a = NUM2DBL(argv[1]);
    b = NUM2DBL(argv[2]);
    ws = gsl_vector_alloc(q->Q);  
    flag = 1;
    break;
  case 4:
    type = FIX2INT(argv[0]);
    a = NUM2DBL(argv[1]);
    b = NUM2DBL(argv[2]);
    Data_Get_Struct(argv[3], gsl_vector, ws);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 3 or 4)", argc);
  }
  status = jac_quadrature_zwd(q, type, a, b, ws->data);
  if (flag == 1) gsl_vector_free(ws);
  return INT2FIX(status);
}

static VALUE rb_jac_integrate(VALUE obj, VALUE ff)
{
  jac_quadrature *q;
  gsl_vector *f;
  CHECK_VECTOR(ff);
  Data_Get_Struct(obj, jac_quadrature, q);
  Data_Get_Struct(ff, gsl_vector, f);
  return rb_float_new(jac_integrate(q, f->data));
}

static VALUE rb_jac_interpolate(int argc, VALUE *argv, VALUE obj)
{
  jac_quadrature *q;
  gsl_vector *f, *fout;
  VALUE vfout;
  switch (argc) {
  case 1:
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, f);
    fout = gsl_vector_alloc(f->size);
    vfout = Data_Wrap_Struct(VECTOR_ROW_COL(CLASS_OF(argv[0])), 0, gsl_vector_free, fout);
    break;
  case 2:
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, f);
    CHECK_VECTOR(argv[1]);
    Data_Get_Struct(argv[1], gsl_vector, fout);
    vfout = argv[1];
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 1 or 2)", argc);
  }
  Data_Get_Struct(obj, jac_quadrature, q);
  jac_interpolate(q, f->data, fout->data);
  return vfout;
}

static VALUE rb_jac_differentiate(int argc, VALUE *argv, VALUE obj)
{
  jac_quadrature *q;
  gsl_vector *f, *fout;
  VALUE vfout;
  switch (argc) {
  case 1:
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, f);
    fout = gsl_vector_alloc(f->size);
    vfout = Data_Wrap_Struct(VECTOR_ROW_COL(CLASS_OF(argv[0])), 0, gsl_vector_free, fout);
    break;
  case 2:
    CHECK_VECTOR(argv[0]);
    Data_Get_Struct(argv[0], gsl_vector, f);
    CHECK_VECTOR(argv[1]);
    Data_Get_Struct(argv[1], gsl_vector, fout);
    vfout = argv[1];
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 1 or 2)", argc);
  }
  Data_Get_Struct(obj, jac_quadrature, q);
  jac_differentiate(q, f->data, fout->data);
  return vfout;
}

/*****/
static VALUE rb_jac_qeval(int argc, VALUE *argv, 
    int (*f)(double*, double*, const int, double, double, double*))
{
  gsl_vector *z, *D, *ws;
  int Q;
  int flag = 0, status;
  double alpha, beta;
  VALUE vD;
  if (argc < 3) rb_raise(rb_eArgError, "Too few arguments (%d for >= 3)", argc);
  CHECK_VECTOR(argv[0]);
  Data_Get_Struct(argv[0], gsl_vector, z);
  argc -= 1;
  argv += 1;
  if (VECTOR_P(argv[argc-1])) {
    Data_Get_Struct(argv[argc-1], gsl_vector, ws);
    argc -= 1;
  } else {
    ws = gsl_vector_alloc(z->size);
    flag = 1;
  }  
  switch (argc) {
  case 2:
    Q = z->size;
    D = gsl_vector_alloc(Q*Q);
    alpha = NUM2DBL(argv[0]);
    beta = NUM2DBL(argv[1]);
    vD = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, D);
    break;
  case 4:
    Data_Get_Struct(argv[0], gsl_vector, D);
    vD = argv[0];
    Q = FIX2INT(argv[1]);
    alpha = NUM2DBL(argv[2]);
    beta = NUM2DBL(argv[3]);
    break;
  case 3:
    if (VECTOR_P(argv[0])) {
      Q = z->size;
      Data_Get_Struct(argv[0], gsl_vector, D);
      vD = argv[0];
    } else {
      Q = FIX2INT(argv[0]);      
      D = gsl_vector_alloc(Q*Q);
      vD = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, D);
    }
    alpha = NUM2DBL(argv[1]);
     beta = NUM2DBL(argv[2]);
     break;
   default:
     rb_raise(rb_eArgError, "Wrong number of arguments.");
  }
  status = (*f)(z->data, D->data, Q, alpha, beta, ws->data);
  if (flag == 1) gsl_vector_free(ws);
  if (status != GSL_SUCCESS) rb_raise(rb_eRuntimeError, "Something wrong.");
  return vD;
}

static VALUE rb_jac_diffmat_gj(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_qeval(argc, argv, jac_diffmat_gj);
}

static VALUE rb_jac_diffmat_glj(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_qeval(argc, argv, jac_diffmat_glj);
}

static VALUE rb_jac_diffmat_grjm(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_qeval(argc, argv, jac_diffmat_grjm);
}

static VALUE rb_jac_diffmat_grjp(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_qeval(argc, argv, jac_diffmat_grjp);
}

static VALUE rb_jac_weights_gj(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_qeval(argc, argv, jac_weights_gj);
}

static VALUE rb_jac_weights_glj(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_qeval(argc, argv, jac_weights_glj);
}

static VALUE rb_jac_weights_grjm(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_qeval(argc, argv, jac_weights_grjm);
}

static VALUE rb_jac_weights_grjp(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_qeval(argc, argv, jac_weights_grjp);
}

static VALUE rb_jac_zeros_gj(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_zeros_eval(argc, argv, module, jac_zeros_gj);
}
static VALUE rb_jac_zeros_glj(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_zeros_eval(argc, argv, module, jac_zeros_glj);
}
static VALUE rb_jac_zeros_grjm(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_zeros_eval(argc, argv, module, jac_zeros_grjm);
}
static VALUE rb_jac_zeros_grjp(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_zeros_eval(argc, argv, module, jac_zeros_grjp);
}

static VALUE rb_jac_lagrange_eval(int argc, VALUE *argv,
    double (*f)(int, double, int, double*, double, double))
{
  gsl_vector *z;
  int i, Q;
  double alpha, beta, zz;

  switch (argc) {
  case 5:
    i = FIX2INT(argv[0]);
    zz = NUM2DBL(argv[1]);
    CHECK_VECTOR(argv[2]);
    Data_Get_Struct(argv[2], gsl_vector, z);
    Q = z->size;
    alpha = NUM2DBL(argv[3]);
    beta = NUM2DBL(argv[4]);
  break;
  case 6:
    i = FIX2INT(argv[0]);
    zz = NUM2DBL(argv[1]);
    Q = FIX2INT(argv[2]);
    CHECK_VECTOR(argv[3]);
    Data_Get_Struct(argv[3], gsl_vector, z);
    alpha = NUM2DBL(argv[4]);
    beta = NUM2DBL(argv[5]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 5 or 6)", argc);
  }
  return rb_float_new((*f)(i, zz, Q, z->data, alpha, beta));
}

static VALUE rb_jac_lagrange_gj(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_lagrange_eval(argc, argv, jac_lagrange_gj);
}
static VALUE rb_jac_lagrange_glj(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_lagrange_eval(argc, argv, jac_lagrange_glj);
}
static VALUE rb_jac_lagrange_grjm(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_lagrange_eval(argc, argv, jac_lagrange_grjm);
}

static VALUE rb_jac_lagrange_grjp(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_lagrange_eval(argc, argv, jac_lagrange_grjp);
}

static VALUE rb_jac_interpmat_eval(int argc, VALUE *argv,
  int (*f) (double*, double*, int, double*, int, double, double))
{
  gsl_vector *imat, *zp, *z;
  double alpha, beta;
  int np, Q, status;
  VALUE vimat;
  
  if (argc < 3) rb_raise(rb_eArgError, "Too few arguments (%d for >= 3)", argc);
  
  CHECK_VECTOR(argv[0]);
  if (VECTOR_P(argv[1])) {
    Data_Get_Struct(argv[0], gsl_vector, imat);
    Data_Get_Struct(argv[1], gsl_vector, zp);    
    vimat = argv[0];
    if (FIXNUM_P(argv[2])) np = FIX2INT(argv[2]);
    argc -= 3;
    argv += 3;
  } else {
    Data_Get_Struct(argv[0], gsl_vector, zp);
    if (FIXNUM_P(argv[1])) np = FIX2INT(argv[1]);    
    else np = zp->size;
    imat = gsl_vector_alloc(np);
    vimat = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, imat);
    argc -= 2;
    argv += 2;    
  }
  CHECK_VECTOR(argv[0]);
  Data_Get_Struct(argv[0], gsl_vector, z);
  argc -= 1;
  argv += 1;
  switch (argc) {
  case 3:
    Q = FIX2INT(argv[0]);
    alpha = NUM2DBL(argv[1]);
    beta = NUM2DBL(argv[2]);
    break;
  case 2:
    Q = z->size;
    alpha = NUM2DBL(argv[0]);
    beta = NUM2DBL(argv[1]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments.");
  }
  (*f)(imat->data, zp->data, np, z->data, Q, alpha, beta);
  if (status != GSL_SUCCESS) rb_raise(rb_eRuntimeError, "Some error.");
  return vimat;
}

static VALUE rb_jac_interpmat_gj(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_interpmat_eval(argc, argv, jac_interpmat_gj);
}

static VALUE rb_jac_interpmat_glj(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_interpmat_eval(argc, argv, jac_interpmat_glj);
}
static VALUE rb_jac_interpmat_grjm(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_interpmat_eval(argc, argv, jac_interpmat_grjm);
}
static VALUE rb_jac_interpmat_grjp(int argc, VALUE *argv, VALUE module)
{
  return rb_jac_interpmat_eval(argc, argv, jac_interpmat_grjp);
}
void Init_jacobi(VALUE module)
{
  VALUE mjac, cjacq;
  
  mjac = rb_define_module("Jac");
  jac_define_const(mjac);
  cjacq = rb_define_class_under(mjac, "Quadrature", cGSL_Object);
  
  rb_define_module_function(mjac, "jacobi_P0_e", rb_jac_jacobi_P0_e, 3);
  rb_define_module_function(mjac, "jacobi_P0", rb_jac_jacobi_P0, 3);  
  rb_define_module_function(mjac, "jacobi_P1_e", rb_jac_jacobi_P1_e, 3);
  rb_define_module_function(mjac, "jacobi_P1", rb_jac_jacobi_P1, 3);    
  rb_define_module_function(mjac, "jacobi", rb_jac_jacobi, -1);      
  rb_define_module_function(mjac, "djacobi_P0_e", rb_jac_djacobi_P0_e, 3);
  rb_define_module_function(mjac, "djacobi_P0", rb_jac_djacobi_P0, 3);  
  rb_define_module_function(mjac, "djacobi_P1_e", rb_jac_djacobi_P1_e, 3);
  rb_define_module_function(mjac, "djacobi_P1", rb_jac_djacobi_P1, 3);      
  rb_define_module_function(mjac, "djacobi", rb_jac_djacobi, 4);        
  
  rb_define_module_function(mjac, "jacobi_zeros", rb_jac_jacobi_zeros, -1);          

  /*****/
  rb_define_singleton_method(cjacq, "alloc", rb_jac_quadrature_alloc, 1);
  rb_define_method(cjacq, "Q", rb_jac_quadrature_Q, 0);
  rb_define_method(cjacq, "type", rb_jac_quadrature_type, 0);  
  rb_define_method(cjacq, "alpha", rb_jac_quadrature_alpha, 0);  
  rb_define_method(cjacq, "beta", rb_jac_quadrature_beta, 0);      
  rb_define_method(cjacq, "x", rb_jac_quadrature_x, 0);
  rb_define_method(cjacq, "w", rb_jac_quadrature_w, 0);
  rb_define_method(cjacq, "D", rb_jac_quadrature_D, 0);
  rb_define_method(cjacq, "xp", rb_jac_quadrature_xp, 0);    
  rb_define_method(cjacq, "interpmat_alloc", rb_jac_interpmat_alloc, -1);        
  rb_define_method(cjacq, "interpmat_free", rb_jac_interpmat_free, 0);
  rb_define_method(cjacq, "zwd", rb_jac_quadrature_zwd, -1);
  rb_define_method(cjacq, "integrate", rb_jac_integrate, 1);
  rb_define_method(cjacq, "interpolate", rb_jac_interpolate, -1);
  rb_define_method(cjacq, "differentiate", rb_jac_differentiate, -1);  
  /*****/
  rb_define_module_function(mjac, "diffmat_gj", rb_jac_diffmat_gj, -1);
  rb_define_module_function(mjac, "diffmat_glj", rb_jac_diffmat_glj, -1);  
  rb_define_module_function(mjac, "diffmat_grjm", rb_jac_diffmat_grjm, -1);  
  rb_define_module_function(mjac, "diffmat_grjp", rb_jac_diffmat_grjp, -1);
  
  rb_define_module_function(mjac, "weights_gj", rb_jac_weights_gj, -1);  
  rb_define_module_function(mjac, "weights_glj", rb_jac_weights_glj, -1);  
  rb_define_module_function(mjac, "weights_grjm", rb_jac_weights_grjm, -1);  
  rb_define_module_function(mjac, "weights_grjp", rb_jac_weights_grjp, -1);  

  rb_define_module_function(mjac, "zeros_gj", rb_jac_zeros_gj, -1);  
  rb_define_module_function(mjac, "zeros_glj", rb_jac_zeros_glj, -1);  
  rb_define_module_function(mjac, "zeros_grjm", rb_jac_zeros_grjm, -1);  
  rb_define_module_function(mjac, "zeros_grjp", rb_jac_zeros_grjp, -1);  
  
  rb_define_module_function(mjac, "lagrange_gj", rb_jac_lagrange_gj, -1);    
  rb_define_module_function(mjac, "lagrange_glj", rb_jac_lagrange_glj, -1);  
  rb_define_module_function(mjac, "lagrange_grjm", rb_jac_lagrange_grjm, -1);  
  rb_define_module_function(mjac, "lagrange_grjp", rb_jac_lagrange_grjp, -1);  

  rb_define_module_function(mjac, "interpmat_gj", rb_jac_interpmat_gj, -1);  
  rb_define_module_function(mjac, "interpmat_glj", rb_jac_interpmat_glj, -1);  
  rb_define_module_function(mjac, "interpmat_grjm", rb_jac_interpmat_grjm, -1);  
  rb_define_module_function(mjac, "interpmat_grjp", rb_jac_interpmat_grjp, -1);  
  
}

#endif

