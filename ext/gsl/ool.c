#ifdef HAVE_OOL_OOL_VERSION_H
#include "include/rb_gsl.h"
#include "include/rb_gsl_array.h"
#include <ool/ool_conmin.h>

static VALUE cool_conmin_function;
static VALUE cool_conmin_constraint;
static VALUE cool_conmin_pgrad;
static VALUE cool_conmin_spg;
static VALUE cool_conmin_gencan;
static VALUE cool_conmin_pgrad_parameters;
static VALUE cool_conmin_spg_parameters;
static VALUE cool_conmin_gencan_parameters;

#ifndef CHECK_OOL_CONMIN_FUNCTION
#define CHECK_OOL_CONMIN_FUNCTION(x) if(CLASS_OF(x)!=cool_conmin_function)\
      rb_raise(rb_eTypeError,\
      "wrong argument type %s (OOL::Conmin::Function expected)",\
      rb_class2name(CLASS_OF(x)));
#endif

enum enum_ool_conmin_minimizer_type {
  OOL_CONMIN_PGRAD,
  OOL_CONMIN_SPG,
  OOL_CONMIN_GENCAN,
};

static const ool_conmin_minimizer_type* get_minimizer_type(VALUE t)
{
  char name[64];

  switch (TYPE(t)) {
  case T_STRING:
    strcpy(name, STR2CSTR(t));
    if (str_tail_grep(name, "pgrad") == 0) {
      return ool_conmin_minimizer_pgrad;
    } else if (str_tail_grep(name, "spg") == 0) {
      return ool_conmin_minimizer_spg;
    } else if (str_tail_grep(name, "gencan") == 0) {
      return ool_conmin_minimizer_gencan;
    } else {
      rb_raise(rb_eTypeError, "%s: unknown minimizer type", name);
    }
    break;
  case T_FIXNUM:
    switch (FIX2INT(t)) {
      case OOL_CONMIN_PGRAD:
        return ool_conmin_minimizer_pgrad;
        break;
      case OOL_CONMIN_SPG:
        return ool_conmin_minimizer_spg;
        break;
      case OOL_CONMIN_GENCAN:
        return ool_conmin_minimizer_gencan;
        break;
      default:
         rb_raise(rb_eTypeError, "%d: unknown minimizer type", FIX2INT(t));
        break;
    }
    break;
  default:
    if (t == cool_conmin_pgrad) return ool_conmin_minimizer_pgrad;
    else if (t == cool_conmin_spg) return ool_conmin_minimizer_spg;
    else if (t == cool_conmin_gencan) return ool_conmin_minimizer_gencan;
    else rb_raise(rb_eTypeError, "type is given by a String or a Fixnum");
    break;
  }
}

static void def_const(VALUE module)
{
  rb_define_const(module, "CONTINUE", INT2FIX(OOL_CONTINUE));
  rb_define_const(module, "SUCCESS", INT2FIX(OOL_SUCCESS));
}


static VALUE rb_ool_conmin_minimizer_set(int argc, VALUE *argv, VALUE obj);
static VALUE rb_ool_conmin_minimizer_alloc(int argc, VALUE *argv, VALUE klass)
{
  ool_conmin_minimizer *m;
  VALUE obj;
  if (argc < 2) rb_raise(rb_eArgError, "Too few arguments (%d for >= 2)", argc);
  m = ool_conmin_minimizer_alloc(get_minimizer_type(argv[0]), FIX2INT(argv[1]));
  obj = Data_Wrap_Struct(klass, 0, ool_conmin_minimizer_free, m);
  if (argc > 2) rb_ool_conmin_minimizer_set(argc-2, argv+2, obj);
  return obj;
}

static void* get_parameter(const ool_conmin_minimizer_type *T, ool_conmin_pgrad_parameters *Pp,
  ool_conmin_spg_parameters *Ps, ool_conmin_gencan_parameters *Pg, VALUE ary);
static VALUE rb_ool_conmin_minimizer_set(int argc, VALUE *argv, VALUE obj)
{
  ool_conmin_minimizer *m;
  ool_conmin_function *F;
  ool_conmin_constraint *C;
  gsl_vector *v;
  ool_conmin_pgrad_parameters Pp;
  ool_conmin_spg_parameters Ps;
  ool_conmin_gencan_parameters Pg;
  void *P;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  switch (argc) {
  case 3:
    if (CLASS_OF(argv[0]) != cool_conmin_function)
      rb_raise(rb_eTypeError, "Wrong argument type 0 (OOL::Conmin::Function expected)");
    if (CLASS_OF(argv[1]) != cool_conmin_constraint)
      rb_raise(rb_eTypeError, "Wrong argument type 1 (OOL::Conmin::Constraint expected)");
    if (!VECTOR_P(argv[2]))
      rb_raise(rb_eTypeError, "Wrong argument type 2 (GSL::Vector expected)");
    Data_Get_Struct(argv[0], ool_conmin_function, F);
    Data_Get_Struct(argv[1], ool_conmin_constraint, C);
    Data_Get_Struct(argv[2], gsl_vector, v);
    P = get_parameter(m->type, &Pp, &Ps, &Pg, Qnil);
    ool_conmin_minimizer_set(m, F, C, v, P);
    break;
  case 4:
    if (CLASS_OF(argv[0]) != cool_conmin_function)
      rb_raise(rb_eTypeError, "Wrong argument type 0 (OOL::Conmin::Function expected)");
    if (CLASS_OF(argv[1]) != cool_conmin_constraint)
      rb_raise(rb_eTypeError, "Wrong argument type 1 (OOL::Conmin::Constraint expected)");
    if (!VECTOR_P(argv[2]))
      rb_raise(rb_eTypeError, "Wrong argument type 2 (GSL::Vector expected)");
    if (!rb_obj_is_kind_of(argv[3], rb_cArray) && argv[3] != Qnil)
      rb_raise(rb_eTypeError, "Wrong argument type 3 (Array expected)");
    Data_Get_Struct(argv[0], ool_conmin_function, F);
    Data_Get_Struct(argv[1], ool_conmin_constraint, C);
    Data_Get_Struct(argv[2], gsl_vector, v);
    P = get_parameter(m->type, &Pp, &Ps, &Pg, argv[3]);
    ool_conmin_minimizer_set(m, F, C, v, P);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 3 or 4)", argc);
  }
  return obj;
}

static void* get_parameter(const ool_conmin_minimizer_type *T, ool_conmin_pgrad_parameters *Pp,
  ool_conmin_spg_parameters *Ps, ool_conmin_gencan_parameters *Pg, VALUE ary)
{
  if (T == ool_conmin_minimizer_pgrad) {
    if (ary == Qnil) {
      ool_conmin_parameters_default(T, (void*) Pp);
    } else {
      Pp->fmin = NUM2DBL(rb_ary_entry(ary, 0));
      Pp->tol = NUM2DBL(rb_ary_entry(ary, 1));
      Pp->alpha = NUM2DBL(rb_ary_entry(ary, 2));
      Pp->sigma1 = NUM2DBL(rb_ary_entry(ary, 3));
      Pp->sigma2 = NUM2DBL(rb_ary_entry(ary, 4));
    }
    return (void*) Pp;
  } else if (T == ool_conmin_minimizer_spg) {
    if (ary == Qnil) {
      ool_conmin_parameters_default(T, (void*) Ps);
    } else {
      Ps->fmin = NUM2DBL(rb_ary_entry(ary, 0));
      Ps->tol = NUM2DBL(rb_ary_entry(ary, 1));
      Ps->M = NUM2DBL(rb_ary_entry(ary, 2));
      Ps->alphamin = NUM2DBL(rb_ary_entry(ary, 3));
      Ps->alphamax = NUM2DBL(rb_ary_entry(ary, 4));
      Ps->gamma = NUM2DBL(rb_ary_entry(ary, 5));
      Ps->sigma2 = NUM2DBL(rb_ary_entry(ary, 6));
      Ps->sigma2 = NUM2DBL(rb_ary_entry(ary, 7));
    }
    return (void*) Ps;
  } else {
    if (ary == Qnil) {
      ool_conmin_parameters_default(T, (void*) Pg);
    } else {
      Pg->epsgpen = NUM2DBL(rb_ary_entry(ary, 0));
      Pg->epsgpsn = NUM2DBL(rb_ary_entry(ary, 1));
      Pg->fmin = NUM2DBL(rb_ary_entry(ary, 2));
      Pg->udelta0 = NUM2DBL(rb_ary_entry(ary, 3));
      Pg->ucgmia = NUM2DBL(rb_ary_entry(ary, 4));
      Pg->ucgmib = NUM2DBL(rb_ary_entry(ary, 5));
      Pg->cg_scre = FIX2INT(rb_ary_entry(ary, 6));
      Pg->cg_gpnf = NUM2DBL(rb_ary_entry(ary, 7));
      Pg->cg_epsi = NUM2DBL(rb_ary_entry(ary, 8));
      Pg->cg_epsf = NUM2DBL(rb_ary_entry(ary, 9));
      Pg->cg_epsnqmp = NUM2DBL(rb_ary_entry(ary, 10));
      Pg->cg_maxitnqmp = (size_t) FIX2INT(rb_ary_entry(ary, 11));
      Pg->nearlyq = FIX2INT(rb_ary_entry(ary, 12));
      Pg->nint = NUM2DBL(rb_ary_entry(ary, 13));
      Pg->next = NUM2DBL(rb_ary_entry(ary, 14));
      Pg->mininterp = (size_t) FIX2INT(rb_ary_entry(ary, 15));
      Pg->maxextrap = (size_t) FIX2INT(rb_ary_entry(ary, 16));
      Pg->trtype = FIX2INT(rb_ary_entry(ary, 17));
      Pg->eta = NUM2DBL(rb_ary_entry(ary, 18));
      Pg->delmin = NUM2DBL(rb_ary_entry(ary, 19));
      Pg->lspgmi = NUM2DBL(rb_ary_entry(ary, 20));
      Pg->lspgma = NUM2DBL(rb_ary_entry(ary, 21));
      Pg->theta = NUM2DBL(rb_ary_entry(ary, 22));
      Pg->gamma = NUM2DBL(rb_ary_entry(ary, 23));
      Pg->beta = NUM2DBL(rb_ary_entry(ary, 24));
      Pg->sigma1 = NUM2DBL(rb_ary_entry(ary, 25));
      Pg->sigma2 = NUM2DBL(rb_ary_entry(ary, 26));
      Pg->epsrel = NUM2DBL(rb_ary_entry(ary, 27));
      Pg->epsabs = NUM2DBL(rb_ary_entry(ary, 28));
      Pg->infrel = NUM2DBL(rb_ary_entry(ary, 29));
      Pg->infabs = NUM2DBL(rb_ary_entry(ary, 30));
    }
    return (void*) Pg;
  }
}

static VALUE create_parameters_ary_pgrad(ool_conmin_pgrad_parameters *Pp)
{
  VALUE ary;
  ary = rb_ary_new2(5);
  rb_ary_store(ary, 0, rb_float_new(Pp->fmin));
  rb_ary_store(ary, 1, rb_float_new(Pp->tol));
  rb_ary_store(ary, 2, rb_float_new(Pp->alpha));
  rb_ary_store(ary, 3, rb_float_new(Pp->sigma1));
  rb_ary_store(ary, 4, rb_float_new(Pp->sigma2));
  return ary;
}

static VALUE create_parameters_ary_spg(ool_conmin_spg_parameters *Ps)
{
  VALUE ary;
    ary = rb_ary_new2(8);
    rb_ary_store(ary, 0, rb_float_new(Ps->fmin));
    rb_ary_store(ary, 1, rb_float_new(Ps->tol));
    rb_ary_store(ary, 2, rb_float_new(Ps->M));
    rb_ary_store(ary, 3, rb_float_new(Ps->alphamin));
    rb_ary_store(ary, 4, rb_float_new(Ps->alphamax));
    rb_ary_store(ary, 5, rb_float_new(Ps->gamma));
    rb_ary_store(ary, 6, rb_float_new(Ps->sigma2));
    rb_ary_store(ary, 7, rb_float_new(Ps->sigma2));
    return ary;
}

static VALUE create_parameters_ary_gencan(ool_conmin_gencan_parameters *Pg)
{
  VALUE ary;
    ary = rb_ary_new2(31);
    rb_ary_store(ary, 0, rb_float_new(Pg->epsgpen));
    rb_ary_store(ary, 1, rb_float_new(Pg->epsgpsn));
    rb_ary_store(ary, 2, rb_float_new(Pg->fmin));
    rb_ary_store(ary, 3, rb_float_new(Pg->udelta0));
    rb_ary_store(ary, 4, rb_float_new(Pg->ucgmia));
    rb_ary_store(ary, 5, rb_float_new(Pg->ucgmib));
    rb_ary_store(ary, 6, INT2FIX(Pg->cg_scre));
    rb_ary_store(ary, 7, rb_float_new(Pg->cg_gpnf));
    rb_ary_store(ary, 8, rb_float_new(Pg->cg_epsi));
    rb_ary_store(ary, 9, rb_float_new(Pg->cg_epsf));
    rb_ary_store(ary, 10, rb_float_new(Pg->cg_epsnqmp));
    rb_ary_store(ary, 11, INT2FIX((int) Pg->cg_maxitnqmp));
    rb_ary_store(ary, 12, INT2FIX(Pg->nearlyq));
    rb_ary_store(ary, 13, rb_float_new(Pg->nint));
    rb_ary_store(ary, 14, rb_float_new(Pg->next));
    rb_ary_store(ary, 15, INT2FIX((int)Pg->mininterp));
    rb_ary_store(ary, 16, INT2FIX((int)Pg->maxextrap));
    rb_ary_store(ary, 17, INT2FIX(Pg->trtype));
    rb_ary_store(ary, 18, rb_float_new(Pg->eta));
    rb_ary_store(ary, 19, rb_float_new(Pg->delmin));
    rb_ary_store(ary, 20, rb_float_new(Pg->lspgmi));
    rb_ary_store(ary, 21, rb_float_new(Pg->lspgma));
    rb_ary_store(ary, 22, rb_float_new(Pg->theta));
    rb_ary_store(ary, 23, rb_float_new(Pg->gamma));
    rb_ary_store(ary, 24, rb_float_new(Pg->beta));
    rb_ary_store(ary, 25, rb_float_new(Pg->sigma1));
    rb_ary_store(ary, 26, rb_float_new(Pg->sigma2));
    rb_ary_store(ary, 27, rb_float_new(Pg->epsrel));
    rb_ary_store(ary, 28, rb_float_new(Pg->epsabs));
    rb_ary_store(ary, 29, rb_float_new(Pg->infrel));
    rb_ary_store(ary, 30, rb_float_new(Pg->infabs));
  return ary;
}
static VALUE rb_ool_conmin_minimizer_parameters_get(VALUE obj)
{
  ool_conmin_minimizer *m;
  ool_conmin_pgrad_parameters *Pp;
  ool_conmin_spg_parameters *Ps;
  ool_conmin_gencan_parameters *Pg;
  void *P;
  VALUE ary;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  ool_conmin_parameters_get(m, P);
  if (m->type ==   ool_conmin_minimizer_pgrad) {
    Pp = (ool_conmin_pgrad_parameters*) P;
    ary = create_parameters_ary_pgrad(Pp);
  } else if (m->type == ool_conmin_minimizer_spg) {
    Ps = (ool_conmin_spg_parameters*) P;
    ary = create_parameters_ary_spg(Ps);
  } else {
    Pg = (ool_conmin_gencan_parameters*) P;
    ary = create_parameters_ary_gencan(Pg);
  }
  return ary;
}

static VALUE rb_ool_conmin_minimizer_parameters_set(VALUE obj, VALUE params)
{
  ool_conmin_minimizer *m;
  ool_conmin_pgrad_parameters *Pp;
  ool_conmin_spg_parameters *Ps;
  ool_conmin_gencan_parameters *Pg;
  void *P;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  P = get_parameter(m->type, Pp, Ps, Pg, params);
  ool_conmin_parameters_set(m, P);
  return params;
}

static VALUE rb_ool_conmin_minimizer_name(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return rb_str_new2(ool_conmin_minimizer_name(m));
}
static VALUE rb_ool_conmin_minimizer_f(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return rb_float_new(m->f);
}
static VALUE rb_ool_conmin_minimizer_x(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return Data_Wrap_Struct(cgsl_vector, 0, NULL, m->x);
}
static VALUE rb_ool_conmin_minimizer_gradient(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return Data_Wrap_Struct(cgsl_vector, 0, NULL, m->gradient);
}
static VALUE rb_ool_conmin_minimizer_minimum(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return rb_float_new(ool_conmin_minimizer_minimum(m));
}
static VALUE rb_ool_conmin_minimizer_dx(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return Data_Wrap_Struct(cgsl_vector, 0, NULL, m->dx);
}
static VALUE rb_ool_conmin_minimizer_size(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return rb_float_new(ool_conmin_minimizer_size(m));
}
static VALUE rb_ool_conmin_minimizer_fcount(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return INT2FIX((int) ool_conmin_minimizer_fcount(m));
}
static VALUE rb_ool_conmin_minimizer_gcount(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return INT2FIX((int) ool_conmin_minimizer_gcount(m));
}
static VALUE rb_ool_conmin_minimizer_hcount(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return INT2FIX((int) ool_conmin_minimizer_hcount(m));
}
static VALUE rb_ool_conmin_minimizer_is_optimal(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return INT2FIX((int) ool_conmin_is_optimal(m));
}
static VALUE rb_ool_conmin_minimizer_is_optimal2(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  if (ool_conmin_is_optimal(m)) return Qtrue;
  else return Qfalse;
}
static VALUE rb_ool_conmin_minimizer_iterate(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return INT2FIX((int) ool_conmin_minimizer_iterate(m));
}
static VALUE rb_ool_conmin_minimizer_restart(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  return INT2FIX((int) ool_conmin_minimizer_restart(m));
}

static VALUE rb_ool_conmin_pgrad_parameters_default(VALUE klass);
static VALUE rb_ool_conmin_spg_parameters_default(VALUE klass);
static VALUE rb_ool_conmin_gencan_parameters_default(VALUE klass);

static VALUE rb_ool_conmin_minimizer_parameters_default(VALUE obj)
{
  ool_conmin_minimizer *m;
  Data_Get_Struct(obj, ool_conmin_minimizer, m);
  if (m->type == ool_conmin_minimizer_spg) {
      return rb_ool_conmin_spg_parameters_default(cool_conmin_spg);
  } else if (m->type == ool_conmin_minimizer_pgrad) {
      return rb_ool_conmin_pgrad_parameters_default(cool_conmin_pgrad);
  } else if (m->type == ool_conmin_minimizer_gencan) {
      return rb_ool_conmin_gencan_parameters_default(cool_conmin_gencan);
  } else {
      rb_raise(rb_eRuntimeError, "Unkowm minimizer type.");
  }
  return Qnil;   /* never reaches here */
}

/***/

static void rb_ool_conmin_function_mark(ool_conmin_function *F)
{
  rb_gc_mark((VALUE) F->params);
}

static double rb_ool_conmin_function_f(const gsl_vector *x, void *p);
static void rb_ool_conmin_function_df(const gsl_vector *x, void *p, gsl_vector *g);
static void rb_ool_conmin_function_fdf(const gsl_vector *x, void *p,
              double *f, gsl_vector *g);
static void rb_ool_conmin_function_Hv(const gsl_vector *X, void *params,
      const gsl_vector *V, gsl_vector *hv);
static void set_functions(int argc, VALUE *argv, ool_conmin_function *F);
static void set_params(ool_conmin_function *F, VALUE p);
static VALUE rb_ool_conmin_function_set_n(VALUE obj, VALUE nn);
static VALUE rb_ool_conmin_function_set_functions(int argc, VALUE *argv, VALUE obj);
static VALUE rb_ool_conmin_function_set(int argc, VALUE *argv, VALUE obj);

static VALUE rb_ool_conmin_function_alloc(int argc, VALUE *argv, VALUE klass)
{
  ool_conmin_function *F = NULL;
  VALUE ary, obj;

  F = ALLOC(ool_conmin_function);
  F->f = &rb_ool_conmin_function_f;
  F->df = &rb_ool_conmin_function_df;
  F->fdf = &rb_ool_conmin_function_fdf;
  F->Hv = &rb_ool_conmin_function_Hv;
  F->n = 0;
  ary = rb_ary_new2(5);

  F->params = (void *) ary;
  rb_ary_store(ary, 0, Qnil);  /* proc f */
  rb_ary_store(ary, 1, Qnil);  /* proc df */
  rb_ary_store(ary, 2, Qnil);  /* proc fdf */
  rb_ary_store(ary, 3, Qnil);  /* proc Hv */
  rb_ary_store(ary, 4, Qnil);  /* params */
//  set_functions(argc, argv, F);
  obj = Data_Wrap_Struct(klass, rb_ool_conmin_function_mark, free, F);
  rb_ool_conmin_function_set(argc, argv, obj);
  return obj;
}

static VALUE rb_ool_conmin_function_set(int argc, VALUE *argv, VALUE obj)
{
  ool_conmin_function *F;
  Data_Get_Struct(obj, ool_conmin_function, F);
  switch (argc) {
  case 0:
      break;
  case 1:
    rb_ool_conmin_function_set_n(obj, argv[0]);
    break;
  case 4:
    set_functions(argc, argv, F);
    break;
  case 5:
    if (FIXNUM_P(argv[0])) {
      rb_ool_conmin_function_set_n(obj, argv[0]);
      set_functions(argc-1, argv+1, F);
    } else {
      set_functions(argc-1, argv, F);
      set_params(F, argv[argc-1]);
    }
    break;
  case 6:
    rb_ool_conmin_function_set_n(obj, argv[0]);
    set_functions(argc-2, argv+1, F);
    set_params(F, argv[argc-1]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments.");
  }
  return obj;
}

static VALUE rb_ool_conmin_function_set_n(VALUE obj, VALUE nn)
{
  ool_conmin_function *F = NULL;
  if (FIXNUM_P(nn)) {
    Data_Get_Struct(obj, ool_conmin_function, F);
    F->n = (size_t) FIX2INT(nn);
  } else {
      rb_raise(rb_eArgError, "Wrong argument type %s (Fixnum expected)",
              rb_class2name(CLASS_OF(nn)));
  }
  return nn;
}

static VALUE rb_ool_conmin_function_n(VALUE obj)
{
  ool_conmin_function *F = NULL;
  Data_Get_Struct(obj, ool_conmin_function, F);
  return INT2FIX((int) F->n);
}
static double rb_ool_conmin_function_f(const gsl_vector *x, void *p)
{
  VALUE vx, proc, vp, result, ary;
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  ary = (VALUE) p;
  proc = rb_ary_entry(ary, 0);
  vp = rb_ary_entry(ary, RARRAY_LEN(ary)-1);
  if (NIL_P(vp)) result = rb_funcall(proc, RBGSL_ID_call, 1, vx);
  else result = rb_funcall(proc, RBGSL_ID_call, 2, vx, vp);
  return NUM2DBL(result);
}

static void rb_ool_conmin_function_df(const gsl_vector *x, void *p, gsl_vector *g)
{
  VALUE vx, vg, proc, vp, ary;
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  vg = Data_Wrap_Struct(cgsl_vector, 0, NULL, g);
  ary = (VALUE) p;
  proc = rb_ary_entry(ary, 1);
  vp = rb_ary_entry(ary, RARRAY_LEN(ary)-1);
  if (NIL_P(vp)) {
    rb_funcall(proc, RBGSL_ID_call, 2, vx, vg);
  } else {
    rb_funcall(proc, RBGSL_ID_call, 3, vx, vp, vg);
  }
}

static void rb_ool_conmin_function_fdf(const gsl_vector *x, void *p,
              double *f, gsl_vector *g)
{
  VALUE vx, vf, vg, proc_fdf, proc_f, proc_df, vp, ary, result;
  vx = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector *) x);
  vg = Data_Wrap_Struct(cgsl_vector, 0, NULL, g);
  vf = rb_float_new(*f);
  ary = (VALUE) p;
  proc_f = rb_ary_entry(ary, 0);
  proc_df = rb_ary_entry(ary, 1);
  proc_fdf = rb_ary_entry(ary, 2);
  vp = rb_ary_entry(ary, RARRAY_LEN(ary)-1);
  if (NIL_P(vp)) {
    result = rb_funcall(proc_f, RBGSL_ID_call, 1, vx);
    rb_funcall(proc_df, RBGSL_ID_call, 2, vx, vg);
  } else {
    result = rb_funcall(proc_f, RBGSL_ID_call, 2, vx, vp);
    rb_funcall(proc_df, RBGSL_ID_call, 3, vx, vp, vg);
  }
  *f = NUM2DBL(result);
}

static void rb_ool_conmin_function_Hv(const gsl_vector *X, void *params,
      const gsl_vector *V, gsl_vector *hv)
{
  VALUE vX, vV, vHv, ary, proc_Hv, vp;
  vX = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector*) X);
  vV = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector*) V);
  vHv = Data_Wrap_Struct(cgsl_vector, 0, NULL, (gsl_vector*) hv);
  ary = (VALUE) params;
  proc_Hv = rb_ary_entry(ary, 3);
  vp = rb_ary_entry(ary, RARRAY_LEN(ary)-1);
  if (NIL_P(vp)) {
    rb_funcall(proc_Hv, RBGSL_ID_call, 3, vX, vV, vHv);
  } else {
    rb_funcall(proc_Hv, RBGSL_ID_call, 4, vX, vp, vV, vHv);
  }
}

static VALUE rb_ool_conmin_function_set_functions(int argc, VALUE *argv, VALUE obj)
{
  ool_conmin_function *F;
  Data_Get_Struct(obj, ool_conmin_function, F);
  set_functions(argc, argv, F);
  return obj;
}

static void set_functions(int argc, VALUE *argv, ool_conmin_function *F)
{
  VALUE ary;
  size_t i;
  if (F->params == NULL) {
    ary = rb_ary_new2(5);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  switch (argc) {
  case 4:
    for (i = 0; i < argc; i++) rb_ary_store(ary, i, argv[i]);
    break;
  default:
    rb_raise(rb_eArgError,"Wrong number of arguments (%d for 4)", argc);
  }
}

static VALUE rb_ool_conmin_function_set_f(VALUE obj, VALUE proc)
{
  ool_conmin_function *F;
  VALUE ary;
  Data_Get_Struct(obj, ool_conmin_function, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(5);
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  rb_ary_store(ary, 0, proc);
  return proc;
}

static VALUE rb_ool_conmin_function_set_df(VALUE obj, VALUE proc)
{
  ool_conmin_function *F;
  VALUE ary;
  Data_Get_Struct(obj, ool_conmin_function, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(5);
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  rb_ary_store(ary, 1, proc);
  return proc;
}

static VALUE rb_ool_conmin_function_set_fdf(VALUE obj, VALUE proc)
{
  ool_conmin_function *F;
  VALUE ary;
  Data_Get_Struct(obj, ool_conmin_function, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(5);
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  rb_ary_store(ary, 2, proc);
  return proc;
}

static VALUE rb_ool_conmin_function_set_Hv(VALUE obj, VALUE proc)
{
  ool_conmin_function *F;
  VALUE ary;
  Data_Get_Struct(obj, ool_conmin_function, F);
  if (F->params == NULL) {
    ary = rb_ary_new2(5);
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  rb_ary_store(ary, 3, proc);
  return proc;
}

static VALUE rb_ool_conmin_function_set_params(VALUE obj, VALUE p)
{
  ool_conmin_function *F;
  Data_Get_Struct(obj, ool_conmin_function, F);
  set_params(F, p);
  return p;
}

static void set_params(ool_conmin_function *F, VALUE p)
{
  VALUE ary;
  if (F->params == NULL) {
    ary = rb_ary_new2(5);
    /*    (VALUE) F->params = ary;*/
    F->params = (void *) ary;
  } else {
    ary = (VALUE) F->params;
  }
  rb_ary_store(ary, 4, p);
}

static VALUE rb_ool_conmin_function_params(VALUE obj)
{
  ool_conmin_function *F;
  Data_Get_Struct(obj, ool_conmin_function, F);
  return rb_ary_entry((VALUE) F->params, 4);;
}

static VALUE rb_ool_conmin_constraint_set(int argc, VALUE *argv, VALUE obj);
static VALUE rb_ool_conmin_constraint_set_n(VALUE obj, VALUE n);
static VALUE rb_ool_conmin_constraint_alloc(int argc, VALUE *argv, VALUE klass)
{
  ool_conmin_constraint *C;
  VALUE obj;
  C = ALLOC(ool_conmin_constraint);
  C->n = 0;
  C->L = NULL;
  C->U = NULL;

  obj = Data_Wrap_Struct(klass, 0, free, C);
  rb_ool_conmin_constraint_set(argc, argv, obj);
  return obj;
}

static VALUE rb_ool_conmin_constraint_set_n(VALUE obj, VALUE n)
{
  ool_conmin_constraint *C;
  if (!FIXNUM_P(n)) rb_raise(rb_eArgError, "Wrong argument type %s (Fixnum expected)",
    rb_class2name(CLASS_OF(n)));
  Data_Get_Struct(obj, ool_conmin_constraint, C);
  C->n = (size_t) FIX2INT(n);
  return n;
}

static VALUE rb_ool_conmin_constraint_set_L(VALUE obj, VALUE vL)
{
  ool_conmin_constraint *C;
  gsl_vector *L;
  CHECK_VECTOR(vL);
  Data_Get_Struct(obj, ool_conmin_constraint, C);
  Data_Get_Struct(vL, gsl_vector, L);
  C->L = L;
  return vL;
}

static VALUE rb_ool_conmin_constraint_set_U(VALUE obj, VALUE vU)
{
  ool_conmin_constraint *C;
  gsl_vector *U;
  CHECK_VECTOR(vU);
  Data_Get_Struct(obj, ool_conmin_constraint, C);
  Data_Get_Struct(vU, gsl_vector, U);
  C->U = U;
  return vU;
}

static VALUE rb_ool_conmin_constraint_set_LU(VALUE obj, VALUE vL, VALUE vU)
{
  rb_ool_conmin_constraint_set_L(obj, vL);
  rb_ool_conmin_constraint_set_U(obj, vU);
  return obj;
}

static VALUE rb_ool_conmin_constraint_set(int argc, VALUE *argv, VALUE obj)
{
  ool_conmin_constraint *C;
  Data_Get_Struct(obj, ool_conmin_constraint, C);
  switch (argc) {
  case 0:
    break;
  case 1:
    rb_ool_conmin_constraint_set_n(obj, argv[0]);
    break;
  case 2:
    rb_ool_conmin_constraint_set_LU(obj, argv[0], argv[1]);
    break;
  case 3:
    rb_ool_conmin_constraint_set_n(obj, argv[0]);
    rb_ool_conmin_constraint_set_LU(obj, argv[1], argv[2]);
    break;
  default:
    rb_raise(rb_eArgError, "Wrong number of arguments (%d for 1-3)", argc);
  }
  return obj;
}

static VALUE rb_ool_conmin_pgrad_parameters_default(VALUE klass)
{
  ool_conmin_pgrad_parameters P;
  VALUE ary;
  ool_conmin_parameters_default(ool_conmin_minimizer_pgrad, (void*) &P);
  ary = create_parameters_ary_pgrad(&P);
  RBGSL_SET_CLASS(ary, cool_conmin_pgrad_parameters);
  return ary;
}

static VALUE rb_ool_conmin_spg_parameters_default(VALUE klass)
{
  ool_conmin_spg_parameters P;
  VALUE ary;
  ool_conmin_parameters_default(ool_conmin_minimizer_spg, (void*) &P);
  ary = create_parameters_ary_spg(&P);
  RBGSL_SET_CLASS(ary, cool_conmin_spg_parameters);
  return ary;
}

static VALUE rb_ool_conmin_gencan_parameters_default(VALUE klass)
{
  ool_conmin_gencan_parameters P;
  VALUE ary;
  ool_conmin_parameters_default(ool_conmin_minimizer_gencan, (void*) &P);
  ary = create_parameters_ary_gencan(&P);
  RBGSL_SET_CLASS(ary, cool_conmin_gencan_parameters);
  return ary;
}

/*************************************************/
void Init_ool(VALUE module)
{
  VALUE mOOL, mConmin;
  VALUE cool_conmin_minimizer;

  mOOL = rb_define_module("OOL");
  mConmin = rb_define_module_under(mOOL, "Conmin");
  cool_conmin_function = rb_define_class_under(mConmin, "Function", cgsl_function);
  cool_conmin_constraint = rb_define_class_under(mConmin, "Constraint", cGSL_Object);
   cool_conmin_minimizer = rb_define_class_under(mConmin, "Minimizer", cGSL_Object);
   cool_conmin_pgrad = rb_define_class_under(cool_conmin_minimizer, "Pgrad", cGSL_Object);
   cool_conmin_spg = rb_define_class_under(cool_conmin_minimizer, "Spg", cGSL_Object);
   cool_conmin_gencan = rb_define_class_under(cool_conmin_minimizer, "Gencan", cGSL_Object);

   def_const(mOOL);

   rb_define_singleton_method(cool_conmin_minimizer, "alloc", rb_ool_conmin_minimizer_alloc, -1);
   rb_define_method(cool_conmin_minimizer, "set", rb_ool_conmin_minimizer_set, -1);
   rb_define_method(cool_conmin_minimizer, "parameters_default", rb_ool_conmin_minimizer_parameters_default, 0);
   rb_define_method(cool_conmin_minimizer, "name", rb_ool_conmin_minimizer_name, 0);
   rb_define_method(cool_conmin_minimizer, "size", rb_ool_conmin_minimizer_size, 0);
   rb_define_method(cool_conmin_minimizer, "f", rb_ool_conmin_minimizer_f, 0);
   rb_define_method(cool_conmin_minimizer, "x", rb_ool_conmin_minimizer_x, 0);
   rb_define_method(cool_conmin_minimizer, "dx", rb_ool_conmin_minimizer_dx, 0);
   rb_define_method(cool_conmin_minimizer, "gradient", rb_ool_conmin_minimizer_gradient, 0);
   rb_define_method(cool_conmin_minimizer, "minimum", rb_ool_conmin_minimizer_minimum, 0);
   rb_define_method(cool_conmin_minimizer, "fcount", rb_ool_conmin_minimizer_fcount, 0);
   rb_define_method(cool_conmin_minimizer, "gcount", rb_ool_conmin_minimizer_gcount, 0);
   rb_define_method(cool_conmin_minimizer, "hcount", rb_ool_conmin_minimizer_hcount, 0);
   rb_define_method(cool_conmin_minimizer, "is_optimal", rb_ool_conmin_minimizer_is_optimal, 0);
   rb_define_method(cool_conmin_minimizer, "is_optimal?", rb_ool_conmin_minimizer_is_optimal2, 0);
   rb_define_method(cool_conmin_minimizer, "iterate", rb_ool_conmin_minimizer_iterate, 0);
   rb_define_method(cool_conmin_minimizer, "restart", rb_ool_conmin_minimizer_restart, 0);
   rb_define_method(cool_conmin_minimizer, "parameters_get", rb_ool_conmin_minimizer_parameters_get, 0);
   rb_define_method(cool_conmin_minimizer, "parameters_set", rb_ool_conmin_minimizer_parameters_set, 1);

   rb_define_singleton_method(cool_conmin_function, "alloc", rb_ool_conmin_function_alloc, -1);
   rb_define_method(cool_conmin_function, "set", rb_ool_conmin_function_set, -1);
   rb_define_method(cool_conmin_function, "set_n", rb_ool_conmin_function_set_n, 1);
  rb_define_alias(cool_conmin_function, "n=", "set_n");
   rb_define_method(cool_conmin_function, "n", rb_ool_conmin_function_n, 0);
   rb_define_method(cool_conmin_function, "params", rb_ool_conmin_function_params, 0);
   rb_define_method(cool_conmin_function, "set_params", rb_ool_conmin_function_set_params, 1);
  rb_define_alias(cool_conmin_function, "params=", "set_params");
   rb_define_method(cool_conmin_function, "set_functions", rb_ool_conmin_function_set_functions, 1);
  rb_define_alias(cool_conmin_function, "functions=", "set_functions");
   rb_define_method(cool_conmin_function, "set_f", rb_ool_conmin_function_set_f, 1);
  rb_define_alias(cool_conmin_function, "f=", "set_f");
   rb_define_method(cool_conmin_function, "set_df", rb_ool_conmin_function_set_df, 1);
  rb_define_alias(cool_conmin_function, "df=", "set_df");
   rb_define_method(cool_conmin_function, "set_fdf", rb_ool_conmin_function_set_fdf, 1);
  rb_define_alias(cool_conmin_function, "fdf=", "set_fdf");
   rb_define_method(cool_conmin_function, "set_Hv", rb_ool_conmin_function_set_Hv, 1);
  rb_define_alias(cool_conmin_function, "Hv=", "set_Hv");

  rb_define_singleton_method(cool_conmin_constraint, "alloc", rb_ool_conmin_constraint_alloc,
    -1);
  rb_define_method(cool_conmin_constraint, "set", rb_ool_conmin_constraint_set, -1);
  rb_define_method(cool_conmin_constraint, "set_n", rb_ool_conmin_constraint_set_n, 1);
  rb_define_alias(cool_conmin_constraint, "n=", "set_n");
  rb_define_method(cool_conmin_constraint, "set_L", rb_ool_conmin_constraint_set_L, 1);
  rb_define_alias(cool_conmin_constraint, "L=", "set_L");
  rb_define_method(cool_conmin_constraint, "set_U", rb_ool_conmin_constraint_set_U, 1);
  rb_define_alias(cool_conmin_constraint, "U=", "set_U");
  rb_define_method(cool_conmin_constraint, "set_LU", rb_ool_conmin_constraint_set_LU, 2);
  rb_define_alias(cool_conmin_constraint, "LU=", "set_LU");

  cool_conmin_pgrad_parameters = rb_define_class_under(cool_conmin_pgrad, "Parameters",
      rb_cArray);
  cool_conmin_spg_parameters = rb_define_class_under(cool_conmin_spg, "Parameters",
      rb_cArray);
  cool_conmin_gencan_parameters = rb_define_class_under(cool_conmin_gencan, "Parameters",
      rb_cArray);
  rb_define_singleton_method(cool_conmin_pgrad, "parameters_default",
    rb_ool_conmin_pgrad_parameters_default, 0);
  rb_define_singleton_method(cool_conmin_spg, "parameters_default",
    rb_ool_conmin_spg_parameters_default, 0);
  rb_define_singleton_method(cool_conmin_gencan, "parameters_default",
    rb_ool_conmin_gencan_parameters_default, 0);
}

#endif
