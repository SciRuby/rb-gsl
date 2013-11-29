
#include "rb_gsl.h"
#ifdef GSL_1_9_LATER
#include "gsl/gsl_bspline.h"

static VALUE cBSWS;

static VALUE rb_gsl_bspline_alloc(VALUE klass, VALUE k, VALUE n)
{
	gsl_bspline_workspace *w;
	w = gsl_bspline_alloc(FIX2INT(k), FIX2INT(n));
	return Data_Wrap_Struct(klass, 0, gsl_bspline_free, w);
}

static VALUE rb_gsl_bspline_ncoeffs(VALUE obj)
{
	gsl_bspline_workspace *w;
	Data_Get_Struct(obj, gsl_bspline_workspace, w);
	return INT2FIX((int)gsl_bspline_ncoeffs(w));
}
static VALUE rb_gsl_bspline_order(VALUE obj)
{
	gsl_bspline_workspace *w;
	Data_Get_Struct(obj, gsl_bspline_workspace, w);
	return INT2FIX((int)gsl_bspline_order(w));
}
static VALUE rb_gsl_bspline_nbreak(VALUE obj)
{
	gsl_bspline_workspace *w;
	Data_Get_Struct(obj, gsl_bspline_workspace, w);
	return INT2FIX((int)gsl_bspline_nbreak(w));
}
static VALUE rb_gsl_bspline_breakpoint(VALUE obj, VALUE i)
{
	gsl_bspline_workspace *w;
	Data_Get_Struct(obj, gsl_bspline_workspace, w);
	return rb_float_new(gsl_bspline_breakpoint(FIX2INT(i), w));
}
static VALUE rb_gsl_bspline_knots(VALUE obj, VALUE b)
{
	gsl_bspline_workspace *w;
	gsl_vector *bpts;
	CHECK_VECTOR(b);
	Data_Get_Struct(obj, gsl_bspline_workspace, w);
	Data_Get_Struct(b, gsl_vector, bpts);
	gsl_bspline_knots(bpts, w);
	return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, w->knots);
}
static VALUE rb_gsl_bspline_knots_uniform(int argc, VALUE *argv, VALUE obj)
{
	gsl_bspline_workspace *w;
	int argc2;
	switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (!rb_obj_is_kind_of(argv[argc-1], cBSWS)) {
      rb_raise(rb_eTypeError, "Wrong argument type %s (GSL::BSpline expected)",
        rb_class2name(CLASS_OF(argv[argc-1])));
    }    
    Data_Get_Struct(argv[argc-1], gsl_bspline_workspace, w);
    argc2 = argc-1;
    break;
  default:
    Data_Get_Struct(obj, gsl_bspline_workspace, w);  	
    argc2 = argc;
	}
  if (argc2 != 2) rb_raise(rb_eArgError, "Wrong number of arguments.");
	gsl_bspline_knots_uniform(NUM2DBL(argv[0]), NUM2DBL(argv[1]), w);
	return Data_Wrap_Struct(cgsl_vector_view_ro, 0, NULL, w->knots);	
}
static VALUE rb_gsl_bspline_eval(int argc, VALUE *argv, VALUE obj)
{
	gsl_bspline_workspace *w;
	double x;
	gsl_vector *B;
	VALUE vB;

	Data_Get_Struct(obj, gsl_bspline_workspace, w);

	switch (argc) {
	case 2:
		CHECK_VECTOR(argv[1]);
		Data_Get_Struct(argv[1], gsl_vector, B);
		vB = argv[1];
		x = NUM2DBL(argv[0]);	
		break;
	case 1:
		x = NUM2DBL(argv[0]);
		B = gsl_vector_alloc(w->nbreak+w->k-2);
		vB = Data_Wrap_Struct(cgsl_vector, 0, gsl_vector_free, B);
		break;
	default:
		rb_raise(rb_eArgError, "Wrong number of arguments (%d for 1 or 2)", argc);
	}

	gsl_bspline_eval(x, B, w);

	return vB;
}
#ifdef GSL_1_13_LATER
static VALUE rb_gsl_bspline_greville_abscissa(VALUE obj, VALUE i)
{
  gsl_bspline_workspace *w;
  Data_Get_Struct(obj, gsl_bspline_workspace, w);
  return rb_float_new(gsl_bspline_greville_abscissa(i, w));
}
#endif

void Init_bspline(VALUE module)
{
  cBSWS = rb_define_class_under(module, "BSpline", cGSL_Object);
  
  rb_define_singleton_method(cBSWS, "alloc", rb_gsl_bspline_alloc, 2);
  
  rb_define_method(cBSWS, "ncoeffs", rb_gsl_bspline_ncoeffs, 0);
  rb_define_method(cBSWS, "order", rb_gsl_bspline_order, 0);
  rb_define_method(cBSWS, "nbreak", rb_gsl_bspline_nbreak, 0);
  rb_define_method(cBSWS, "breakpoint", rb_gsl_bspline_breakpoint, 1);
  rb_define_method(cBSWS, "knots", rb_gsl_bspline_knots, 1);
  rb_define_method(cBSWS, "knots_uniform", rb_gsl_bspline_knots_uniform, -1);
  rb_define_singleton_method(cBSWS, "knots_uniform", rb_gsl_bspline_knots_uniform, -1);	
  rb_define_method(cBSWS, "eval", rb_gsl_bspline_eval, -1);	

#ifdef GSL_1_13_LATER
  rb_define_method(cBSWS, "greville_abscissa", rb_gsl_bspline_greville_abscissa, 1);
#endif

}
#endif
