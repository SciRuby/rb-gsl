#include "rb_gsl_dirac.h"

static VALUE cgsl_matrix_complex_const;
static VALUE cPauli;
static VALUE cAlpha;
static VALUE cGamma;
static VALUE cLambda;

static gsl_matrix_complex *Pauli[3];
static gsl_matrix_complex *Alpha[3];
static gsl_matrix_complex *Beta;
static gsl_matrix_complex *Gamma[5];
static gsl_matrix_complex *Eye2, *Eye4;
static gsl_matrix_complex *IEye2, *IEye4;
static gsl_matrix_complex *Lambda[8];

static VALUE VPauli[3];
static VALUE VAlpha[3];
static VALUE VGamma[5];
static VALUE VEye2, VEye4, VIEye2, VIEye4;
static VALUE VLambda[8];

static void Init_gsl_dirac_const(VALUE module);
static void Init_gsl_dirac_common(VALUE module);
static void Init_gsl_dirac_test(VALUE module);

static VALUE rb_dirac_refuse_set(int argc, VALUE *argv, VALUE obj)
{
  rb_raise(rb_eRuntimeError, "Cannot modify this object.");
}

static VALUE rb_dirac_commute(VALUE obj, VALUE mm1, VALUE mm2)
{
  gsl_matrix_complex *m1, *m2;
  gsl_matrix_complex *mnew1, *mnew2;
  CHECK_MATRIX_COMPLEX(mm1);
  CHECK_MATRIX_COMPLEX(mm2);
  Data_Get_Struct(mm1, gsl_matrix_complex, m1);
  Data_Get_Struct(mm2, gsl_matrix_complex, m2);
  mnew1 = gsl_matrix_complex_alloc(m1->size1, m1->size2);
  mnew2 = gsl_matrix_complex_alloc(m1->size1, m1->size2);
  gsl_matrix_complex_mul(mnew1, m1, m2);
  gsl_matrix_complex_mul(mnew2, m2, m1);
  gsl_matrix_complex_sub(mnew1, mnew2);
  gsl_matrix_complex_free(mnew2);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, 
			  mnew1);
}

static VALUE rb_dirac_anticommute(VALUE obj, VALUE mm1, VALUE mm2)
{
  gsl_matrix_complex *m1, *m2;
  gsl_matrix_complex *mnew1, *mnew2;
  CHECK_MATRIX_COMPLEX(mm1);
  CHECK_MATRIX_COMPLEX(mm2);
  Data_Get_Struct(mm1, gsl_matrix_complex, m1);
  Data_Get_Struct(mm2, gsl_matrix_complex, m2);
  mnew1 = gsl_matrix_complex_alloc(m1->size1, m1->size2);
  mnew2 = gsl_matrix_complex_alloc(m1->size1, m1->size2);
  gsl_matrix_complex_mul(mnew1, m1, m2);
  gsl_matrix_complex_mul(mnew2, m2, m1);
  gsl_matrix_complex_add(mnew1, mnew2);
  gsl_matrix_complex_free(mnew2);
  return Data_Wrap_Struct(cgsl_matrix_complex, 0, gsl_matrix_complex_free, 
			  mnew1);
}

static void Init_gsl_dirac_common(VALUE module)
{
  /*  VALUE cBeta;*/
  rb_define_singleton_method(module, "commute", rb_dirac_commute, 2);
  rb_define_singleton_method(module, "anticommute", rb_dirac_anticommute, 2);

  cgsl_matrix_complex_const = rb_define_class_under(module, "Const",
						    cgsl_matrix_complex);
  rb_define_method(cgsl_matrix_complex_const, "set", rb_dirac_refuse_set, -1);

  cPauli = rb_define_class_under(module, "Pauli", cgsl_matrix_complex_const);
  cAlpha = rb_define_class_under(module, "Alpha", cgsl_matrix_complex_const);
  /*  cBeta = rb_define_class_under(module, "BetaMatrix", cgsl_matrix_complex_const);*/
  cGamma = rb_define_class_under(module, "Gamma", cgsl_matrix_complex_const);
  cLambda = rb_define_class_under(module, "Lambda", cgsl_matrix_complex_const);
}


static void define_eye(VALUE module)
{
  gsl_complex z;

  Eye2 = gsl_matrix_complex_calloc(2, 2);
  VEye2 = Data_Wrap_Struct(cgsl_matrix_complex_const, 0, 
			   gsl_matrix_complex_free, Eye2);
  z.dat[0] = 1; z.dat[1] = 0;
  gsl_matrix_complex_set(Eye2, 0, 0, z);
  gsl_matrix_complex_set(Eye2, 1, 1, z);
  rb_define_const(module, "Eye2", VEye2);

  Eye4 = gsl_matrix_complex_calloc(4, 4);
  VEye4 = Data_Wrap_Struct(cgsl_matrix_complex_const, 0, 
			   gsl_matrix_complex_free, Eye4);
  z.dat[0] = 1; z.dat[1] = 0;
  gsl_matrix_complex_set(Eye4, 0, 0, z);
  gsl_matrix_complex_set(Eye4, 1, 1, z);
  gsl_matrix_complex_set(Eye4, 2, 2, z);
  gsl_matrix_complex_set(Eye4, 3, 3, z);
  rb_define_const(module, "Eye4", VEye4);

  IEye2 = gsl_matrix_complex_calloc(2, 2);
  VIEye2 = Data_Wrap_Struct(cgsl_matrix_complex_const, 0, 
			    gsl_matrix_complex_free, IEye2);
  z.dat[0] = 0; z.dat[1] = 1;
  gsl_matrix_complex_set(IEye2, 0, 0, z);
  gsl_matrix_complex_set(IEye2, 1, 1, z);
  rb_define_const(module, "IEye2", VIEye2);

  IEye4 = gsl_matrix_complex_calloc(4, 4);
  VIEye4 = Data_Wrap_Struct(cgsl_matrix_complex_const, 0, 
			   gsl_matrix_complex_free, IEye4);
  gsl_matrix_complex_set(IEye4, 0, 0, z);
  gsl_matrix_complex_set(IEye4, 1, 1, z);
  gsl_matrix_complex_set(IEye4, 2, 2, z);
  gsl_matrix_complex_set(IEye4, 3, 3, z);
  rb_define_const(module, "IEye4", VIEye4);
}

static void define_pauli(VALUE module)
{
  gsl_complex z;

  Pauli[0] = gsl_matrix_complex_calloc(2, 2);
  VPauli[0] = Data_Wrap_Struct(cPauli, 0, 
			   gsl_matrix_complex_free, Pauli[0]);
  z.dat[0] = 1; z.dat[1] = 0;
  gsl_matrix_complex_set(Pauli[0], 0, 1, z);
  gsl_matrix_complex_set(Pauli[0], 1, 0, z);
  rb_define_const(module, "Pauli1", VPauli[0]);

  Pauli[1] = gsl_matrix_complex_calloc(2, 2);
  VPauli[1] = Data_Wrap_Struct(cPauli, 0, 
			   gsl_matrix_complex_free, Pauli[1]);
  z.dat[0] = 0; z.dat[1] = -1;
  gsl_matrix_complex_set(Pauli[1], 0, 1, z);
  z.dat[0] = 0; z.dat[1] = 1;
  gsl_matrix_complex_set(Pauli[1], 1, 0, z);
  rb_define_const(module, "Pauli2", VPauli[1]);

  Pauli[2] = gsl_matrix_complex_calloc(2, 2);
  VPauli[2] = Data_Wrap_Struct(cPauli, 0, 
			   gsl_matrix_complex_free, Pauli[2]);
  z.dat[0] = 1; z.dat[1] = 0;
  gsl_matrix_complex_set(Pauli[2], 0, 0, z);
  z.dat[0] = -1; z.dat[1] = 0;
  gsl_matrix_complex_set(Pauli[2], 1, 1, z);
  rb_define_const(module, "Pauli3", VPauli[2]);
}

static void define_beta(VALUE module)
{
  gsl_complex z;

  Beta = gsl_matrix_complex_calloc(4, 4);
  VGamma[0] = Data_Wrap_Struct(cGamma, 0, 
			   gsl_matrix_complex_free, Beta);
  z.dat[0] = 1; z.dat[1] = 0;
  gsl_matrix_complex_set(Beta, 0, 0, z);
  gsl_matrix_complex_set(Beta, 1, 1, z);
  z.dat[0] = -1; z.dat[1] = 0;
  gsl_matrix_complex_set(Beta, 2, 2, z);
  gsl_matrix_complex_set(Beta, 3, 3, z);
  rb_define_const(module, "Beta", VGamma[0]);
  rb_define_const(module, "Gamma0", VGamma[0]);
}

static void define_alpha(VALUE module)
{
  size_t i, j, k;
  char name[7];
  for (i = 0; i < 3; i++) {
    Alpha[i] = gsl_matrix_complex_calloc(4, 4);

    for (j = 2; j < 4; j++) {
      for (k = 0; k < 2; k++) {
	gsl_matrix_complex_set(Alpha[i], j, k, 
			       gsl_matrix_complex_get(Pauli[i], j-2, k));
      }
    }
    for (j = 0; j < 2; j++) {
      for (k = 2; k < 4; k++) {
	gsl_matrix_complex_set(Alpha[i], j, k, 
			       gsl_matrix_complex_get(Pauli[i], j, k-2));
      }
    }
    VAlpha[i] = Data_Wrap_Struct(cAlpha, 0, 
			    gsl_matrix_complex_free, Alpha[i]);
    sprintf(name, "Alpha%d", (int) i+1);
    rb_define_const(module, name, VAlpha[i]);
  }
  
}

static void define_gamma(VALUE module)
{
  size_t i;
  char name[7];
  gsl_complex z;
  for (i = 1; i <= 3; i++) {
    Gamma[i] = gsl_matrix_complex_calloc(4, 4);
    gsl_matrix_complex_mul(Gamma[i], Beta, Alpha[i-1]);
    VGamma[i] = Data_Wrap_Struct(cGamma, 0, 
			    gsl_matrix_complex_free, Gamma[i]);
    sprintf(name, "Gamma%d", (int) i);
    rb_define_const(module, name, VGamma[i]);
  }
  Gamma[4] = gsl_matrix_complex_calloc(4, 4);
  z.dat[0] = 1.0; z.dat[1] = 0.0;
  gsl_matrix_complex_set(Gamma[4], 0, 2, z);
  gsl_matrix_complex_set(Gamma[4], 1, 3, z);
  gsl_matrix_complex_set(Gamma[4], 2, 0, z);
  gsl_matrix_complex_set(Gamma[4], 3, 1, z);
  VGamma[4] = Data_Wrap_Struct(cGamma, 0, 
			  gsl_matrix_complex_free, Gamma[4]);
  rb_define_const(module, "Gamma5", VGamma[4]);
}

static void define_lambda(VALUE module)
{
  gsl_complex z1, zm1, zi, zmi;
  size_t i;
  char name[8];
  double sqrt3 = sqrt(3.0);
  z1.dat[0] = 1; z1.dat[1] = 0;
  zm1.dat[0] = -1; zm1.dat[1] = 0;
  zi.dat[0] = 0; zi.dat[1] = 1;
  zmi.dat[0] = 0; zmi.dat[1] = -1;
  for (i = 0; i < 8; i++) {
    Lambda[i] = gsl_matrix_complex_calloc(3, 3);
    VLambda[i] = Data_Wrap_Struct(cLambda, 0, 
			    gsl_matrix_complex_free, Lambda[i]);
    sprintf(name, "Lambda%d", (int) i+1);
    rb_define_const(module, name, VLambda[i]);
  }
  gsl_matrix_complex_set(Lambda[0], 0, 1, z1);
  gsl_matrix_complex_set(Lambda[0], 1, 0, z1);
  gsl_matrix_complex_set(Lambda[1], 0, 1, zmi);
  gsl_matrix_complex_set(Lambda[1], 1, 0, zi);
  gsl_matrix_complex_set(Lambda[2], 0, 0, z1);
  gsl_matrix_complex_set(Lambda[2], 1, 1, zm1);
  gsl_matrix_complex_set(Lambda[3], 0, 2, z1);
  gsl_matrix_complex_set(Lambda[3], 2, 0, z1);
  gsl_matrix_complex_set(Lambda[4], 0, 2, zmi);
  gsl_matrix_complex_set(Lambda[4], 2, 0, zi);
  gsl_matrix_complex_set(Lambda[5], 1, 2, z1);
  gsl_matrix_complex_set(Lambda[5], 2, 1, z1);
  gsl_matrix_complex_set(Lambda[6], 1, 2, zmi);
  gsl_matrix_complex_set(Lambda[6], 2, 1, zi);
  gsl_matrix_complex_set(Lambda[7], 0, 0, gsl_complex_mul_real(z1, 1.0/sqrt3));
  gsl_matrix_complex_set(Lambda[7], 1, 1, gsl_complex_mul_real(z1, 1.0/sqrt3));
  gsl_matrix_complex_set(Lambda[7], 2, 2, gsl_complex_mul_real(z1, -2.0/sqrt3));
}

static void Init_gsl_dirac_const(VALUE module)
{
  define_eye(module);
  define_pauli(module);
  define_beta(module);
  define_alpha(module);
  define_gamma(module);
  define_lambda(module);
}

static int matrix_is_equal(gsl_matrix_complex *m1, gsl_matrix_complex *m2, gsl_complex *c)
{
  gsl_complex a, b, ab, absave;
  double eps = 1e-6;
  size_t i, j;
  absave.dat[0] = 99999;
  absave.dat[1] = 99999;
  if (m1->size1 != m2->size1 ||  m1->size2 != m2->size2) return 0;
  for (i = 0; i < m1->size1; i++) {
    for (j = 0; j < m1->size2; j++) {
      a = gsl_matrix_complex_get(m1, i, j);
      b = gsl_matrix_complex_get(m2, i, j);
      if (!gsl_fcmp(gsl_complex_abs(b), 0.0, eps)) continue;
      ab = gsl_complex_div(a, b);
      if (!gsl_fcmp(gsl_complex_abs(ab), 0.0, eps)) continue;
      if ((int) absave.dat[0] == 99999) absave = ab;
      if (gsl_fcmp(ab.dat[0], absave.dat[0], eps)) return 0;
      if (gsl_fcmp(ab.dat[1], absave.dat[1], eps)) return 0;
    }
  }
  if ((int) absave.dat[0] == 99999) return 0;
  *c = ab;
  return 1;
}

static VALUE rb_Dirac_matrix_is_equal(int argc, VALUE *argv, VALUE obj)
{
  gsl_complex ztmp, *z;
  gsl_matrix_complex *m1, *m2;
  VALUE vz;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    CHECK_MATRIX_COMPLEX(argv[0]);
    CHECK_MATRIX_COMPLEX(argv[1]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, m1);
    Data_Get_Struct(argv[1], gsl_matrix_complex, m2);
    if (matrix_is_equal(m1, m2, &ztmp)) {
      vz = Data_Make_Struct(cgsl_complex, gsl_complex, 0, free, z);
      *z = ztmp;
      return vz;
    } else {
      return Qfalse;
    }
    break;
  default:
    CHECK_MATRIX_COMPLEX(argv[0]);
    Data_Get_Struct(obj, gsl_matrix_complex, m1);
    Data_Get_Struct(argv[0], gsl_matrix_complex, m2);
    if (matrix_is_equal(m1, m2, &ztmp)) {
      vz = Data_Make_Struct(cgsl_complex, gsl_complex, 0, free, z);
      *z = ztmp;
      return vz;
    } else {
      return Qfalse;
    }
    break;
  }
}

#define NUM 20
static VALUE rb_Dirac_matrix_whoami(int argc, VALUE *argv, VALUE obj)
{
  VALUE array[NUM] = {VPauli[0], VPauli[1], VPauli[2], 
		      VGamma[0], VGamma[1], VGamma[2], VGamma[3],
		      VGamma[4], VEye2, VEye4, VIEye2, VIEye4,
		      VLambda[0], VLambda[1], VLambda[2], VLambda[3],
		      VLambda[4], VLambda[5], VLambda[6], VLambda[7]};

  const char *name[NUM] = {"Pauli1", "Pauli2", "Pauli3", 
		     "Gamma0", "Gamma1", "Gamma2", "Gamma3", "Gamma5",
		     "Eye2", "Eye4", "IEye2", "IEye4", "Lambda1", "Lambda2",
		     "Lambda3", "Lambda4", "Lambda5", "Lambda6",
		     "Lambda7", "Lambda8"};
  gsl_matrix_complex *m1, *m2;
  VALUE vz;
  gsl_complex ztmp, *z;
  size_t i;
  switch (TYPE(obj)) {
  case T_MODULE:
  case T_CLASS:
  case T_OBJECT:
    if (argc != 1) rb_raise(rb_eArgError, "matrix not given");
    CHECK_MATRIX_COMPLEX(argv[0]);
    Data_Get_Struct(argv[0], gsl_matrix_complex, m1);
    break;
  default:
    Data_Get_Struct(obj, gsl_matrix_complex, m1);
    break;
  }
  for (i = 0; i < NUM; i++) {
    Data_Get_Struct(array[i], gsl_matrix_complex, m2);
    if(matrix_is_equal(m1, m2, &ztmp)) {
      vz = Data_Make_Struct(cgsl_complex, gsl_complex, 0, free, z);
      *z = ztmp;
      return rb_ary_new3(3, array[i], rb_str_new2(name[i]), vz);
    }
  }
  return Qfalse;
}
#undef NUM

static void Init_gsl_dirac_test(VALUE module)
{
  rb_define_singleton_method(module, "is_equal?", rb_Dirac_matrix_is_equal, -1);
  rb_define_method(cgsl_matrix_complex, "is_equal?", rb_Dirac_matrix_is_equal, -1);
  rb_define_singleton_method(module, "whatisthis", rb_Dirac_matrix_whoami, -1);
  rb_define_method(cgsl_matrix_complex, "whoami", rb_Dirac_matrix_whoami, -1);
}

void Init_gsl_dirac(VALUE module)
{
  VALUE mDirac;
  mDirac = rb_define_module_under(module, "Dirac");
  Init_gsl_dirac_common(mDirac);
  Init_gsl_dirac_const(mDirac);
  Init_gsl_dirac_test(mDirac);
}
