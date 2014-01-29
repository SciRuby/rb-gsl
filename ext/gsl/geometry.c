#include "rb_gsl.h"

/*!
  Counter-clockwise rotation around the X-axis 

  / x' \    /  1     0      0    \ / x \
  | y' |  = |  0  cos_th -sin_th | | y |
  \ z' /    \  0  sin_th  cos_th / \ z /

 */
void vector3_rotateX(const double x[3] /*!< Input */, 
                        double theta /*!< Rotation angle */, 
                        double xout[3]  /*!< Output */)
{
  double a, b, c;
  double costheta, sintheta;
  
  costheta = cos(theta);
  sintheta = sin(theta);

  a = x[0];
  b = x[1]*costheta - x[2]*sintheta;
  c = x[1]*sintheta + x[2]*costheta;

  xout[0] = a;
  xout[1] = b;
  xout[2] = c;
}

/*!
  Counter-clockwise rotation around the Y-axis 
  (Note the sign of the matrix element.)

  / x' \    /  cos_th     0   sin_th \ / x \
  | y' |  = |    0        1      0   | | y |
  \ z' /    \ -sin_th     0   cos_th / \ z /

 */
void vector3_rotateY(const double x[3] /*!< Input */, 
                        double theta /*!< Rotation angle */,
                        double xout[3] /*!< Output */)
{
  double a, b, c;
  double costheta, sintheta;

  costheta = cos(theta);
  sintheta = sin(theta);

  a =  x[0]*costheta + x[2]*sintheta;
  b =  x[1];
  c = -x[0]*sintheta + x[2]*costheta;

  xout[0] = a;
  xout[1] = b;
  xout[2] = c;
}

/*!
  Counter-clockwise rotation around the Z-axis 

  / x' \    /  cos_th  -sin_th  0 \ / x \
  | y' |  = |  sin_th   cos_th  0 | | y |
  \ z' /    \    0        0     1 / \ z /

 */
void vector3_rotateZ(const double x[3] /*!< Input */, 
                        double theta  /*!< Rotation angle */, 
                        double xout[3] /*!< Output */)
{
  double a, b, c;
  double costheta, sintheta;

  costheta = cos(theta);
  sintheta = sin(theta);

  a = x[0]*costheta - x[1]*sintheta;
  b = x[0]*sintheta + x[1]*costheta;
  c = x[2];

  xout[0] = a;
  xout[1] = b;
  xout[2] = c;
}

/*!
  Rotate a 3-vector.

  If flag != 1, cos/sin values stored in the static variables are used 
  for efficiency.
  The input and the output vectors (x, xout) can be the same pointer.

  / x' \     /   cos_phi -sin_phi   0  \ /  cos_th   0   sin_th \  / x \
  | y' |  =  |   sin_phi  cos_phi   0  | |     0     1          |  | y |
  \ z' /     \      0        0      1  / \ -sin_th   0   cos_th /  \ z /

 */
void vector3_rotate(const double x[3], double theta, double phi,
 double xout[3])
{
  vector3_rotateY(x, theta,xout);
  vector3_rotateZ(xout, phi, xout);
}

static VALUE rb_gsl_vector_rotateX(VALUE obj, VALUE angle)
{
  gsl_vector *v;
  Data_Get_Struct(obj, gsl_vector, v);
  vector3_rotateX(v->data, NUM2DBL(angle), v->data);
  return obj;
}
static VALUE rb_gsl_vector_rotateY(VALUE obj, VALUE angle)
{
  gsl_vector *v;
  Data_Get_Struct(obj, gsl_vector, v);
  vector3_rotateY(v->data, NUM2DBL(angle), v->data);
  return obj;
}
static VALUE rb_gsl_vector_rotateZ(VALUE obj, VALUE angle)
{
  gsl_vector *v;
  Data_Get_Struct(obj, gsl_vector, v);
  vector3_rotateZ(v->data, NUM2DBL(angle), v->data);
  return obj;
}
static VALUE rb_gsl_vector_rotate(VALUE obj, VALUE theta, VALUE phi)
{
  gsl_vector *v;
  Data_Get_Struct(obj, gsl_vector, v);
  vector3_rotate(v->data, NUM2DBL(theta), NUM2DBL(phi), v->data);
  return obj;
}

void Init_geometry(VALUE module)
{
  rb_define_method(cgsl_vector, "rotateX", rb_gsl_vector_rotateX, 1);
  rb_define_method(cgsl_vector, "rotateY", rb_gsl_vector_rotateY, 1);
  rb_define_method(cgsl_vector, "rotateZ", rb_gsl_vector_rotateZ, 1);
  rb_define_method(cgsl_vector, "rotate", rb_gsl_vector_rotate, 2);
}
