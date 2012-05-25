// Not exactly sure this function will work. I translated it from Ruby::Rational modulo.
inline TYPE mod4(INT xn, const INT xd, INT yn, const INT yd) {
  LONG_TYPE prod;
  LONG_INT floor_div;

  if (yn < 0) {
    xn = -xn;
    yn = -yn;
  }

  floor_div = (xn * yn) / (xd * yd);
  prod = mul4(yn, yd, floor_div, 1);

  return sub4(xn, yn, prod.n, prod.d);
}
