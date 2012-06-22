inline TYPE add4(INT xn, INT xd, INT yn, INT yd) {
  TYPE result;

  INT ig   = gcf(xd, yd);
  INT x    = xn * (yd / ig);
  INT y    = yn * (xd / ig);
  INT z    = x+y;

  y        = xd / ig;
  ig       = gcf(xd, ig);
  result.n = z / ig;
  x        = yd / ig;
  result.d = x*y;

  return result;
}
