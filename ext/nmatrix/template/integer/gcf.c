inline TYPE gcf(TYPE x, TYPE y) {
  TYPE t;

  if (x < 0) x = -x;
  if (y < 0) y = -y;

  if (x == 0) return y;
  if (y == 0) return x;

  while (x > 0) {
    t = x;
    x = y % x;
    y = t;
  }

  return y;
}