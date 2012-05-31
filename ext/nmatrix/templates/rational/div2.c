inline TYPE div2(const TYPE x, const TYPE y) {
  if (y.n < 0)
    return mul4(-x.n, x.d, -y.n, y.d);
  else
    return mul4(x.n, x.d, y.n, y.d);
}
