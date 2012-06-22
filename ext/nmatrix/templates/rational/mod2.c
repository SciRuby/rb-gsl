inline TYPE mod2(const TYPE x, const TYPE y) {
  return mod4(x.n, x.d, y.n, y.d);
}
