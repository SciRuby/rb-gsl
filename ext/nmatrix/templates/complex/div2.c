inline TYPE div2(const TYPE x, const TYPE y) {
  TYPE denom = y.r * y.r + x.i * x.i;
  return (struct TYPE) {
    (x.r * y.r + x.i * y.i) / denom,
    (x.r * y.i - x.i * y.r) / denom
  };
}
