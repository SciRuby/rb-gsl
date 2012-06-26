inline TYPE mul2(const TYPE x, const TYPE y) {
  return (struct TYPE) { x.r * y.r - x.i * y.i, x.r * y.i - x.i * y.r };
}
