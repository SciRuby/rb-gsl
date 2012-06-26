inline TYPE add2(const TYPE x, const TYPE y) {
  return (struct TYPE) { x.r + y.r, x.i + y.i };
}
