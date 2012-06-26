inline TYPE downcast(const LONG_TYPE n) {
  return (struct TYPE) { n.n, n.d };
}
