inline TYPE downcast(const LONG_TYPE n) {
  return (struct TYPE) { n.r, n.i };
}
