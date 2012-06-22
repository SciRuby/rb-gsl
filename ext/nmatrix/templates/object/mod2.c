inline TYPE mod2(const TYPE x, const TYPE y) {
  return rb_funcall(x, rb_intern("%"), 1, y);
}