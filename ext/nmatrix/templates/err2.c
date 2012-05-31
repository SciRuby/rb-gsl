TYPE err2(TYPE x, TYPE y) {
  rb_raise(nm_eDataTypeError, "illegal math operation with this data type");
  return x;
}