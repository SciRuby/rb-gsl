// Given two vectors of some type, determine if contents are numerically equal.
//
// This function is specific for complex, float, and double. Other types are numeqeq.
inline bool eqeq(const TYPE* x, const TYPE* y, const int len, const int unused) {
  int p;
  for (p = 0; p < len; ++p) {
    if (x[p] != y[p]) return false;
  }
  return true;
}
