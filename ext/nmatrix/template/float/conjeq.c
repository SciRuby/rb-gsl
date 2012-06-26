// Same function pointer type as eqeq (in shared templates)
inline bool conjeq(const TYPE* x, const TYPE* y, const int len, const int unused) {
  int p;
  for (p = 0; p < len; ++p) {
    if (x[p] != -y[p]) return false;
  }
  return true;
}