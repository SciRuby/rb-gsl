// Given two vectors of some type, determine if contents are numerically equal.
//
// This function is specific for complex, float, and double. Other types are numeqeq.
inline bool %%TYPE_ABBREV%%eqeq(const %%TYPE%%* x, const %%TYPE%%* y, const size_t len, const size_t unused) {
  size_t p;
  for (p = 0; p < len; ++p) {
    if (%%TYPE x[p] != y[p]%%) return false;
  }
  return true;
}
