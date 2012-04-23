
inline %%TYPE%% BOOL2%%= dtype.id.to_s.upcase%%(bool expr) {
  %%TYPE%% result;
  result.r = expr;
  result.i = 0;
  return result;
}



inline %%TYPE%% %%TYPE_ABBREV%%_mul(double ar, double ai, double br, double bi) {
  %%TYPE%% result;

  result.r = ar * br - ai * bi;
  result.i = ar * bi + br * ai;

  return result;
}


inline %%TYPE%% %%TYPE_ABBREV%%_div(double ar, double ai, double br, double bi) {
  %%TYPE%% result;
  double denom = br * br + bi * bi;

  result.r = (ar * br + ai * bi) / denom;
  result.i = (ai * br - ar * bi) / denom;

  return result;
}


inline %%TYPE%% %%TYPE_ABBREV%%_add(%%= dtype.sym == :complex128 ? "double ar, double ai, double br, double bi" : "float ar, float ai, float br, float bi"%%) {
  %%TYPE%% result;

  result.r = ar + br;
  result.i = ai + bi;

  return result;
}


inline %%TYPE%% %%TYPE_ABBREV%%_sub(%%= dtype.sym == :complex128 ? "double ar, double ai, double br, double bi" : "float ar, float ai, float br, float bi"%%) {
  %%TYPE%% result;

  result.r = ar - br;
  result.i = ai - bi;

  return result;
}


inline %%TYPE%% %%TYPE_ABBREV%%_mod(%%= dtype.sym == :complex128 ? "double ar, double ai, double br, double bi" : "float ar, float ai, float br, float bi"%%) {
  %%TYPE%% result;
  rb_raise(rb_eNotImpError, "modulo arithmetic for complex numbers not yet implemented");

  return result;
}


inline %%TYPE%% %%TYPE_ABBREV%%_bang(%%= dtype.sym == :complex128 ? "double ar, double ai" : "float ar, float ai"%%) {
  %%TYPE%% result = {!ar, 0};
  return result;
}

inline %%TYPE%% %%TYPE_ABBREV%%_negate(%%= dtype.sym == :complex128 ? "double ar, double ai" : "float ar, float ai"%%) {
  %%TYPE%% result;
  result.r = -ar;
  result.i = -ai;
  return result;
}

// Same function pointer type as eqeq (in shared templates)
inline bool %%TYPE_ABBREV%%conjeq(const %%TYPE%%* x, const %%TYPE%%* y, const size_t len, const size_t unused) {
  size_t p;
  for (p = 0; p < len; ++p) {
    if (x[p].r != y[p].r || x[p].i != -y[p].i) return false;
  }
  return true;
}
