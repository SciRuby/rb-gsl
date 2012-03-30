
inline static %%TYPE%% BOOL2%%= dtype.id.to_s.upcase%%(bool expr) {
  %%TYPE%% result;
  result.r = expr;
  result.i = 0;
  return result;
}



inline static %%TYPE%% %%TYPE_ABBREV%%_mul(double ar, double ai, double br, double bi) {
  %%TYPE%% result;

  result.r = ar * br - ai * bi;
  result.i = ar * bi + br * ai;

  return result;
}


inline static %%TYPE%% %%TYPE_ABBREV%%_div(double ar, double ai, double br, double bi) {
  %%TYPE%% result;
  double denom = br * br + bi * bi;

  result.r = (ar * br + ai * bi) / denom;
  result.i = (ai * br - ar * bi) / denom;

  return result;
}


inline static %%TYPE%% %%TYPE_ABBREV%%_add(%%= dtype == :complex128 ? "double ar, double ai, double br, double bi" : "float ar, float ai, float br, float bi"%%) {
  %%TYPE%% result;

  result.r = ar + br;
  result.i = ai + bi;

  return result;
}


inline static %%TYPE%% %%TYPE_ABBREV%%_sub(%%= dtype == :complex128 ? "double ar, double ai, double br, double bi" : "float ar, float ai, float br, float bi"%%) {
  %%TYPE%% result;

  result.r = ar - br;
  result.i = ai - bi;

  return result;
}


inline static %%TYPE%% %%TYPE_ABBREV%%_mod(%%= dtype == :complex128 ? "double ar, double ai, double br, double bi" : "float ar, float ai, float br, float bi"%%) {
  %%TYPE%% result;
  rb_raise(rb_eNotImpError, "modulo arithmetic for complex numbers not yet implemented");

  return result;
}


inline static %%TYPE%% %%TYPE_ABBREV%%_bang(%%= dtype == :complex128 ? "double ar, double ai" : "float ar, float ai"%%) {
  %%TYPE%% result = {!ar, 0};
  return result;
}

inline static %%TYPE%% %%TYPE_ABBREV%%_negate(%%= dtype == :complex128 ? "double ar, double ai" : "float ar, float ai"%%) {
  %%TYPE%% result;
  result.r = -ar;
  result.i = -ai;
  return result;
}

