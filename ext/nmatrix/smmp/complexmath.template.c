

inline static %%TYPE%% %%TYPE_ABBREV%%_mul(double ar, double ai, double br, double bi) {
  %%TYPE%% result;

  result.r = ar * br - ai * bi;
  result.i = ar * bi + br * ai;

  return result;
}


inline static %%TYPE%% %%TYPE_ABBREV%%_add(double ar, double ai, double br, double bi) {
  %%TYPE%% result;

  result.r = ar + br;
  result.i = ai + bi;

  return result;
}
