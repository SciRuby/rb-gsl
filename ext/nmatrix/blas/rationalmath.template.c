

%%TYPE%% %%TYPE_ABBREV%%_muldiv(int64_t anum, int64_t aden, int64_t bnum, int64_t bden, char k) {
  %%TYPE%% result;
  int64_t t, g1, g2;

  if (k == '/') { // Switch numerator and denominator for division (and move sign)
    if (bnum < 0) {
      anum = -anum;
      bnum = -bnum;
    }
    t = bnum;
    bnum = bden;
    bden = t;
  }

  g1 = nmrb_gcd(anum, bden);
  g2 = nmrb_gcd(aden, bnum);

  result.n = (anum / g1) * (bnum / g2);
  result.d = (aden / g2) * (bden / g1);

  return result;
}


%%TYPE%% %%TYPE_ABBREV%%_addsub(int64_t anum, int64_t aden, int64_t bnum, int64_t bden, char k) {
  %%TYPE%% result;

  int64_t ig = nmrb_gcd(aden, bden);
  int64_t a  = anum * (bden / ig);
  int64_t b  = bnum * (aden / ig);
  int64_t c;

  if (k == '+') c=a+b;
  else          c=a-b;

  b        = aden / ig;
  ig       = nmrb_gcd(aden, ig);
  result.n = c / ig;
  a        = bden / ig;
  result.d = a*b;

  return result;
}


%%TYPE%% %%TYPE_ABBREV%%_left_shift(%%= dtype == :rational128 ? "int64_t anum, int64_t aden, int64_t bnum" : (dtype == :rational64 ? "int32_t anum, int32_t aden, int32_t bnum" : "int16_t anum, int16_t aden, int16_t bnum")%%)
{
  %%TYPE%% c;
  c.n = anum << bnum;
  c.d = aden << bnum; // assume its an integer cast to a rational
  return c;
}


%%TYPE%% %%TYPE_ABBREV%%_right_shift(%%= dtype == :rational128 ? "int64_t anum, int64_t aden, int64_t bnum" : (dtype == :rational64 ? "int32_t anum, int32_t aden, int32_t bnum" : "int16_t anum, int16_t aden, int16_t bnum")%%)
{
  %%TYPE%% c;
  c.n = anum >> bnum;
  c.d = aden >> bnum; // assume its an integer cast to a rational
  return c;
}
