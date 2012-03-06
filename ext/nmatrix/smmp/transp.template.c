
void %%INT_ABBREV%%_%%REAL_ABBREV%%_transp_(y_size_t n, y_size_t m, u_%%INT%% *ia, u_%%INT%% *ja, bool diaga, %%REAL%% *a, u_%%INT%% *ib, u_%%INT%% *jb, %%REAL%% *b, bool move)
{
  u_%%INT%% i, j, index;

  // Clear B
  for (i = 0; i < m+1; ++i)
    ib[i] = 0;
  if (move) {
    for (i = 0; i < m+1; ++i)
      b[i] = 0;
  }

  if (diaga) ib[0] = m + 1;
  else       ib[0] = 0;

/*       count indices for each column */

  for (i = 0; i < n; ++i) {
    for (j = ia[i]; j < ia[i+1]; ++j)
        ++(ib[ja[j]+1]);
  }

  for (i = 0; i < m; ++i)
    ib[i+1] = ib[i] + ib[i+1];

/*       now make jb */

  for (i = 0; i < n; ++i) {

    for (j = ia[i]; j < ia[i+1]; ++j) {
      index = ja[j];
      jb[ib[index]] = i;

      if (move)
        b[ib[index]] = a[j];

      ++(ib[index]);
    }
  }

/*       now fixup ib */

  for (i = m; i >= 1; --i)
    ib[i] = ib[i-1];


  if (diaga) {
    if (move) {
      j = SMMP_MIN(n,m);

      for (i = 0; i < j; ++i)
        b[i] = a[i];
    }
    ib[0] = m + 1;

  } else {
    ib[0] = 0;
  }
}
