// numeric matrix multiply c=a*b
void %%INT_ABBREV%%_%%REAL_ABBREV%%_numbmm_(u_%%INT%% n, u_%%INT%% m, u_%%INT%% l,
    u_%%INT%% *ia, u_%%INT%% *ja, bool diaga, %%REAL%% *a,
    u_%%INT%% *ib, u_%%INT%% *jb,	bool diagb, %%REAL%% *b,
    u_%%INT%% *ic, u_%%INT%% *jc, bool diagc,	%%REAL%% *c)
{
  %%REAL%% temp[SMMP_MAX_THREE(l,m,n) * sizeof(%%REAL%%)];

  /* Local variables */
  u_%%INT%% i, j, k, jj;
  %%REAL%% ajj;
  u_%%INT%% minlm, minln, minmn, maxlmn;

  maxlmn = SMMP_MAX_THREE(l,m,n);

  for (i = 0; i < maxlmn; ++i) // initialize scratch array to 0
    temp[i] = 0;

  minlm = SMMP_MIN(l,m);
  minln = SMMP_MIN(l,n);
  minmn = SMMP_MIN(m,n);


  for (i = 0; i < n; ++i) { // walk down the rows

    for (jj = ia[i]; jj <= ia[i+1]; ++jj) { // walk through columns in each row

      if (jj == ia[i+1]) { // if we're in the last entry for this row:
        if (!diaga || i >= minmn) continue;
        j   = i;      // if it's a new Yale matrix, and last entry, get the diagonal position (j) and entry (ajj)
        ajj = a[i];
      } else {
        j   = ja[jj]; // if it's not the last entry for this row, get the column (j) and entry (ajj)
        ajj = a[jj];
      }

      if (diagb && j < minlm) temp[j] += ajj * b[j];

      for (k = ib[j]; k <= ib[j+1]-1; ++k) // walk through columns in row j of B
        temp[jb[k]] += ajj * b[k];
    }

    if (diagc && i < minln) {
      c[i]    = temp[i];
      temp[i] = 0;
    }

    for (j = ic[i]; j <= ic[i+1]-1; ++j) {
      c[j]        = temp[jc[j]];
      temp[jc[j]] = 0;
    }
  }
} /* numbmm_ */





