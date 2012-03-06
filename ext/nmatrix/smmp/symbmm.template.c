
// Symbolic matrix multiply c=a*b
void %%INT_ABBREV%%_symbmm_(y_size_t n, y_size_t m, u_%%INT%% *ia, u_%%INT%% *ja, bool diaga, u_%%INT%% *ib, u_%%INT%% *jb, bool diagb, u_%%INT%% *ic, u_%%INT%% *jc, bool diagc)
{
  u_%%INT%% mask[m];
  u_%%INT%% i, j, k, kk, jj, minmn, ndnz = n; /* Local variables */

  for (j = 0; j < m; ++j)
    mask[j] = U%%INT_MAX%%;

  if (diagc)  ic[0] = n+1;
  else        ic[0] = 0;

  minmn = SMMP_MIN(m,n);

  for (i = 0; i < n; ++i) { // MAIN LOOP: through rows

    for (jj = ia[i]; jj <= ia[i+1]; ++jj) { // merge row lists, walking through columns in each row

      // j <- column index given by JA[jj], or handle diagonal.
      if (jj == ia[i+1]) { // Don't really do it the last time -- just handle diagonals in a new yale matrix.
        if (!diaga || i >= minmn) continue;
        j = i;
      } else j = ja[jj];

      for (kk = ib[j]; kk <= ib[j+1]; ++kk) { // Now walk through columns of row J in matrix B.
        if (kk == ib[j+1]) {
          if (!diagb || j >= minmn) continue;
          k = j;
        } else k = jb[kk];

        if (mask[k] != i) {
          mask[k] = i;
          ++ndnz;
        }
      }
    }

    if (diagc && !mask[i]) --ndnz;

    ic[i+1] = ndnz;
  }
} /* symbmm_ */





