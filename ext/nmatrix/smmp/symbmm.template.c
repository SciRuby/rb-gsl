
// Symbolic matrix multiply c=a*b
void %%INT_ABBREV%%_symbmm_(u_%%INT%% n, u_%%INT%% m, u_%%INT%% l,
    u_%%INT%% *ia, u_%%INT%% *ja, bool diaga,
    u_%%INT%% *ib, u_%%INT%% *jb, bool diagb,
	  u_%%INT%% *ic, u_%%INT%% *jc, bool diagc)
{
  u_%%INT%% index[(1+SMMP_MAX_THREE(l,m,n)) * sizeof(u_%%INT%%)];
  u_%%INT%% i, j, k, jj, minlm, minmn, maxlmn, istart, length; /* Local variables */

  /* Function Body */
  /* Computing MAX */
  maxlmn = SMMP_MAX_THREE(l,m,n);

  for (j = 0; j <= maxlmn; ++j) {
    index[j] = 0;
  }

  if (diagc)  ic[0] = n+1;
  else        ic[0] = 0;

  minlm = SMMP_MIN(l,m);
  minmn = SMMP_MIN(m,n);

  for (i = 0; i < n; ++i) { // MAIN LOOP: through rows
    istart = U%%INT_MAX%%;
    length = 0;

    for (jj = ia[i]; jj <= ia[i+1]; ++jj) { // merge row lists, walking through columns in each row

      // j <- column index given by JA[jj], or handle diagonal.
      if (jj == ia[i+1]) { // Don't really do it the last time -- just handle diagonals in a new yale matrix.
        if (!diaga || i >= minmn) continue;
        j = i;
      } else {
        j = ja[jj];
      }

      if (!index[j+1] && diagb && j < minlm) { // Is there anything at all in column j?
        index[j+1] = istart;
        istart            = j;
        ++length;
      }

      for (k = ib[j]; k <= ib[j+1]-1; ++k) { // Now walk through columns of row J in matrix B.
        if (!index[jb[k]+1]) {
          index[jb[k]+1] = istart;
          istart         = jb[k]; // current column in matrix 2
          ++length;
        }
      }
    }

    if (diagc && index[i+1])
      --length;

    ic[i+1] = ic[i] + length;
    fprintf(stderr, "set row length for row i=%u to %u\n", i, length);

    for (j = ic[i]; j <= ic[i+1]-1; ++j) {
      if (diagc && istart == i) {
        istart                 = index[istart+1];
        index[istart+1] = 0;
      }
      jc[j]                 = istart;

      istart                = index[istart+1];
      index[jc[j]+1] = 0;
    }
    index[i+1] = 0;
  }
} /* symbmm_ */





