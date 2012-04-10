
// Symbolic matrix multiply c=a*b
void %%TYPE_ABBREV%%_symbmm_(y_size_t n, y_size_t m, YALE_PARAM A, YALE_PARAM B, YALE_PARAM C)
{
  u_%%TYPE%% mask[m];
  u_%%TYPE%% i, j, k, kk, jj, minmn, ndnz = n; /* Local variables */

  u_%%TYPE%% *ia = (u_%%TYPE%%*)(A.ia),
            *ja = (u_%%TYPE%%*)(A.ja),
            *ib = (u_%%TYPE%%*)(B.ia),
            *jb = (u_%%TYPE%%*)(B.ja),
            *ic = (u_%%TYPE%%*)(C.ia);

  for (j = 0; j < m; ++j)
    mask[j] = U%%TYPE_MAX%%;

  if (C.diag)  ic[0] = n+1;
  else        ic[0] = 0;

  minmn = SMMP_MIN(m,n);

  for (i = 0; i < n; ++i) { // MAIN LOOP: through rows

    for (jj = ia[i]; jj <= ia[i+1]; ++jj) { // merge row lists, walking through columns in each row

      // j <- column index given by JA[jj], or handle diagonal.
      if (jj == ia[i+1]) { // Don't really do it the last time -- just handle diagonals in a new yale matrix.
        if (!A.diag || i >= minmn) continue;
        j = i;
      } else j = ja[jj];

      for (kk = ib[j]; kk <= ib[j+1]; ++kk) { // Now walk through columns of row J in matrix B.
        if (kk == ib[j+1]) {
          if (!B.diag || j >= minmn) continue;
          k = j;
        } else k = jb[kk];

        if (mask[k] != i) {
          mask[k] = i;
          ++ndnz;
        }
      }
    }

    if (C.diag && !mask[i]) --ndnz;

    ic[i+1] = ndnz;
  }
} /* symbmm_ */





