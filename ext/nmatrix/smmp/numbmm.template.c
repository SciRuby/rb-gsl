// numeric matrix multiply c=a*b
void %%INT_ABBREV%%_%%REAL_ABBREV%%_numbmm_(u_%%INT%% n, u_%%INT%% m,
    u_%%INT%% *ia, u_%%INT%% *ja, bool diaga, %%REAL%% *a,
    u_%%INT%% *ib, u_%%INT%% *jb,	bool diagb, %%REAL%% *b,
    u_%%INT%% *ic, u_%%INT%% *jc, bool diagc,	%%REAL%% *c)
{
  u_%%INT%% next[m];
  %%REAL%% sums[m];

  %%REAL%% v;

  u_%%INT%% head, length, temp, ndnz = 0;
  u_%%INT%% jj_start, jj_end, kk_start, kk_end;
  u_%%INT%% i, j, k, kk, jj;
  u_%%INT%% minmn = SMMP_MIN(m,n);

  for (i = 0; i < m; ++i) { // initialize scratch arrays
    next[i] = U%%INT_MAX%%;
    sums[i] = 0;
  }

  for (i = 0; i < n; ++i) { // walk down the rows
    head = U%%INT_MAX%%-1; // head gets assigned as whichever column of B's row j we last visited
    length = 0;

    jj_start = ia[i];
    jj_end   = ia[i+1];

    fprintf(stderr, "*** part A\n");
    for (jj = jj_start; jj <= jj_end; ++jj) { // walk through entries in each row

      if (jj == jj_end) { // if we're in the last entry for this row:
        if (!diaga || i >= minmn) continue;
        j   = i;      // if it's a new Yale matrix, and last entry, get the diagonal position (j) and entry (ajj)
        v   = a[i];
      } else {
        j   = ja[jj]; // if it's not the last entry for this row, get the column (j) and entry (ajj)
        v   = a[jj];
      }

      kk_start = ib[j];   // Find the first entry of row j of matrix B
      kk_end   = ib[j+1];
      for (kk = kk_start; kk <= kk_end; ++kk) {

        if (kk == kk_end) { // Get the column id for that entry
          if (!diagb || j >= minmn) continue;
          k  = j;
          fprintf(stderr, "i=%u, jj=%u, j=%u, DIAG, k=%u, v=%f: assigning partial sum v*b[k]+=%f\n", i, jj, j, k, v, v*b[k]);
          sums[k] += v*b[k];
        } else {
          k  = jb[kk];
          fprintf(stderr, "i=%u, jj=%u, j=%u, kk=%u, k=%u, v=%f: assigning partial sum v*b[kk]+=%f\n", i, jj, j, kk, k, v, v*b[kk]);
          sums[k] += v*b[kk];
        }

        if (next[k] == U%%INT_MAX%%) {
          next[k] = head;
          head    = k;
          ++length;
          fprintf(stderr, "\tk=%u, next[k]=head=%u, new head=%u, length=%u\n", k, next[k], head, length);
        }
      }
    }

/*    if (diagc && i < minln) {
      c[i]    = temp[i];
      temp[i] = 0;
    }

    for (j = ic[i]; j <= ic[i+1]-1; ++j) {
      c[j]        = temp[jc[j]];
      temp[jc[j]] = 0;
    }*/
    fprintf(stderr, "*** part B\n");
    for (jj = 0; jj < length; ++jj) {
      if (sums[head]) {
        if (diagc && head == i) {
          fprintf(stderr, "diag: assigning c[head=%u] <- %f\n", head, sums[head]);
          c[head]    = sums[head];
        } else {
          fprintf(stderr, "nd: assigning c[n+ndnz=%u+%u], head=%u <- %f\n", n,ndnz, head, sums[head]);
          jc[n+1+ndnz] = head;
          c[ n+1+ndnz] = sums[head];
          ++ndnz;
        }
      }

      temp = head;
      head = next[head];

      next[temp] = U%%INT_MAX%%;
      sums[temp] = 0;
    }

    fprintf(stderr, "ic[i+1] at i=%u being assigned %u+%u\n", i, n, ndnz);
    ic[i+1] = n+1+ndnz;
  }
} /* numbmm_ */





