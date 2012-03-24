// numeric matrix multiply c=a*b
void %%INT_ABBREV%%_%%TYPE_ABBREV%%_numbmm_(y_size_t n, y_size_t m, YALE_PARAM A, YALE_PARAM B, YALE_PARAM C)
{
  u_%%INT%% next[m];
  %%TYPE%% sums[m];

  %%TYPE%% v;
  %%= if [:rational,:complex].include?(dtype.type); "#{dtype.long_dtype.sizeof} temp1;"; end%%

  u_%%INT%% head, length, temp, ndnz = 0;
  u_%%INT%% jj_start, jj_end, kk_start, kk_end;
  u_%%INT%% i, j, k, kk, jj;
  u_%%INT%% minmn = SMMP_MIN(m,n);

  u_%%INT%% *ia = (u_%%INT%%*)(A.ia),
            *ja = (u_%%INT%%*)(A.ja),
            *ib = (u_%%INT%%*)(B.ia),
            *jb = (u_%%INT%%*)(B.ja),
            *ic = (u_%%INT%%*)(C.ia),
            *jc = (u_%%INT%%*)(C.ja);
  %%TYPE%% *a = A.a,
           *b = B.a,
           *c = C.a;

  for (i = 0; i < m; ++i) { // initialize scratch arrays
    next[i] = U%%INT_MAX%%;
    %%TYPE sums[i] = 0%%
  }

  for (i = 0; i < n; ++i) { // walk down the rows
    head = U%%INT_MAX%%-1; // head gets assigned as whichever column of B's row j we last visited
    length = 0;

    jj_start = ia[i];
    jj_end   = ia[i+1];

    for (jj = jj_start; jj <= jj_end; ++jj) { // walk through entries in each row

      if (jj == jj_end) { // if we're in the last entry for this row:
        if (!A.diag || i >= minmn) continue;
        j   = i;      // if it's a new Yale matrix, and last entry, get the diagonal position (j) and entry (ajj)
        %%TYPE v = a[i]%%
      } else {
        j   = ja[jj]; // if it's not the last entry for this row, get the column (j) and entry (ajj)
        %%TYPE v = a[jj]%%
      }

      kk_start = ib[j];   // Find the first entry of row j of matrix B
      kk_end   = ib[j+1];
      for (kk = kk_start; kk <= kk_end; ++kk) {

        if (kk == kk_end) { // Get the column id for that entry
          if (!B.diag || j >= minmn) continue;
          k  = j;
          %%TYPE sums[k] += v*b[k]%%
        } else {
          k  = jb[kk];
          %%TYPE sums[k] += v*b[kk]%%
        }

        if (next[k] == U%%INT_MAX%%) {
          next[k] = head;
          head    = k;
          ++length;
        }
      }
    }

    for (jj = 0; jj < length; ++jj) {
      if (%%TYPE sums[head] != 0%%) {
        if (C.diag && head == i) {
          %%TYPE c[head] = sums[head]%%
        } else {
          jc[n+1+ndnz] = head;
          %%TYPE c[n+1+ndnz] = sums[head]%%
          ++ndnz;
        }
      }

      temp = head;
      head = next[head];

      next[temp] = U%%INT_MAX%%;
      %%TYPE sums[temp] = 0%%
    }

    ic[i+1] = n+1+ndnz;
  }
} /* numbmm_ */





