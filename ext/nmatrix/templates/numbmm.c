// numeric matrix multiply c=a*b
void numbmm(const unsigned int n, const unsigned int m,
            const UINT* ia, const UINT* ja, const TYPE* a, const bool diaga,
            const UINT* ib, const UINT* jb, const TYPE* b, const bool diagb,
            UINT* ic, UINT* jc, TYPE* c, const bool diagc) {
  UINT next[m];
  TYPE sums[m];

  TYPE v;

  UINT head, length, temp, ndnz = 0;
  UINT jj_start, jj_end, kk_start, kk_end;
  UINT i, j, k, kk, jj;
  UINT minmn = NM_MIN(m,n);

  for (i = 0; i < m; ++i) { // initialize scratch arrays
    next[i] = -1;
    sums[i] = 0;
  }

  for (i = 0; i < n; ++i) { // walk down the rows
    head = -2; // head gets assigned as whichever column of B's row j we last visited
    length = 0;

    jj_start = ia[i];
    jj_end   = ia[i+1];

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
          sums[k] += v*b[k];
        } else {
          k  = jb[kk];
          sums[k] += v*b[kk];
        }

        if (next[k] == -1) {
          next[k] = head;
          head    = k;
          ++length;
        }
      }
    }

    for (jj = 0; jj < length; ++jj) {
      if (sums[head] != 0) {
        if (diagc && head == i) {
          c[head] = sums[head];
        } else {
          jc[n+1+ndnz] = head;
          c[n+1+ndnz]  = sums[head];
          ++ndnz;
        }
      }

      temp = head;
      head = next[head];

      next[temp] = -1;
      sums[temp] = 0;
    }

    ic[i+1] = n+1+ndnz;
  }

  // Not normally a part of numbmm. Easier to call it here than to create a new 2D function
  // pointer array just for smmp_sort_columns.
  smmp_sort_columns(n, ic, jc, c);
} /* numbmm_ */





