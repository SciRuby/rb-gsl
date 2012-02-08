
// Symbolic matrix multiply c=a*b
void %%INT_ABBREV%%_symbmm_(%%INT%% *n, %%INT%% *m, %%INT%% *l, %%INT%% *ia,
	%%INT%% *ja, %%INT%% *diaga, %%INT%% *ib, %%INT%% *jb, %%INT%% *diagb,
	 %%INT%% *ic, %%INT%% *jc, %%INT%% *diagc, %%INT%% *index)
{
  static %%INT%% i, j, k, jj, minlm, minmn, length, maxlmn, istart; /* Local variables */

  /* Parameter adjustments (F2C) */
  --index;
  --jc;
  --ic;
  --jb;
  --ib;
  --ja;
  --ia;

  /* Function Body */
  /* Computing MAX */
  maxlmn = SMMP_MAX_THREE(*l,*m,*n);

  for (i = 1; i <= maxlmn; ++i) // initialize scratch array to 0
    index[i] = 0; /* L10: */

  if (*diagc == 0)  ic[1] = 1;
  else              ic[1] = *n + 2;

  minlm = SMMP_MIN(*l,*m);
  minmn = SMMP_MIN(*m,*n);


  for (i = 1; i <= *n; ++i) { // MAIN LOOP
    istart = -1;
    length = 0;

    for (jj = ia[i]; jj <= ia[i+1]; ++jj) { // merge row lists
      /*    a = d + ... */
      if (jj == ia[i + 1]) {
        if (*diaga == 0 || i > minmn) goto L30; // is this a continue?
        j = i;
      } else {
        j = ja[jj];
      }

      /*    b = d + ... */
      if (index[j] == 0 && *diagb == 1 && j <= minlm) {
        index[j] = istart;
        istart = j;
        ++length;
      }

      for (k = ib[j]; k <= ib[j+1]-1; ++k) {
        if (index[jb[k]] == 0) {
          index[jb[k]] = istart;
          istart = jb[k];
          ++length;
        }
        /* L20: */
      }
L30:
      ;
    }

    /*   row i of jc */
    if (*diagc == 1 && index[i] != 0) --length;

    ic[i+1] = ic[i] + length;

    for (j = ic[i]; j <= ic[i+1]-1; ++j) {
      if (*diagc == 1 && istart == i) {
        istart = index[istart];
        index[i] = 0;
      }
      jc[j] = istart;
      istart = index[istart];
      index[jc[j]] = 0;
      /* L40: */
    }
    index[i] = 0;
  /* L50: */
  }
} /* symbmm_ */





