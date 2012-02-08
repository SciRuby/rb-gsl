// numeric matrix multiply c=a*b
void %%INT_ABBREV%%_%%REAL_ABBREV%%_numbmm_(%%INT%% *n, %%INT%% *m, %%INT%% *l, %%INT%% *ia,
	%%INT%% *ja, %%INT%% *diaga, %%REAL%% *a, %%INT%% *ib, %%INT%% *jb,
	%%INT%% *diagb, %%REAL%% *b, %%INT%% *ic, %%INT%% *jc, %%INT%% *diagc,
	%%REAL%% *c, %%REAL%% *temp)
{

  /* Local variables */
  static %%INT%% i, j, k, jj;
  static %%REAL%% ajj;
  static %%INT%% minlm, minln, minmn, maxlmn;

  /* Parameter adjustments (F2C) */
  --temp;
  --c;
  --jc;
  --ic;
  --b;
  --jb;
  --ib;
  --a;
  --ja;
  --ia;

  /* Function Body */
  /* Computing MAX */
  maxlmn = SMMP_MAX_THREE(*l,*m,*n);

  for (i = 1; i <= maxlmn; ++i) // initialize scratch array to 0
    temp[i] = (%%REAL%%)(0); /* L10: */

  minlm = SMMP_MIN(*l,*m);
  minln = SMMP_MIN(*l,*n);
  minmn = SMMP_MIN(*m,*n);

/*   c = a*b */

  for (i = 1; i <= *n; ++i) {

    for (jj = ia[i]; jj <= ia[i+1]; ++jj) {
      /*    a = d + ... */
      if (jj == ia[i + 1]) {
        if (*diaga == 0 || i > minmn) goto L30;
        j = i;
        ajj = a[i];
      } else {
        j = ja[jj];
        ajj = a[jj];
      }

      /*    b = d + ... */
      if (*diagb == 1 && j <= minlm) temp[j] += ajj * b[j];

      for (k = ib[j]; k <= ib[j+1]-1; ++k)
        temp[jb[k]] += ajj * b[k]; /* L20: */

L30:
      ;
    }
    /*    c = d + ... */
    if (*diagc == 1 && i <= minln) {
      c[i] = temp[i];
      temp[i] = (%%REAL%%)(0);
    }

    for (j = ic[i]; j <= ic[i+1]-1; ++j) {
      c[j] = temp[jc[j]];
      /* L40: */
      temp[jc[j]] = (%%REAL%%)(0);
    }
    /* L50: */
  }
} /* numbmm_ */





