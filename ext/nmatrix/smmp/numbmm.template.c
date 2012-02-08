
/* Subroutine */ int %%INT_ABBREV%%_%%REAL_ABBREV%%_numbmm_(%%INT%% *n, %%INT%% *m, %%INT%% *l, %%INT%% *ia,
	%%INT%% *ja, %%INT%% *diaga, %%REAL%% *a, %%INT%% *ib, %%INT%% *jb,
	%%INT%% *diagb, %%REAL%% *b, %%INT%% *ic, %%INT%% *jc, %%INT%% *diagc,
	%%REAL%% *c__, %%REAL%% *temp)
{
    /* System generated locals */
    %%INT%% i__1, i__2, i__3;

    /* Local variables */
    static %%INT%% i__, j, k, jj;
    static %%REAL%% ajj;
    static %%INT%% minlm, minln, minmn, maxlmn;




/*       numeric matrix multiply c=a*b */

    /* Parameter adjustments */
    --temp;
    --c__;
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
    i__1 = max(*l,*m);
    maxlmn = max(i__1,*n);
    i__1 = maxlmn;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	temp[i__] = 0.f;
    }
    minlm = min(*l,*m);
    minln = min(*l,*n);
    minmn = min(*m,*n);

/*   c = a*b */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1];
	for (jj = ia[i__]; jj <= i__2; ++jj) {
/*    a = d + ... */
	    if (jj == ia[i__ + 1]) {
		if (*diaga == 0 || i__ > minmn) {
		    goto L30;
		}
		j = i__;
		ajj = a[i__];
	    } else {
		j = ja[jj];
		ajj = a[jj];
	    }
/*    b = d + ... */
	    if (*diagb == 1 && j <= minlm) {
		temp[j] += ajj * b[j];
	    }
	    i__3 = ib[j + 1] - 1;
	    for (k = ib[j]; k <= i__3; ++k) {
/* L20: */
		temp[jb[k]] += ajj * b[k];
	    }
L30:
	    ;
	}
/*    c = d + ... */
	if (*diagc == 1 && i__ <= minln) {
	    c__[i__] = temp[i__];
	    temp[i__] = 0.f;
	}
	i__2 = ic[i__ + 1] - 1;
	for (j = ic[i__]; j <= i__2; ++j) {
	    c__[j] = temp[jc[j]];
/* L40: */
	    temp[jc[j]] = 0.f;
	}
/* L50: */
    }
    return 0;
} /* numbmm_ */





