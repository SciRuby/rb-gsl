// Bank-Smith to Yale
/* Subroutine */ int %%INT_ABBREV%%_%%REAL_ABBREV%%_bstoy_(%%INT%% *n, %%INT%% *ia, %%INT%% *ja, %%INT%% *
	syma, %%REAL%% *a, %%INT%% *ib, %%INT%% *jb, %%INT%% *diagb, %%REAL%% *b,
	%%INT%% *move)
{
    /* System generated locals */
    %%INT%% i__1, i__2;

    /* Local variables */
    static %%INT%% i__, j, jj, icor, lshift;




/*       create the yale data structures b from the */
/*       corresponding bank-smith data structures a */

/*       compute ib */

    /* Parameter adjustments */
    --b;
    --jb;
    --ib;
    --a;
    --ja;
    --ia;

    /* Function Body */
    if (*diagb == 1) {
	ib[1] = *n + 2;
	icor = 0;
	if (*move == 1) {
	    lshift = 0;
	    if (*syma == 0) {
		lshift = ia[*n + 1] - ia[1];
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L2: */
		b[i__] = a[i__];
	    }
	}
    } else {
	ib[1] = 1;
	icor = 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	ib[i__ + 1] = ia[i__ + 1] - ia[i__] + icor;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__]; j <= i__2; ++j) {
	    ++ib[ja[j] + 1];
/* L20: */
	}
/* L30: */
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	ib[i__ + 1] += ib[i__];
    }
    if (*diagb == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jb[ib[i__]] = i__;
	    if (*move == 1) {
		b[ib[i__]] = a[i__];
	    }
/* L45: */
	    ++ib[i__];
	}
    }

/*       now compute jb */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (jj = ia[i__]; jj <= i__2; ++jj) {
	    j = ja[jj];
	    jb[ib[j]] = i__;
	    jb[ib[i__]] = j;
	    if (*move == 1) {
		b[ib[j]] = a[jj];
		b[ib[i__]] = a[jj + lshift];
	    }
	    ++ib[i__];
	    ++ib[j];
/* L50: */
	}
/* L60: */
    }

/*       fixup ib */

    for (i__ = *n; i__ >= 2; --i__) {
/* L70: */
	ib[i__] = ib[i__ - 1];
    }
    if (*diagb == 1) {
	ib[1] = *n + 2;
    } else {
	ib[1] = 1;
    }
    return 0;
} /* bstoy_ */
