
/* Subroutine */ void %%INT_ABBREV%%_%%REAL_ABBREV%%_transp_(%%INT%% *n, %%INT%% *m, %%INT%% *ia, %%INT%% *ja,
	 %%INT%% *diaga, %%REAL%% *a, %%INT%% *ib, %%INT%% *jb, %%REAL%% *b, %%INT%% *
	move)
{
    /* System generated locals */
    %%INT%% i__1, i__2;

    /* Local variables */
    static %%INT%% i__, j, index;




/*       compute b = a(transpose) */

/*       first make ib */

    /* Parameter adjustments */
    --b;
    --jb;
    --ib;
    --a;
    --ja;
    --ia;

    /* Function Body */
    i__1 = *m + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	ib[i__] = 0;
    }
    if (*move == 1) {
	i__1 = *m + 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L15: */
	    b[i__] = 0.f;
	}
    }
    if (*diaga == 1) {
	ib[1] = *m + 2;
    } else {
	ib[1] = 1;
    }

/*       count indices for each column */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__]; j <= i__2; ++j) {
	    ++ib[ja[j] + 1];
/* L20: */
	}
/* L30: */
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	ib[i__ + 1] = ib[i__] + ib[i__ + 1];
    }

/*       now make jb */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (j = ia[i__]; j <= i__2; ++j) {
	    index = ja[j];
	    jb[ib[index]] = i__;
	    if (*move == 1) {
		b[ib[index]] = a[j];
	    }
	    ++ib[index];
/* L50: */
	}
/* L60: */
    }

/*       now fixup ib */

    for (i__ = *m; i__ >= 2; --i__) {
/* L70: */
	ib[i__] = ib[i__ - 1];
    }
    if (*diaga == 1) {
	if (*move == 1) {
	    j = SMMP_MIN(*n,*m);
	    i__1 = j;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L80: */
		b[i__] = a[i__];
	    }
	}
	ib[1] = *m + 2;
    } else {
	ib[1] = 1;
    }
} /* transp_ */
