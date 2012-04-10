// Yale to Bank-Smith
/* Subroutine */ int %%INT_ABBREV%%_%%REAL_ABBREV%%_ytobs_(%%INT%% *n, %%INT%% *ia, %%INT%% *ja, %%INT%% *
	diaga, %%INT%% *syma, %%REAL%% *a, %%INT%% *ib, %%INT%% *jb, %%REAL%% *b,
	%%INT%% *move)
{
    /* System generated locals */
    %%INT%% i__1, i__2, i__3;

    /* Local variables */
    static %%INT%% i__, j, k, ii, jj, lshift;




/*       create the bank-smith data structures b from the */
/*       corresponding yale data structures a */

    /* Parameter adjustments */
    --b;
    --jb;
    --ib;
    --a;
    --ja;
    --ia;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	ib[i__ + 1] = ia[i__ + 1] - ia[i__];
    }

/*       look for upper triangular entries and duplicate entries */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (jj = ia[i__]; jj <= i__2; ++jj) {
	    j = ja[jj];
	    if (i__ == j) {
		--ib[i__ + 1];
		ja[jj] = -j;
	    }
	    if (j > i__) {
		--ib[i__ + 1];
		++ib[j + 1];

/*       check for duplicates */

		i__3 = ia[j + 1] - 1;
		for (k = ia[j]; k <= i__3; ++k) {
		    if (ja[k] == i__) {
			--ib[j + 1];
			ja[jj] = -j;
			goto L30;
		    }
/* L20: */
		}
L30:
		;
	    }
/* L40: */
	}
/* L50: */
    }

/*       compute ib */

    ib[1] = *n + 2;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	ib[i__ + 1] += ib[i__];
    }

/*       initialize b if move = 1 */

    if (*move == 1) {
	lshift = 0;
	if (*syma == 0) {
	    lshift = ib[*n + 1] - ib[1];
	}
	i__1 = ib[*n + 1] + lshift - 1;
	for (ii = 1; ii <= i__1; ++ii) {
/* L62: */
	    b[ii] = 0.f;
	}
	if (*diaga == 1) {
	    i__1 = *n;
	    for (ii = 1; ii <= i__1; ++ii) {
/* L64: */
		b[ii] = a[ii];
	    }
	}
    }

/*       compute jb */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (jj = ia[i__]; jj <= i__2; ++jj) {
	    j = ja[jj];
	    if (j > i__) {
		jb[ib[j]] = i__;
		if (*move == 1) {
		    b[ib[j]] = a[jj];
		}
		++ib[j];
	    } else {
		if (j <= 0) {
		    ja[jj] = -j;
		    if (*move == 1 && i__ == -j) {
			b[i__] = a[jj];
		    }
		} else {
		    jb[ib[i__]] = j;
		    if (*move == 1) {
			b[ib[i__] + lshift] = a[jj];
		    }
		    ++ib[i__];
		}
	    }
/* L70: */
	}
/* L80: */
    }

/*       fixup ib */

    for (i__ = *n; i__ >= 2; --i__) {
/* L90: */
	ib[i__] = ib[i__ - 1];
    }
    ib[1] = *n + 2;
    return 0;
} /* ytobs_ */
