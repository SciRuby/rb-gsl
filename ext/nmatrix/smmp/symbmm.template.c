
/* Subroutine */ int %%INT_ABBREV%%_symbmm_(%%INT%% *n, %%INT%% *m, %%INT%% *l, %%INT%% *ia,
	%%INT%% *ja, %%INT%% *diaga, %%INT%% *ib, %%INT%% *jb, %%INT%% *diagb,
	 %%INT%% *ic, %%INT%% *jc, %%INT%% *diagc, %%INT%% *index)
{
    /* System generated locals */
    %%INT%% i__1, i__2, i__3;

    /* Local variables */
    static %%INT%% i__, j, k, jj, minlm, minmn, length, maxlmn, istart;



/*       symbolic matrix multiply c=a*b */

    /* Parameter adjustments */
    --index;
    --jc;
    --ic;
    --jb;
    --ib;
    --ja;
    --ia;

    /* Function Body */
/* Computing MAX */
    i__1 = max(*l,*m);
    maxlmn = max(i__1,*n);
    i__1 = maxlmn;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	index[i__] = 0;
    }
    if (*diagc == 0) {
	ic[1] = 1;
    } else {
	ic[1] = *n + 2;
    }
    minlm = min(*l,*m);
    minmn = min(*m,*n);

/*    main loop */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	istart = -1;
	length = 0;

/*    merge row lists */

	i__2 = ia[i__ + 1];
	for (jj = ia[i__]; jj <= i__2; ++jj) {
/*    a = d + ... */
	    if (jj == ia[i__ + 1]) {
		if (*diaga == 0 || i__ > minmn) {
		    goto L30;
		}
		j = i__;
	    } else {
		j = ja[jj];
	    }
/*    b = d + ... */
	    if (index[j] == 0 && *diagb == 1 && j <= minlm) {
		index[j] = istart;
		istart = j;
		++length;
	    }
	    i__3 = ib[j + 1] - 1;
	    for (k = ib[j]; k <= i__3; ++k) {
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

	if (*diagc == 1 && index[i__] != 0) {
	    --length;
	}
	fprintf(stderr, "smmp:i__=%d\n", i__);
	ic[i__ + 1] = ic[i__] + length;
	i__2 = ic[i__ + 1] - 1;
	for (j = ic[i__]; j <= i__2; ++j) {
	    if (*diagc == 1 && istart == i__) {
		istart = index[istart];
		index[i__] = 0;
	    }
	    jc[j] = istart;
	    istart = index[istart];
	    index[jc[j]] = 0;
/* L40: */
	}
	index[i__] = 0;
/* L50: */
    }
    return 0;
} /* symbmm_ */





