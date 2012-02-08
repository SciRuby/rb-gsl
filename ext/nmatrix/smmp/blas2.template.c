
void %%INT_ABBREV%%_%%REAL_ABBREV%%_numbmm(
  %%INT%%  n,    // # rows in A / C
  %%INT%% m,    // # columns in A / rows in B
  %%INT%% l,    // # columns in B / C
  %%INT%% *ia,  // the IA array for A
  %%INT%% *ja,  // the JA array for A
  %%INT%% diaga,// 1 for new yale, 0 for old yale
  %%REAL%% *a,  // the A array for matrix A
  %%INT%% *ib,  // the IA array for B
  %%INT%% *jb,  // the JA array for B
  %%INT%% diagb,// 1 for new yale, 0 for old yale
  %%REAL%% *b,  // the A array for matrix B
  %%INT%% *ic,  // the IA array for result
  %%INT%% *jc,  // the JA array for result
  %%INT%% diagc,// 1 for new yale, 0 for old yale
  %%REAL%% *c,  // the A array for matrix C
  %%REAL%% *temp// scratch vector -- set to NULL if you want automatic allocation
) {
  bool alloc_temp;

  if (!temp) {
    temp = malloc( SMMP_MAX_THREE(l,m,n) * sizeof(%%REAL%%) );
    alloc_temp = true;
  } else alloc_temp = false;

  %%INT_ABBREV%%_%%REAL_ABBREV%%_numbmm_(&n, &m, &l, ia, ja, &diaga, a, ib, jb, &diagb, b, ic, jc, &diagc, c, temp);

  if (alloc_temp) free(temp);
}

// Perform both the symbolic and numeric steps together.
void %%INT_ABBREV%%_%%REAL_ABBREV%%_smmp(
  %%INT%%  n,    // # rows in A / C
  %%INT%% m,    // # columns in A / rows in B
  %%INT%% l,    // # columns in B / C
  %%INT%% *ia,  // the IA array for A
  %%INT%% *ja,  // the JA array for A
  %%INT%% diaga,// 1 for new yale, 0 for old yale
  %%REAL%% *a,  // the A array for matrix A
  %%INT%% *ib,  // the IA array for B
  %%INT%% *jb,  // the JA array for B
  %%INT%% diagb,// 1 for new yale, 0 for old yale
  %%REAL%% *b,  // the A array for matrix B
  %%INT%% *ic,  // the IA array for result
  %%INT%% *jc,  // the JA array for result
  %%INT%% diagc,// 1 for new yale, 0 for old yale
  %%REAL%% *c   // the A array for matrix C
) {
  char scratch[SMMP_MAX_THREE(l,m,n) * SMMP_MAX(sizeof(%%INT%%), sizeof(%%REAL%%))];
  fprintf(stderr, "allocated scratch=%p, size %u\n", scratch, SMMP_MAX_THREE(l,m,n) * SMMP_MAX(sizeof(%%INT%%), sizeof(%%REAL%%)));

  %%INT_ABBREV%%_symbmm_(&n, &m, &l, ia, ja, &diaga, ib, jb, &diagb, ic, jc, &diagc, (%%INT%%*)(&scratch));
  %%INT_ABBREV%%_%%REAL_ABBREV%%_numbmm_(&n, &m, &l, ia, ja, &diaga, a, ib, jb, &diagb, b, ic, jc, &diagc, c, (%%REAL%%*)(&scratch));
}
