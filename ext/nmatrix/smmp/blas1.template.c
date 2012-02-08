
void %%INT_ABBREV%%_symbmm(%%INT%% n,   // # rows in A / C
       %%INT%% m,   // # columns in A / rows in B
       %%INT%% l,   // # columns in B / C
       %%INT%%* ia, // the IA array for A
       %%INT%%* ja, // the JA array for A
       %%INT%% diaga, // 1 for new yale, 0 for old yale
       %%INT%%* ib, // the IA array for B
       %%INT%%* jb, // the JA array for B
       %%INT%% diagb, // 1 for new yale, 0 for old yale
       %%INT%%* ic, // the IA array for result
       %%INT%%* jc, // the JA array for result
       %%INT%% diagc, // 1 for new yale, 0 for old yale
       %%INT%%* index // scratch vector -- set to NULL if you want automatic allocation
       )
{
  bool alloc_index;

  // Do we need to allocate a scratch vector? Make it so it automatically frees when we finish.
  if (!index) {
    fprintf(stderr, "allocating index to size %d * %d\n", SMMP_MAX_THREE(l,m,m), sizeof(%%INT%%));
    index = malloc( SMMP_MAX_THREE(l,m,n) * sizeof(%%INT%%) );
    alloc_index = true;
  } else alloc_index = false;

  %%INT_ABBREV%%_symbmm_(&n, &m, &l, ia, ja, &diaga, ib, jb, &diagb, ic, jc, &diagc, index);

  if (alloc_index) free(index);
}
