
void %%INT_ABBREV%%_%%REAL_ABBREV%%_numbmm(y_size_t n, y_size_t m, u_%%INT%% *ia, u_%%INT%% *ja, bool diaga, %%REAL%% *a, u_%%INT%% *ib, u_%%INT%% *jb, bool diagb, %%REAL%% *b, u_%%INT%% *ic, u_%%INT%% *jc, bool diagc, %%REAL%% *c)
{
  if (diaga && ia != ja) {
    fprintf(stderr, "diaga=1, but ia!=ja. For new yale, ia must equal ja.");
    return;
  }

  if (diagb && ib != jb) {
    fprintf(stderr, "diagb=1, but ib!=jb. For new yale, ib must equal jb.");
    return;
  }

  if (diagc && ic != jc) {
    fprintf(stderr, "diagc=1, but ic!=jc. For new yale, ic must equal jc.");
    return;
  }

  %%INT_ABBREV%%_%%REAL_ABBREV%%_numbmm_(n, m, ia, ja, diaga, a, ib, jb, diagb, b, ic, jc, diagc, c);
}

// Perform both the symbolic and numeric steps together.
void %%INT_ABBREV%%_%%REAL_ABBREV%%_smmp(y_size_t n, y_size_t m, u_%%INT%% *ia, u_%%INT%% *ja, bool diaga, %%REAL%% *a, u_%%INT%% *ib, u_%%INT%% *jb, bool diagb, %%REAL%% *b, u_%%INT%% *ic, u_%%INT%% *jc, bool diagc, %%REAL%% *c)
{
  if (diaga && ia != ja) {
    fprintf(stderr, "diaga=1, but ia!=ja. For new yale, ia must equal ja.");
    return;
  }

  if (diagb && ib != jb) {
    fprintf(stderr, "diagb=1, but ib!=jb. For new yale, ib must equal jb.");
    return;
  }

  if (diagc && ic != jc) {
    fprintf(stderr, "diagc=1, but ic!=jc. For new yale, ic must equal jc.");
    return;
  }

  %%INT_ABBREV%%_symbmm_(n, m, ia, ja, diaga, ib, jb, diagb, ic, jc, diagc);
  %%INT_ABBREV%%_%%REAL_ABBREV%%_numbmm_(n, m, ia, ja, diaga, a, ib, jb, diagb, b, ic, jc, diagc, c);
  %%INT_ABBREV%%_%%REAL_ABBREV%%_smmp_sort_columns_(n, ic, jc, c);
}
