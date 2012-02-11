
void %%INT_ABBREV%%_symbmm(u_%%INT%% n,   // # rows in A / C
       u_%%INT%% m,   // # columns in A / rows in B
       u_%%INT%% l,   // # columns in B / C
       u_%%INT%%* ia, // the IA array for A
       u_%%INT%%* ja, // the JA array for A
       bool diaga, // 1 for new yale, 0 for old yale
       u_%%INT%%* ib, // the IA array for B
       u_%%INT%%* jb, // the JA array for B
       bool diagb, // 1 for new yale, 0 for old yale
       u_%%INT%%* ic, // the IA array for result
       u_%%INT%%* jc, // the JA array for result
       bool diagc) // 1 for new yale, 0 for old yale
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


  %%INT_ABBREV%%_symbmm_(n, m, l, ia, ja, diaga, ib, jb, diagb, ic, jc, diagc);
}
