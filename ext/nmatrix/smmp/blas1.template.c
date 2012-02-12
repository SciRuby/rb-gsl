
void %%INT_ABBREV%%_symbmm(u_%%INT%% n, u_%%INT%% m, u_%%INT%% l, u_%%INT%%* ia, u_%%INT%%* ja, bool diaga, u_%%INT%%* ib, u_%%INT%%* jb, bool diagb, u_%%INT%%* ic, u_%%INT%%* jc, bool diagc)
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
}
