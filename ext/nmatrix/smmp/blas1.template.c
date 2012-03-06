
void %%INT_ABBREV%%_symbmm(y_size_t n, y_size_t m, YALE_PARAM A, YALE_PARAM B, YALE_PARAM C)
{
  if (A.diag && A.ia != A.ja) {
    fprintf(stderr, "A.diag=true, but ia!=ja. For new yale, ia must equal ja.");
    return;
  }

  if (B.diag && B.ia != B.ja) {
    fprintf(stderr, "B.diag=true, but ia!=ja. For new yale, ia must equal ja.");
    return;
  }

  if (C.diag && C.ia != C.ja) {
    fprintf(stderr, "C.diag=true, but ia!=ja. For new yale, ia must equal ja.");
    return;
  }


  %%INT_ABBREV%%_symbmm_(n, m, A, B, C);
}
