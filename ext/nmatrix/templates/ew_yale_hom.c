int ew_yale_hom(const TYPE* A, const TYPE* B, TYPE* C, const UINT* ija, const UINT* ijb, UINT* ijc, const unsigned int n, const unsigned int m, const enum MathHomOps op)
{
  unsigned int i, c_jj, a_jj, b_jj;
  TYPE left, right;

  for (i = 0; i < n; ++i) {
    // do operation on diagonals first:
    if (i < m) {
      C[i] = MathHomOps[op](A[i], B[i]);
    }

    c_jj = ijc[i];
    a_jj = ija[i];
    if (ijb) b_jj = ijb[i];

    while (c_jj < ijc[i+1]) {
      // one or the other has to be the correct column (j = ijc[c_jj] ==? ija[a_jj])

      if (a_jj < ija[i+1] && ija[a_jj] == ijc[c_jj]) {
        left = A[a_jj];
        ++a_jj;
      } else {
        left = 0;
      }

      if (b_jj < ijb[i+1] && ijb[b_jj] == ijc[c_jj]) {
        right = B[b_jj];
        ++b_jj;
      } else {
        right = 0;
      }

      C[i] = MathHomOps[op](A[i], B[i]);

      ++c_jj;
    }
  }

  return 0;
}
