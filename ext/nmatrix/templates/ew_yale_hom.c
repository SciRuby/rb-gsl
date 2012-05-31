int ew_yale_hom(const TYPE* a, const TYPE* b, TYPE* c,
            const u_%%INT%%* ija, const u_%%INT%%* ijb, const u_%%INT%%* ijc,
            unsigned int n, unsigned int m, enum NMatrix_Ops op )
{
  unsigned int i, c_jj, a_jj, b_jj;
  TYPE left, right;

  for (i = 0; i < n; ++i) {
    // do operation on diagonals first:
    if (i < m) {
      if (b) ew_op_binary(op, &(c[i]), a[i], b[i]);
      else   ew_op_unary(op, &(c[i]), a[i]);
    }

    c_jj = ijc[i];
    a_jj = ija[i];
    if (ijb) b_jj = ijb[i];

    while (c_jj < ijc[i+1]) {
      // one or the other has to be the correct column (j = ijc[c_jj] ==? ija[a_jj])
      if (a_jj < ija[i+1] && ija[a_jj] == ijc[c_jj]) {
        left = a[a_jj];
        ++a_jj;
      } else {
        left = 0;
      }

      if (b) { // some ops don't have a second vector
        if (b_jj < ijb[i+1] && ijb[b_jj] == ijc[c_jj]) {
          right = b[b_jj];
          ++b_jj;
        } else {
          right = 0;
        }

        ew_op_binary(op, &(c[c_jj]), left, right);
      } else {
        ew_op_unary(op,  &(c[c_jj]), left);
      }

      ++c_jj;
    }
  }
  return 0;
}
