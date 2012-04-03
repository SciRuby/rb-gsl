
int %%INT_ABBREV%%_%%TYPE_ABBREV%%_ew(y_size_t n, y_size_t m, enum NMatrix_Ops op, const u_%%INT%%* ija, const u_%%INT%%* ijb, const u_%%INT%%* ijc, %%TYPE%%* a, %%TYPE%%* b, %%TYPE%%* c)
{
  y_size_t i;
  y_size_t c_jj, a_jj, b_jj;
  %%TYPE%% left, right;

  for (i = 0; i < n; ++i) {
    // do operation on diagonals first:
    if (i < m) {
      if (b) %%TYPE_ABBREV%%_ew_op_binary(op, &(c[i]), a[i], b[i]);
      else   %%TYPE_ABBREV%%_ew_op_unary(op, &(c[i]), a[i]);
    }

    c_jj = ijc[i];
    a_jj = ija[i];
    if (ijb) b_jj = ijb[i];

    //fprintf(stderr, "i16_f64_ew: n=%d, i=%d, c_jj=%d\n", (int)(n), (int)(i), (int)(c_jj));
    while (c_jj < ijc[i+1]) {
      // one or the other has to be the correct column (j = ijc[c_jj] ==? ija[a_jj])
      if (a_jj < ija[i+1] && ija[a_jj] == ijc[c_jj]) {
        %%TYPE left = a[a_jj]%%
        ++a_jj;
      } else {
        %%TYPE left = 0%%
      }

      if (b) { // some ops don't have a second vector
        if (b_jj < ijb[i+1] && ijb[b_jj] == ijc[c_jj]) {
          %%TYPE right = b[b_jj]%%
          ++b_jj;
        } else {
          %%TYPE right = 0%%
        }

        %%TYPE_ABBREV%%_ew_op_binary(op, &(c[c_jj]), left, right);
      } else {
        %%TYPE_ABBREV%%_ew_op_unary(op,  &(c[c_jj]), left);
      }

      ++c_jj;
    }
  }
  return 0;
}
