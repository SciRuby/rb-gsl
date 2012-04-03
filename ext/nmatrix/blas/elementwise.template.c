
int nm_d_%%TYPE_ABBREV%%_elementwise(const %%TYPE%%* A, const %%TYPE%%* B, %%TYPE%%* C, size_t n, enum NMatrix_Ops op)
{
  size_t i;
  //fprintf(stderr, "elementwise: n=%d, op=%c\n", n, op);

  switch (op) {
  case '+':
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] + B[i]%%
    }
    break;
  case '-':
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] - B[i]%%
    }
    break;
  case '*':
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] * B[i]%%
    }
    break;
  case '/':
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] / B[i]%%
    }
    break;
  case '%':
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] % B[i]%%
    }
    break;
  case '!':
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = !A[i]%%
    }
    break;
  case NM_OP_NEG:
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = -A[i]%%
    }
    break;
  case NM_OP_EQ:
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] == B[i]%%
    }
    break;
  case '>':
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] > B[i]%%
    }
    break;
  case '<':
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] < B[i]%%
    }
    break;
  case NM_OP_GTE:
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] >= B[i]%%
    }
    break;
  case NM_OP_LTE:
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] <= B[i]%%
    }
    break;
  case '~':
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = ~A[i]%%
    }
    break;
  case '&':
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] & B[i]%%
    }
    break;
  case '|':
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] | B[i]%%
    }
    break;
  case '^':
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] ^ B[i]%%
    }
    break;
  case NM_OP_LSH:
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] << B[i]%%
    }
    break;
  case NM_OP_RSH:
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] >> B[i]%%
    }
    break;
  default:
    rb_raise(rb_eNotImpError, "Unrecognized element-wise operator");
  }
  return 0;
}
