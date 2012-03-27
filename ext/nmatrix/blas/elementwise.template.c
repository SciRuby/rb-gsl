
int nm_%%TYPE_ABBREV_elementwise(const %%TYPE%%* A, const %%TYPE%%* B, %%TYPE%%* C, size_t n, enum NMatrix_Ops op)
{
  size_t i;

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
  case NM_OP_GT:
    for (i = 0; i < n; ++i) {
      %%TYPE C[i] = A[i] > B[i]%%
    }
    break;
  case NM_OP_LT:
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
  case NM_OP_NOT:
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
