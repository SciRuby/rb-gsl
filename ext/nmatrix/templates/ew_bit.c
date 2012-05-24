// Bitwise element-wise operators
int ew_bit(const TYPE* A, const TYPE* B, TYPE* C, const int n, enum NMatrix_Ops op) {
  int i;

  switch(op) {
  case '~':
    for (i = 0; i < n; ++i) {
      C[i] = ~A[i];
    }
    break;
  case '&':
    for (i = 0; i < n; ++i) {
      C[i] = A[i] & B[i];
    }
    break;
  case '|':
    for (i = 0; i < n; ++i) {
      C[i] = A[i] | B[i];
    }
    break;
  case '^':
    for (i = 0; i < n; ++i) {
      C[i] = A[i] ^ B[i];
    }
    break;
  case NM_OP_LSH:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] << B[i];
    }
    break;
  case NM_OP_RSH:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] >> B[i];
    }
    break;
  default:
    rb_raise(rb_eNotImpError, "Unrecognized bitwise element-wise operator");
  }
  return 0;
}