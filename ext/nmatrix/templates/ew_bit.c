// Bitwise element-wise operators
int ew_bit(const TYPE* A, const TYPE* B, TYPE* C, const int n, enum NMatrix_Ops op) {
  int i;

  switch(op) {
  case NM_MATHOP_NOT:
    for (i = 0; i < n; ++i) {
      C[i] = ~A[i];
    }
    break;
  case NM_MATHOP_AND:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] & B[i];
    }
    break;
  case NM_MATHOP_OR:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] | B[i];
    }
    break;
  case NM_MATHOP_XOR:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] ^ B[i];
    }
    break;
  case NM_MATHOP_LSH:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] << B[i];
    }
    break;
  case NM_MATHOP_RSH:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] >> B[i];
    }
    break;
  default:
    rb_raise(rb_eNotImpError, "Unrecognized bitwise element-wise operator");
  }
  return 0;
}