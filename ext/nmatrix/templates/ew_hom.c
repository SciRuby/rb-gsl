// Homogeneous element-wise operators
int ew_hom(const TYPE* A, const TYPE* B, TYPE* C, const int n, enum NMatrix_Ops op) {
  int i;

  switch(op) {
  case '+':
    for (i = 0; i < n; ++i) {
      C[i] = A[i] + B[i];
    }
    break;
  case '-':
    for (i = 0; i < n; ++i) {
      C[i] = A[i] - B[i];
    }
    break;
  case '*':
    for (i = 0; i < n; ++i) {
      C[i] = A[i] * B[i];
    }
    break;
  case '/':
    for (i = 0; i < n; ++i) {
      C[i] = A[i] / B[i];
    }
    break;
  case '%':
    for (i = 0; i < n; ++i) {
      C[i] = A[i] % B[i];
    }
    break;
  default:
    rb_raise(rb_eNotImpError, "Unrecognized homogeneous element-wise operator");
  }
  return 0;
}