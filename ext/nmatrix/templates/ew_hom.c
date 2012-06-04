// Homogeneous element-wise operators
int ew_hom(const TYPE* A, const TYPE* B, TYPE* C, const int n, enum MathHomOps op) {
  int i;

  for (i = 0; i < n; ++i) {
    C[i] = MathHomOps[op](A[i], B[i]);
  }

  return 0;
}