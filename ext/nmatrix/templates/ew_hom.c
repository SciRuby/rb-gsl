// Homogeneous element-wise operators
int ew_hom(const TYPE* A, const TYPE* B, TYPE* C, const int n, enum MathHomOps op) {
  int i;

#RUBY if [:complex,:rational].include?(blueprint.id)

  for (i = 0; i < n; ++i) {
    C[i] = MathHomOps[op](A[i], B[i]);
  }

#RUBY else

  switch(op) {
  case NM_MATHOP_ADD:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] + B[i];
    }
    break;
  case NM_MATHOP_SUB:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] - B[i];
    }
    break;
  case NM_MATHOP_MUL:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] * B[i];
    }
    break;
  case NM_MATHOP_DIV:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] / B[i];
    }
    break;
  case NM_MATHOP_MOD:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] % B[i];
    }
    break;
  default:
    rb_raise(rb_eNotImpError, "Unrecognized homogeneous element-wise operator");
  }

#RUBY end

  return 0;
}