// Boolean element-wise operators
int ew_bool(const TYPE* A, const TYPE* B, bool* C, const int n, const enum MathBoolOps op) {
  int i;

  switch(op) {
// This can't work properly because C[i] is bool, so it doesn't know to look in the blueprint for
// the !@ op. Even if it did, it would produce something like C[i] = (struct rational128) { 0, 1 };
// which doesn't make any sense.
// TODO: Fix !@ bang operator for ew_bool.
/*  case NM_MATHOP_BANG:
    for (i = 0; i < n; ++i) {
      C[i] = !A[i];
    }
    break; */
  case NM_MATHOP_EQEQ:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] == B[i];
    }
    break;
  case NM_MATHOP_NEQ:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] != B[i];
    }
    break;
#RUBY if !blueprint.is_a?(CSquare::Generator) && blueprint.id != :complex
  case NM_MATHOP_GT:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] > B[i];
    }
    break;
  case NM_MATHOP_LT:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] < B[i];
    }
    break;
  case NM_MATHOP_GTE:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] >= B[i];
    }
    break;
  case NM_MATHOP_LTE:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] <= B[i];
    }
    break;
#RUBY end
  default:
    rb_raise(rb_eNotImpError, "Unrecognized boolean element-wise operator");
  }
  return 0;
}