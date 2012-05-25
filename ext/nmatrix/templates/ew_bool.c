// Boolean element-wise operators
int ew_bool(const TYPE* A, const TYPE* B, bool* C, const int n, enum NMatrix_Ops op) {
  int i;

  switch(op) {
  case '!':
    for (i = 0; i < n; ++i) {
      C[i] = !A[i];
    }
    break;
  case NM_OP_EQEQ:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] == B[i];
    }
    break;
  case NM_OP_NEQ:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] != B[i];
    }
    break;
#RUBY if blueprint.id != :complex
  case '>':
    for (i = 0; i < n; ++i) {
      C[i] = A[i] > B[i];
    }
    break;
  case '<':
    for (i = 0; i < n; ++i) {
      C[i] = A[i] < B[i];
    }
    break;
  case NM_OP_GTE:
    for (i = 0; i < n; ++i) {
      C[i] = A[i] >= B[i];
    }
    break;
  case NM_OP_LTE:
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