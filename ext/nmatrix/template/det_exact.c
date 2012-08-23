
// Calculate the determinant for a matrix A of size 2 or 3. Store the result in B.
//
// Recall that lda is the distance between row start positions in the linear arrangement of the matrix A.
void det_exact(const int M, const TYPE* A, const int lda, TYPE* B) {
  LONG_TYPE x, y; // temporaries representing determinants of submatrices (if 3x3)

  if (M == 2) {
    *B = A[0] * A[lda+1] - A[1] * A[lda];

  } else if (M == 3) {
    x = A[lda+1] * A[2*lda+2] - A[lda+2] * A[2*lda+1]; // ei - fh
    y = A[lda] * A[2*lda+2] -   A[lda+2] * A[2*lda];   // fg - di
    x = A[0]*x - A[1]*y ; // a*(ei-fh) - b*(fg-di)

    y = A[lda] * A[2*lda+1] - A[lda+1] * A[2*lda];    // dh - eg
    *B = A[2]*y + x; // c*(dh-eg) + _

  } else if (M < 2) {
    rb_raise(rb_eArgError, "can only calculate exact determinant of a square matrix of size 2 or larger");
  } else {
    rb_raise(rb_eNotImpError, "exact determinant calculation needed for matrices larger than 3x3");
  }
}
