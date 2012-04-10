// Copy-transpose. Not used. Instead see dense_transpose_generic in nmatrix.c.
void %%TYPE_ABBREV%%transp(const unsigned int M, const unsigned int N, const %%TYPE%%* A, const int lda, %%TYPE%%* B, const int ldb) {
  unsigned int i, j;

  for (i = 0; i < N; ++i) {
    for (j = 0; j < M; ++j) {
      %%TYPE B[i*ldb+j] = A[j*lda+i]%%
    }
  }
}
