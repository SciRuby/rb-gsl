// Note that this version of transp isn't currently used, as it's about 3x slower than the
// function pointer version (smmp/transp.template.c).
// See also: https://github.com/mohawkjohn/nmatrix/issues/1
//
// This source file is kept in hopes that one day someone can write an optimized generic
// version which beats out the function pointer version.
//

#ifndef SMMP_C
# define SMMP_C

#include "nmatrix.h"

// extern const int nm_sizeof[NM_TYPES+1];

// Generic version of TRANSP from SMMP (sparse matrix multiplication package).
//
// For a template, see transp.template.c.
//
// itype and dtype are the index dtype and the dtype, respectively. Matrix A and matrix B must be the same type if 'move'
// is true.
void transp(y_size_t n, y_size_t m, void* ia, void* ja, bool diaga, void* a,
            void* ib, void* jb, void* b, bool move,
            int8_t itype, int8_t dtype)
{
  y_size_t i, j, max_j, index, y, x = m+1; // temporary variables for storing unsigned values.

  // Clear B
  memset(ib, 0, nm_sizeof[itype]*m);
  if (move) memset(b, 0, nm_sizeof[dtype]*m);

  if (diaga) SetFuncs[itype][Y_SIZE_T](1, ib, 0, &x, 0); // ib[0] = m+1
  else memset(ib, 0, nm_sizeof[itype]); // ib[0] = 0

  // Count indices for each column
  for (i = 0; i < n; ++i) {
    SetFuncs[Y_SIZE_T][itype](1, &max_j, 0, (char*)(ia)+(i+1)*nm_sizeof[itype], 0); // max_j = ia[i+1]
    for (SetFuncs[Y_SIZE_T][itype](1, &j, 0, (char*)(ia)+i*nm_sizeof[itype], 0) /* j = ia[i]*/ ; j < max_j; ++j) {
      SetFuncs[Y_SIZE_T][itype](1, &x, 0, (char*)(ja)+j*nm_sizeof[itype], 0); // x = ja[j]
      Increment[itype]( (char*)(ib) + (x+1)*nm_sizeof[itype] ); // ++(ib[x+1]);
    }
  }

  SetFuncs[Y_SIZE_T][itype](1, &x, 0, ib, 0); // x = ib[i]
  for (i = 0; i < m; ++i) {
    SetFuncs[Y_SIZE_T][itype](1, &y, 0, (char*)(ib)+(i+1)*nm_sizeof[itype], 0); // y = ib[i+1]
    x += y;
    SetFuncs[itype][Y_SIZE_T](1, (char*)(ib)+(i+1)*nm_sizeof[itype], 0, &x, 0); // ib[i+1] = ib[i] + ib[i+1]
  }

  // Now make jb
  for (i = 0; i < n; ++i) {
    SetFuncs[Y_SIZE_T][itype](1, &max_j, 0, (char*)(ia)+(i+1)*nm_sizeof[itype], 0); // max_j = ia[i+1]
    for (SetFuncs[Y_SIZE_T][itype](1, &j, 0, (char*)(ia)+i*nm_sizeof[itype], 0) /* j = ia[i]*/ ; j < max_j; ++j) {
      SetFuncs[Y_SIZE_T][itype](1, &index, 0, (char*)(ja)+j*nm_sizeof[itype], 0); // index = ja[j]
      SetFuncs[Y_SIZE_T][itype](1, &x, 0, (char*)(ib)+index*nm_sizeof[itype], 0); // x = ib[index]
      SetFuncs[itype][Y_SIZE_T](1, (char*)(jb)+x*nm_sizeof[itype], 0, &i, 0); // jb[ib[index]] = i

      if (move) SetFuncs[dtype][dtype](1, (char*)(b)+x*nm_sizeof[dtype], 0, (char*)(a)+j*nm_sizeof[dtype], 0); // jb[ib[index]] = i

      Increment[itype]((char*)(ib)+index*nm_sizeof[itype]); // ++(ib[index]);
    }
  }

  // Now fix up ib
  for (i = m; i >= 1; --i)
    SetFuncs[itype][itype](1, (char*)(ib)+i*nm_sizeof[itype], 0, (char*)(ib)+(i-1)*nm_sizeof[itype], 0); // ib[i] = ib[i-1]

  if (diaga) {
    if (move) {
      j = SMMP_MIN(n,m);

      for (i = 0; i < j; ++i)
        SetFuncs[dtype][dtype](1, (char*)(b)+i*nm_sizeof[dtype], 0, (char*)(a)+i*nm_sizeof[dtype], 0); // b[i] = a[i]
    }
    SetFuncs[itype][Y_SIZE_T](1, ib, 0, &x, 0); // ib[0] = m+1
  } else {
    memset(ib, 0, nm_sizeof[itype]);
  }
}

#endif