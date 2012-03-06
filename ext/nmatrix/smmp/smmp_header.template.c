/* ======================================================================= *
 * Sparse Matrix Multiplication Package                                    *
 * Randolph E. Bank and Craig C. Douglas                                   *
 * na.bank@na-net.ornl.gov and na.cdouglas@na-net.ornl.gov                 *
 * ======================================================================= */
// This was originally derived from the above paper, but the algorithm they
// give, in Fortran, uses 1-based indexing, and I simply could not make it
// work. So I went back to where I found the link to that paper -- SciPy's
// CSR type -- and looked at what they had done:
//
// https://github.com/scipy/scipy/blob/master/scipy/sparse/sparsetools/csr.h
//
// However, the SciPy version does not use the "new Yale" format, but rather
// "old Yale." Thus, some modification was necessary -- reincorporating some
// stuff from the original Bank & Douglas paper.

#include <stdio.h>
#include "smmp.h"
