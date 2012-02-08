// Part of NMatrix, SciRuby Project.
//
// This file provides API for gemm and symmp.

#include "nmatrix_config.h"
#include "smmp.h"

#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_STDBOOL_H
# include <stdbool.h>
#else
typedef char    bool;
# define true    1;
# define false   0;
#endif

