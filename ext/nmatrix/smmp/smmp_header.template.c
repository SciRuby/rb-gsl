/* ======================================================================= *
 * Sparse Matrix Multiplication Package                                    *
 * Randolph E. Bank and Craig C. Douglas                                   *
 * na.bank@na-net.ornl.gov and na.cdouglas@na-net.ornl.gov                 *
 * ======================================================================= */

/* Translated by f2c (version 20090411) and modified by John Woods.
 * You must link the resulting object file with libf2c:
 *  - on Microsoft Windows system, link with libf2c.lib;
 *  - on Linux or Unix systems, link with .../path/to/libf2c.a -lm
 *  - or, if you install libf2c.a in a standard place, with -lf2c -lm
 *    -- in that order, at the end of the command line, as in
 *
 *		cc *.o -lf2c -lm
 *
 *	The source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,
 *
 *		http://www.netlib.org/f2c/libf2c.zip
 */

#include <stdio.h>
// #include "f2c.h"
#include "smmp.h"