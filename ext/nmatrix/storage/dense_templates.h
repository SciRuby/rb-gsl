/////////////////////////////////////////////////////////////////////
// = NMatrix
//
// A linear algebra library for scientific computation in Ruby.
// NMatrix is part of SciRuby.
//
// NMatrix was originally inspired by and derived from NArray, by
// Masahiro Tanaka: http://narray.rubyforge.org
//
// == Copyright Information
//
// SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
// NMatrix is Copyright (c) 2012, Ruby Science Foundation
//
// Please see LICENSE.txt for additional copyright notices.
//
// == Contributing
//
// By contributing source code to SciRuby, you agree to be bound by
// our Contributor Agreement:
//
// * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
//
// == dense_templates.h
//
// Templates for dense n-dimensional matrix storage.

#ifndef DENSE_TEMPLATES_H
#define DENSE_TEMPLATES_H

/*
 * Standard Includes
 */

#include <stdlib.h>

/*
 * Project Includes
 */

#include "nmatrix.h"
#include "dense.h"

/*
 * Macros
 */

/*
 * Types
 */

/*
 * Data
 */
 

/*
 * Functions
 */

///////////////
// Lifecycle //
///////////////

///////////////
// Accessors //
///////////////

///////////
// Tests //
///////////

#ifdef __cplusplus
extern "C" {
#endif

bool dense_storage_eqeq(const DENSE_STORAGE* left, const DENSE_STORAGE* right);
bool dense_storage_is_symmetric(const DENSE_STORAGE* mat, int lda, bool hermitian);

#ifdef __cplusplus
}
#endif

/////////////
// Utility //
/////////////

/////////////////////////
// Copying and Casting //
/////////////////////////

#endif
