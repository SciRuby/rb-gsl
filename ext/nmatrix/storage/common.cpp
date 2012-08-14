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
// == common.cpp
//
// Code for the STORAGE struct that is common to all storage types.

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

#include "common.h"

/*
 * Macros
 */

/*
 * Global Variables
 */

/*
 * Forward Declarations
 */

/*
 * Functions
 */

/*
 * Calculate the number of elements in the dense storage structure, based on
 * shape and rank.
 */
size_t storage_count_max_elements(const STORAGE* storage) {
  unsigned int i;
  size_t count = 1;
  
  for (i = storage->rank; i-- > 0;) {
		count *= storage->shape[i];
  }
  
  return count;
}

