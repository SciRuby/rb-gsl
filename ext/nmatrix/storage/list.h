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
// == list.h
//
// List-of-lists n-dimensional matrix storage. Uses singly-linked
// lists.

#ifndef LIST_H
#define LIST_H

// Standard Includes

#include <stdlib.h>

// Project Includes

#include "nmatrix.h"

// Macros

// Types

/* Singly-linked ordered list
 * - holds keys and values
 * - no duplicate keys
 * - keys are ordered
 * - values may be lists themselves
 */
typedef struct l_node {
  size_t key;
  void*  val;
  struct l_node* next;
} NODE;

typedef struct {
  NODE* first;
} LIST;

typedef struct {
	int8_t    dtype;
	size_t    rank;
	size_t*   shape;
	size_t*   offset;
	void*     default_val;
	LIST*     rows;
} LIST_STORAGE;

// Functions

#endif
