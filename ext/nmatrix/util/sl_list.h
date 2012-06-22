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
// == sl_list.h
//
// Singly-linked list implementation

#ifndef SL_LIST_H
#define SL_LIST_H

// Standard Includes

#include <stdlib.h>

// Project Includes

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

// Functions

#endif
