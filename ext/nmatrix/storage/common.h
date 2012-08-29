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
// == common.h
//
// Header file for code common to all storage types.

#ifndef STORAGE_COMMON_H
#define STORAGE_COMMON_H

/*
 * Standard Includes
 */

/*
 * Project Includes
 */
#include "data/data.h"
#include "nmatrix.h"

/*
 * Macros
 */
 
extern "C" {

/*
 * Types
 */

// For binary operations involving matrices that need to be casted.
struct STORAGE_PAIR {
  STORAGE* left;
  STORAGE* right;
};

struct SLICE {
  size_t*	coords; // Coordinate of first element
  size_t*	lengths; // Lengths of slice
  bool  	single; // true if all lengths equal to 1 (represents single matrix element)
};

/*
 * Data
 */

/*
 * Functions
 */

  size_t nm_storage_count_max_elements(const STORAGE* storage);

} // end of extern "C" block

namespace nm {

  /*
   * Templated helper function for element-wise operations, used by dense, yale, and list.
   */
  template <ewop_t op, typename LDType, typename RDType>
  inline LDType ew_op_switch(LDType& left, RDType& right) {
    switch (op) {
      case EW_ADD:
        return left + right;

      case EW_SUB:
        return left - right;

      case EW_MUL:
        return left * right;

      case EW_DIV:
        return left / right;

      case EW_MOD:
        rb_raise(rb_eNotImpError, "Element-wise modulo is currently not supported.");
        break;

      default:
        rb_raise(rb_eStandardError, "this should not happen");
    }
    return 0;
  }

  #define EWOP_NO_DIV_BY_ZERO(ltype, rtype)       template <>       \
  inline ltype ew_op_switch<EW_DIV>( ltype & left, rtype & right) { \
    if (right == 0) rb_raise(rb_eZeroDivError, "cannot divide type by 0, would throw SIGFPE");  \
    return left / right;  \
  }

  // Outlaw certain divisions which would abort SIGFPE
  EWOP_NO_DIV_BY_ZERO(int64_t, int64_t)
  EWOP_NO_DIV_BY_ZERO(int32_t, int32_t)
  EWOP_NO_DIV_BY_ZERO(int32_t, int64_t)
  EWOP_NO_DIV_BY_ZERO(int16_t, int16_t)
  EWOP_NO_DIV_BY_ZERO(int16_t, int32_t)
  EWOP_NO_DIV_BY_ZERO(int16_t, int64_t)
  EWOP_NO_DIV_BY_ZERO(int8_t, int8_t)
  EWOP_NO_DIV_BY_ZERO(int8_t, u_int8_t)
  EWOP_NO_DIV_BY_ZERO(int8_t, int16_t)
  EWOP_NO_DIV_BY_ZERO(int8_t, int32_t)
  EWOP_NO_DIV_BY_ZERO(int8_t, int64_t)
  EWOP_NO_DIV_BY_ZERO(uint8_t, uint8_t)
  EWOP_NO_DIV_BY_ZERO(uint8_t, int8_t)
  EWOP_NO_DIV_BY_ZERO(uint8_t, int16_t)
  EWOP_NO_DIV_BY_ZERO(uint8_t, int32_t)
  EWOP_NO_DIV_BY_ZERO(uint8_t, int64_t)
  EWOP_NO_DIV_BY_ZERO(float, int8_t)
  EWOP_NO_DIV_BY_ZERO(float, u_int8_t)
  EWOP_NO_DIV_BY_ZERO(float, int16_t)
  EWOP_NO_DIV_BY_ZERO(float, int32_t)
  EWOP_NO_DIV_BY_ZERO(float, int64_t)
  EWOP_NO_DIV_BY_ZERO(double, int8_t)
  EWOP_NO_DIV_BY_ZERO(double, u_int8_t)
  EWOP_NO_DIV_BY_ZERO(double, int16_t)
  EWOP_NO_DIV_BY_ZERO(double, int32_t)
  EWOP_NO_DIV_BY_ZERO(double, int64_t)

}

#endif // STORAGE_COMMON_H
