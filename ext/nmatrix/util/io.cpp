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
// == io.cpp
//
// Input/output support functions.

#include "io.h"

#include <ruby.h>

namespace nm { namespace io {

  const char* const MATLAB_DTYPE_NAMES[NUM_MATLAB_DTYPES] = {
    "miUNDEFINED0",
    "miINT8",
    "miUINT8",
    "miINT16",
    "miUINT16",
    "miINT32",
    "miUINT32",
    "miSINGLE",
    "miRESERVED8",
    "miDOUBLE",
    "miRESERVED10",
    "miRESERVED11",
    "miINT64",
    "miUINT64",
    "miMATRIX"
  };

  const size_t MATLAB_DTYPE_SIZES[NUM_MATLAB_DTYPES] = {
    1, // undefined
    1, // int8
    1, // uint8
    2, // int16
    2, // uint16
    4, // int32
    4, // uint32
    sizeof(float),
    1, // reserved
    sizeof(double),
    1, // reserved
    1, // reserved
    8, // int64
    8, // uint64
    1  // matlab array?
  };


/*
 * Templated function for converting from MATLAB dtypes to NMatrix dtypes.
 */
template <typename DType, typename MDType>
char* matlab_cstring_to_dtype_string(size_t& result_len, const char* str, size_t bytes) {

  result_len   = bytes / sizeof(DType);
  char* result = ALLOC_N(char, bytes / sizeof(DType));

  if (bytes % sizeof(MDType) != 0) {
    rb_raise(rb_eArgError, "the given string does not divide evenly for the given MATLAB dtype");
  }

  for (size_t i = 0, j = 0; i < bytes; i += sizeof(MDType), j += sizeof(DType)) {
    *reinterpret_cast<DType*>(result+j) = (DType)(*reinterpret_cast<const MDType*>(str + i));
  }

  return result;
}

}} // end of namespace nm::io

extern "C" {

///////////////////////
// Utility Functions //
///////////////////////

/*
 * Converts a string to a data type.
 */
dtype_t nm_dtype_from_rbstring(VALUE str) {

  for (size_t index = 0; index < NM_NUM_DTYPES; ++index) {
  	if (!std::strncmp(RSTRING_PTR(str), DTYPE_NAMES[index], RSTRING_LEN(str))) {
  		return static_cast<dtype_t>(index);
  	}
  }

  rb_raise(rb_eArgError, "Invalid data type string specified.");
}


/*
 * Converts a symbol to a data type.
 */
dtype_t nm_dtype_from_rbsymbol(VALUE sym) {

  for (size_t index = 0; index < NM_NUM_DTYPES; ++index) {
    if (SYM2ID(sym) == rb_intern(DTYPE_NAMES[index])) {
    	return static_cast<dtype_t>(index);
    }
  }

  rb_raise(rb_eArgError, "Invalid data type symbol specified.");
}


/*
 * Converts a symbol to an index type.
 */
itype_t nm_itype_from_rbsymbol(VALUE sym) {

  for (size_t index = 0; index < NM_NUM_ITYPES; ++index) {
    if (SYM2ID(sym) == rb_intern(ITYPE_NAMES[index])) {
    	return static_cast<itype_t>(index);
    }
  }

  rb_raise(rb_eArgError, "Invalid index type specified.");
}

/*
 * Converts a string to a storage type. Only looks at the first three
 * characters.
 */
stype_t nm_stype_from_rbstring(VALUE str) {

  for (size_t index = 0; index < NM_NUM_STYPES; ++index) {
    if (!std::strncmp(RSTRING_PTR(str), STYPE_NAMES[index], 3)) {
    	return static_cast<stype_t>(index);
    }
  }

  rb_raise(rb_eArgError, "Invalid storage type string specified");
  return DENSE_STORE;
}

/*
 * Converts a symbol to a storage type.
 */
stype_t nm_stype_from_rbsymbol(VALUE sym) {

  for (size_t index = 0; index < NM_NUM_STYPES; ++index) {
    if (SYM2ID(sym) == rb_intern(STYPE_NAMES[index])) {
    	return static_cast<stype_t>(index);
    }
  }

  rb_raise(rb_eArgError, "Invalid storage type symbol specified");
  return DENSE_STORE;
}


/*
 * Converts a MATLAB data-type symbol to an enum.
 */
static nm::io::matlab_dtype_t matlab_dtype_from_rbsymbol(VALUE sym) {
  for (size_t index = 0; index < nm::io::NUM_MATLAB_DTYPES; ++index) {
    if (SYM2ID(sym) == rb_intern(nm::io::MATLAB_DTYPE_NAMES[index])) {
    	return static_cast<nm::io::matlab_dtype_t>(index);
    }
  }

  rb_raise(rb_eArgError, "Invalid matlab type specified.");
}


/*
 * Take a string of bytes which represent MATLAB data type values and repack them into a string
 * of bytes representing values of an NMatrix dtype (or itype).
 *
 * Returns what appears to be a Ruby String.
 *
 * Arguments:
 * * str        :: the data
 * * from       :: symbol representing MATLAB data type (e.g., :miINT8)
 * * options    :: hash, give either :itype => some_itype or :dtype => some_dtype, tells function
 *                 what to give as output.
 */
static VALUE nm_rbstring_matlab_repack(VALUE self, VALUE str, VALUE from, VALUE options) {
  nm::io::matlab_dtype_t from_type = matlab_dtype_from_rbsymbol(from);
  uint8_t to_type;

  if (TYPE(options) != T_HASH) {
    rb_raise(rb_eArgError, "third argument to repack must be an options hash");
  }

  if (RB_HASH_HAS_SYMBOL_KEY(options, "dtype")) { // Hash#has_key?(:dtype)

    dtype_t to_dtype  = nm_dtype_from_rbsymbol(rb_hash_aref(options, ID2SYM(rb_intern("dtype"))));
    to_type           = static_cast<int8_t>(to_dtype);

  } else if (RB_HASH_HAS_SYMBOL_KEY(options, "itype")) {

    itype_t to_itype  = nm_itype_from_rbsymbol(rb_hash_aref(options, ID2SYM(rb_intern("itype"))));

    // we're going to cheat and use the DTYPE template table. To do this, we just act like uint8_t
    // is a dtype (both are 0, itype and dtype), or we add 1 to the other itypes and treat them as
    // signed.
    to_type           = static_cast<uint8_t>(to_itype);
    if (to_itype != UINT8) to_type += 1;


  } else {
    rb_raise(rb_eArgError, "third argument must have either :itype or :dtype as a key");
  }

  // For next few lines, see explanation above NM_MATLAB_DTYPE_TEMPLATE_TABLE definition in io.h.
  if (to_type >= static_cast<uint8_t>(COMPLEX64)) {
    rb_raise(rb_eArgError, "can only repack into a simple dtype, no complex/rational/VALUE");
  }

  // Do the actual repacking -- really simple!
  NM_MATLAB_DTYPE_TEMPLATE_TABLE(ttable, nm::io::matlab_cstring_to_dtype_string, char*, size_t& result_len, const char* str, size_t bytes);

  size_t repacked_data_length;
  char* repacked_data = ttable[to_type][from_type](repacked_data_length, RSTRING_PTR(str), RSTRING_LEN(str));

  // Encode as 8-bit ASCII with a length -- don't want to hiccup on \0
  VALUE result = rb_str_new(repacked_data, repacked_data_length);
  free(repacked_data); // Don't forget to free what we allocated!

  return result;
}


/*
 * Take two byte-strings (real and imaginary) and treat them as if they contain
 * a sequence of data of type dtype. Merge them together and return a new string.
 */
static VALUE nm_rbstring_merge(VALUE self, VALUE rb_real, VALUE rb_imaginary, VALUE rb_dtype) {

  // Sanity check.
  if (RSTRING_LEN(rb_real) != RSTRING_LEN(rb_imaginary)) {
    rb_raise(rb_eArgError, "real and imaginary components do not have same length");
  }

  dtype_t dtype = nm_dtype_from_rbsymbol(rb_dtype);
  size_t len    = DTYPE_SIZES[dtype];

  char *real    = RSTRING_PTR(rb_real),
       *imag    = RSTRING_PTR(rb_imaginary);

  char* merge   = ALLOCA_N(char, RSTRING_LEN(rb_real)*2);

  size_t merge_pos = 0;

  // Merge the two sequences
  for (size_t i = 0; i < RSTRING_LEN(rb_real); i += len) {

    // Copy real number
    memcpy(merge + merge_pos, real + i, len);
    merge_pos += len;

    // Copy imaginary number
    memcpy(merge + merge_pos, imag + i, len);
    merge_pos += len;
  }

  return rb_str_new(merge, merge_pos);
}


void nm_init_io() {
  cNMatrix_IO = rb_define_module_under(cNMatrix, "IO");
  cNMatrix_IO_Matlab = rb_define_module_under(cNMatrix_IO, "Matlab");

  rb_define_singleton_method(cNMatrix_IO_Matlab, "repack", (METHOD)nm_rbstring_matlab_repack, 3);
  rb_define_singleton_method(cNMatrix_IO_Matlab, "complex_merge", (METHOD)nm_rbstring_merge, 3);
}



}