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

#include <sstream>
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
std::string matlab_cstring_to_dtype_string(const char* str, size_t bytes) {

  std::ostringstream out;

  for (size_t i = 0; i < bytes / sizeof(MDType); i += sizeof(MDType)) {
    DType val = *reinterpret_cast<const MDType*>(str+i);
    out.write( reinterpret_cast<char*>(&val), sizeof(DType) );
  }

  return out.str();
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
  NM_MATLAB_DTYPE_TEMPLATE_TABLE(ttable, nm::io::matlab_cstring_to_dtype_string, std::string, const char* str, size_t bytes);

  nm::io::matlab_dtype_t from_type = matlab_dtype_from_rbsymbol(from);
  uint8_t to_type;

  if (rb_funcall(options, rb_intern("has_key?"), rb_intern("dtype"))) { // Hash#has_key?(:dtype)
    dtype_t to_dtype = nm_dtype_from_rbsymbol(rb_hash_aref(options, rb_intern("dtype")));
    to_type  = static_cast<int8_t>(to_dtype);
  } else {
    itype_t to_itype = nm_itype_from_rbsymbol(rb_hash_aref(options, rb_intern("itype")));

    // we're going to cheat and use the DTYPE template table. To do this, we just act like uint8_t
    // is a dtype (both are 0, itype and dtype), or we add 1 to the other itypes and treat them as
    // signed.
    to_type = static_cast<uint8_t>(to_itype);
    if (to_itype != UINT8) to_type += 1;
  }

  // For next few lines, see explanation above NM_MATLAB_DTYPE_TEMPLATE_TABLE definition in io.h.
  if (to_type >= static_cast<uint8_t>(COMPLEX64)) {
    rb_raise(rb_eArgError, "can only repack into a simple dtype, no complex/rational/VALUE");
  }

  std::string repacked_data = ttable[to_type][from_type](RSTRING_PTR(str), RSTRING_LEN(str));

  return rb_str_new2(repacked_data.c_str());
}


void nm_init_io() {
  cNMatrix_IO = rb_define_module_under(cNMatrix, "IO");
  cNMatrix_IO_Matlab = rb_define_module_under(cNMatrix_IO, "Matlab");

  rb_define_singleton_method(cNMatrix_IO_Matlab, "repack", (METHOD)nm_rbstring_matlab_repack, 3);
}



}