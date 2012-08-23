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
// == ruby_symbols.cpp
//
// Ruby symbols used throught the NMatrix project.

/*
 * Standard Includes
 */

#include <ruby.h>

/*
 * Project Includes
 */

/*
 * Macros
 */

/*
 * Global Variables
 */

ID	nm_rb_real,
		nm_rb_imag,

		nm_rb_numer,
		nm_rb_denom,

		nm_rb_complex_conjugate,
		nm_rb_transpose,
		nm_rb_no_transpose,

		nm_rb_dense,
		nm_rb_list,
		nm_rb_yale,

		nm_rb_add,
		nm_rb_sub,
		nm_rb_mul,
		nm_rb_div,

		nm_rb_percent,
		nm_rb_gt,
		nm_rb_lt,
		nm_rb_eql,
		nm_rb_neql,
		nm_rb_gte,
		nm_rb_lte;

VALUE cNMatrix,
      cNMatrix_IO,
      cNMatrix_IO_Matlab,
			cNVector,
			cNMatrix_YaleFunctions,
			cNMatrix_BLAS,
			
			nm_eDataTypeError,
			nm_eStorageTypeError;

/*
 * Forward Declarations
 */

/*
 * Functions
 */

void nm_init_ruby_constants(void) {

	nm_rb_real							= rb_intern("real");
	nm_rb_imag							= rb_intern("imag");

	nm_rb_numer							= rb_intern("numerator");
	nm_rb_denom							= rb_intern("denominator");

	nm_rb_complex_conjugate	= rb_intern("complex_conjugate");
	nm_rb_transpose					= rb_intern("transpose");
	nm_rb_no_transpose			= rb_intern("no_transpose");

	nm_rb_dense 						= rb_intern("dense");
	nm_rb_list							= rb_intern("list");
	nm_rb_yale							= rb_intern("yale");

	nm_rb_add								= rb_intern("+");
	nm_rb_sub								= rb_intern("-");
	nm_rb_mul								= rb_intern("*");
	nm_rb_div								= rb_intern("/");

	nm_rb_percent						= rb_intern("%");
	nm_rb_gt								= rb_intern(">");
	nm_rb_lt								= rb_intern("<");
	nm_rb_eql								= rb_intern("==");
	nm_rb_neql							= rb_intern("!=");
	nm_rb_gte								= rb_intern(">=");
	nm_rb_lte								= rb_intern("<=");
	
}
