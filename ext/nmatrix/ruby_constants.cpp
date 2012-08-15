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

ID	rbsym_real,
		rbsym_imag,

		rbsym_numer,
		rbsym_denom,

		rbsym_complex_conjugate,
		rbsym_transpose,
		rbsym_no_transpose,

		rbsym_dense,
		rbsym_list,
		rbsym_yale,

		rbsym_add,
		rbsym_sub,
		rbsym_mul,
		rbsym_div,

		rbsym_percent,
		rbsym_gt,
		rbsym_lt,
		rbsym_eql,
		rbsym_neql,
		rbsym_gte,
		rbsym_lte;

VALUE cNMatrix,
			cNVector,
			cYaleFunctions,
			cBLAS,
			
			nm_eDataTypeError,
			nm_eStorageTypeError;

/*
 * Forward Declarations
 */

/*
 * Functions
 */

void Init_ruby_constants(void) {

	rbsym_real							= rb_intern("real");
	rbsym_imag							= rb_intern("imag");

	rbsym_numer							= rb_intern("numerator");
	rbsym_denom							= rb_intern("denominator");

	rbsym_complex_conjugate	= rb_intern("complex_conjugate");
	rbsym_transpose					= rb_intern("transpose");
	rbsym_no_transpose			= rb_intern("no_transpose");

	rbsym_dense 						= rb_intern("dense");
	rbsym_list							= rb_intern("list");
	rbsym_yale							= rb_intern("yale");

	rbsym_add								= rb_intern("+");
	rbsym_sub								= rb_intern("-");
	rbsym_mul								= rb_intern("*");
	rbsym_div								= rb_intern("/");

	rbsym_percent						= rb_intern("%");
	rbsym_gt								= rb_intern(">");
	rbsym_lt								= rb_intern("<");
	rbsym_eql								= rb_intern("==");
	rbsym_neql							= rb_intern("!=");
	rbsym_gte								= rb_intern(">=");
	rbsym_lte								= rb_intern("<=");
	
}

