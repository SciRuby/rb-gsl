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
// == data.h
//
// Header file for dealing with data types.

#ifndef RUBY_CONSTANTS_H
#define RUBY_CONSTANTS_H

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
 * Types
 */

/*
 * Data
 */

extern ID rbsym_real,
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

extern VALUE	cNMatrix,
							cNVector,
							cYaleFunctions,
							cBLAS,
			
							nm_eDataTypeError,
							nm_eStorageTypeError;

/*
 * Functions
 */

void nm_init_ruby_constants(void);

#endif // RUBY_CONSTANTS_H
