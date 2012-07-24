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
// == ruby_object.h
//
// Functions and classes for dealing with Ruby objects.

#ifndef RUBY_OBJECT_H
#define RUBY_OBJECT_H

/*
 * Standard Includes
 */

#include <ruby.h>

/*
 * Project Includes
 */

#include "ruby_symbols.h"

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
 * Classes and Functions
 */

class RubyObject {
	public:
	VALUE rval;
	
	/*
	 * Default constructor.
	 */
	inline RubyObject(VALUE ref) {
		this->rval = ref;
	}
	
	/*
	 * Copy constructors.
	 */
	inline RubyObject(const RubyObject& other) {
		this->rval = other.rval;
	}
	
	/*
	 * Binary operator definitions.
	 */
	
	inline RubyObject operator+(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->rval, rbsym_add, 1, other.rval));
	}
	
	inline RubyObject operator-(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->rval, rbsym_sub, 1, other.rval));
	}
	
	inline RubyObject operator*(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->rval, rbsym_mul, 1, other.rval));
	}
	
	inline RubyObject operator/(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->rval, rbsym_div, 1, other.rval));
	}
	
	inline RubyObject operator%(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->rval, rbsym_percent, 1, other.rval));
	}
	
	inline bool operator>(const RubyObject& other) const {
		return rb_funcall(this->rval, rbsym_gt, 1, other.rval) == Qtrue;
	}
	
	inline bool operator<(const RubyObject& other) const {
		return rb_funcall(this->rval, rbsym_lt, 1, other.rval) == Qtrue;
	}
	
	inline bool operator==(const RubyObject& other) const {
		return rb_funcall(this->rval, rbsym_eql, 1, other.rval) == Qtrue;
	}
	
	inline bool operator!=(const RubyObject& other) const {
		return rb_funcall(this->rval, rbsym_neql, 1, other.rval) == Qtrue;
	}
	
	inline bool operator>=(const RubyObject& other) const {
		return rb_funcall(this->rval, rbsym_gte, 1, other.rval) == Qtrue;
	}
	
	inline bool operator<=(const RubyObject& other) const {
		return rb_funcall(this->rval, rbsym_lte, 1, other.rval) == Qtrue;
	}
};

#endif
