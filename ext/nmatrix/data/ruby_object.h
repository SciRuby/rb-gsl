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
// == dtype.h
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
	VALUE obj_ref;
	
	/*
	 * Default constructor.
	 */
	inline RubyObject(VALUE ref) {
		this->obj_ref = ref;
	}
	
	/*
	 * Copy constructors.
	 */
	inline RubyObject(const RubyObject& other) {
		this->obj_ref = other.obj_ref;
	}
	
	/*
	 * Binary operator definitions.
	 */
	
	inline RubyObject operator+(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->obj_ref, rb_intern("+"), 1, other.obj_ref));
	}
	
	inline RubyObject operator-(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->obj_ref, rb_intern("-"), 1, other.obj_ref));
	}
	
	inline RubyObject operator*(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->obj_ref, rb_intern("*"), 1, other.obj_ref));
	}
	
	inline RubyObject operator/(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->obj_ref, rb_intern("/"), 1, other.obj_ref));
	}
	
	inline RubyObject operator%(const RubyObject& other) const {
		return RubyObject(rb_funcall(this->obj_ref, rb_intern("%"), 1, other.obj_ref));
	}
	
	// FIXME: Add comparison operators.
};

#endif
