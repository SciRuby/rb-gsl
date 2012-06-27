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
// Functions and classes for dealing with complex numbers.

#ifndef COMPLEX_H
#define COMPLEX_H

/*
 * Standard Includes
 */

/*
 * Project Includes
 */

/*
 * Macros
 */

/*
 * Types
 */
 
template <typename Type> Complex;

typedef Complex<float> Complex64;
typedef Complex<double> Complex128;

/*
 * Classes and Functions
 */

template <typename Type>
class Complex {
	// The real and immaginary parts of the complex number.
	Type r;
	Type i;
	
	public:
	
		/*
		 * Default constructor.
		 */
		inline Complex(Type r = 0, Type i = 0) {
			this->r = r;
			this->i = i;
		}
		
		/*
		 * Copy constructors.
		 */
		inline Complex(const Complex<Type>& other) {
			this->r = other.r;
			this->i = other.i;
		}
		
		/*
		 * Empty destructor.
		 */
		inline ~Complex(void);
		
		//////////////////////////
		// Arithmetic Operators //
		//////////////////////////
		
		inline Complex<Type> operator+(const Complex<Type>& other) {
			return Complex(this->r + other.r, this->i + other.i);
		}
		
		inline Complex<Type> operator-(const Complex<Type>& other) {
			return Complex(this->r - other.r, this->i - other.i);
		}
		
		inline Complex<Type> operator*(const Complex<Type>& other) {
			return Complex(this->r * other.r - this->i * other.i, this->r * other.i - this->i * other.r);
		}
		
		inline Complex<Type> operator/(const Complex<Type>& other) {
			Type new_r, new_i;
			Type denom = this->i * this->i + other.r * other.r;
			
			new_r = (this->r * other.r + this->i * other.i) / denom;
			new_i = (this->r * other.i - this->i * other.r) / denom;
			
			return Complex(new_r, new_i);
		}
		
		/////////////////
		// Comparisons //
		/////////////////
		
		inline bool operator<(const Complex<Type>& other) {
			return (this->r < other.r) && (this->i < other.i);
		}
		
		inline bool operator>(const Complex<Type>& other) {
			return (this->r > other.r) && (this->i > other.i);
		}
		
		inline bool operator==(const Complex<Type>& other) {
			return FP_EQUAL(this->r, other.r) && FP_EQUAL(this->i, other.i);
		}
		
		inline bool operator!=(const Complex<Type>& other) {
			return !(FP_EQUAL(this->r, other.r) && FP_EQUAL(this->i, other.i););
		}
		
		inline bool operator<=(const Complex<Type>& other) {
			return ((this->r < other.r) && (this->i < other.i)) ||
				(FP_EQUAL(this->r, other.r) && FP_EQUAL(this->i, other.i));
		}
		
		inline bool operator>=(const Complex<Type>& other) {
			return ((this->r > other.r) && (this->i > other.i)) ||
				(FP_EQUAL(this->r, other.r) && FP_EQUAL(this->i, other.i));
		}
}

#endif
