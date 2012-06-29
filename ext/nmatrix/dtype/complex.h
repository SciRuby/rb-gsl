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

#define COMPLEX_COMPLEX_OPS(other_type) 																																\
	inline Complex<Type> operator+(const Complex<other_type>& other) const {															\
		return Complex<Type>(this->r + other.r, this->i + other.i);																					\
	}																																																			\
																																																				\
	inline Complex<Type> operator-(const Complex<other_type>& other) const {															\
		return Complex<Type>(this->r - other.r, this->i - other.i);																					\
	}																																																			\
																																																				\
	inline Complex<Type> operator*(const Complex<other_type>& other) const {															\
		return Complex<Type>(this->r * other.r - this->i * other.i, this->r * other.i - this->i * other.r);	\
	}																																																			\
																																																				\
	inline Complex<Type> operator/(const Complex<other_type>& other) const {															\
		Type new_r, new_i;																																									\
		Type denom = this->i * this->i + other.r * other.r;																									\
																																																				\
		new_r = (this->r * other.r + this->i * other.i) / denom;																						\
		new_i = (this->r * other.i - this->i * other.r) / denom;																						\
																																																				\
		return Complex<Type>(new_r, new_i);																																	\
	}																																																			\
																																																				\
	inline bool operator<(const Complex<other_type>& other) const {																				\
		return (this->r < other.r) || ((this->r <= other.r) && (this->i < other.i));												\
	}																																																			\
																																																				\
	inline bool operator>(const Complex<other_type>& other) const {																				\
		return (this->r > other.r) || ((this->r >= other.r) && (this->i > other.i));												\
	}																																																			\
																																																				\
	inline bool operator==(const Complex<other_type>& other) const {																			\
		return FP_EQUAL(this->r, other.r) && FP_EQUAL(this->i, other.i);																		\
	}																																																			\
																																																				\
	inline bool operator!=(const Complex<other_type>& other) const {																			\
		return !(*this == other);																																						\
	}																																																			\
																																																				\
	inline bool operator<=(const Complex<other_type>& other) const {																			\
		return (*this < other) || (*this == other);																													\
	}																																																			\
																																																				\
	inline bool operator>=(const Complex<other_type>& other) const {																			\
		return (*this > other) || (*this == other);																													\
	}																																																			\

#define COMPLEX_NATIVE_OPS(type) 																	\
	inline Complex<Type> operator+(const other_type& other) const {	\
		return Complex<Type>(this->r + other, this->i);								\
	}																																\
																																	\
	inline Complex<Type> operator-(const other_type& other) const {	\
		return Complex<Type>(this->r - other, this->i);								\
	}																																\
																																	\
	inline Complex<Type> operator*(const other_type& other) const {	\
		return Complex<Type>(this->r * other, this->i * other);				\
	}																																\
																																	\
	inline Complex<Type> operator/(const other_type& other) const {	\
		Type new_r, new_i;																						\
		Type denom = this->i * this->i + other * other;								\
																																	\
		new_r =  (this->r * other) / denom;														\
		new_i = -(this->i * other) / denom;														\
																																	\
		return Complex<Type>(new_r, new_i);														\
	}																																\
																																	\
	inline bool operator<(const other_type& other) const {					\
		return this->r < other;																				\
	}																																\
																																	\
	inline bool operator>(const other_type& other) const {					\
		return this->r > other.r;																			\
	}																																\
																																	\
	inline bool operator==(const other_type& other) const {					\
		return this->r == other;																			\
	}																																\
																																	\
	inline bool operator!=(const other_type& other) const {					\
		return this->r != other;																			\
	}																																\
																																	\
	inline bool operator<=(const other_type& other) const {					\
		return this->r <= other;																			\
	}																																\
																																	\
	inline bool operator>=(const other_type& other) const {					\
		return this->r >= other;																			\
	}																																\

/*
 * Types
 */
 
template <typename Type> class Complex;

typedef Complex<float>	Complex64;
typedef Complex<double>	Complex128;

/*
 * Classes and Functions
 */

template <typename Type>
class Complex {
	public:
		// The real and immaginary parts of the complex number.
		Type r;
		Type i;
		
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
		 * Binary operator definitions for varous types.
		 */
		
		COMPLEX_COMPLEX_OPS(float)
		COMPLEX_COMPLEX_OPS(double)
		
		COMPLEX_NATIVE_OPS(char)
		COMPLEX_NATIVE_OPS(unsigned char)
		COMPLEX_NATIVE_OPS(short)
		COMPLEX_NATIVE_OPS(int)
		COMPLEX_NATIVE_OPS(long)
		COMPLEX_NATIVE_OPS(float)
		COMPLEX_NATIVE_OPS(double)
};

#endif
