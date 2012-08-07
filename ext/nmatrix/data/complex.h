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
// == complex.h
//
// Functions and classes for dealing with complex numbers.

#ifndef COMPLEX_H
#define COMPLEX_H

/*
 * Standard Includes
 */

#include <type_traits>

/*
 * Project Includes
 */

#include "types.h"

/*
 * Macros
 */

/*
 * Types
 */
 
template <typename Type> class Complex;

typedef Complex<float32_t> Complex64;
typedef Complex<float64_t> Complex128;

/*
 * Data
 */
 
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
	inline Complex(Type real = 0, Type imaginary = 0) : r(real), i(imaginary) {}
	
	/*
	 * Copy constructors.
	 */
	inline Complex(const Complex<Type>& other) : r(other.r), i(other.i) {}
	
	/*
	 * Binary operator definitions for varous types.
	 */
	
	////////////////////////////////
	// Complex-Complex Operations //
	////////////////////////////////
	
	template <typename OtherType>
	inline Complex<Type> operator+(const Complex<OtherType>& other) const {
		return Complex<Type>(this->r + other.r, this->i + other.i);
	}
	
	template <typename OtherType>
	inline Complex<Type> operator-(const Complex<OtherType>& other) const {
		return Complex<Type>(this->r - other.r, this->i - other.i);
	}
	
	template <typename OtherType>
	inline Complex<Type> operator*(const Complex<OtherType>& other) const {
		return Complex<Type>(this->r * other.r - this->i * other.i, this->r * other.i - this->i * other.r);
	}

	template <typename OtherType>
	inline Complex<Type> operator/(const Complex<OtherType>& other) const {
		Type new_r, new_i;
		Type denom = this->i * this->i + other.r * other.r;
		
		new_r = (this->r * other.r + this->i * other.i) / denom;
		new_i = (this->r * other.i - this->i * other.r) / denom;
		
		return Complex<Type>(new_r, new_i);
	}
	
	template <typename OtherType>
	inline bool operator<(const Complex<OtherType>& other) const {
		return (this->r < other.r) || ((this->r <= other.r) && (this->i < other.i));
	}
	
	template <typename OtherType>
	inline bool operator>(const Complex<OtherType>& other) const {
		return (this->r > other.r) || ((this->r >= other.r) && (this->i > other.i));
	}
	
	template <typename OtherType>
	inline bool operator==(const Complex<OtherType>& other) const {
		return FP_EQUAL(this->r, other.r) && FP_EQUAL(this->i, other.i);
	}
	
	template <typename OtherType>
	inline bool operator!=(const Complex<OtherType>& other) const {
		return !(*this == other);
	}
	
	template <typename OtherType>
	inline bool operator<=(const Complex<OtherType>& other) const {
		return (*this < other) || (*this == other);
	}
	
	template <typename OtherType>
	inline bool operator>=(const Complex<OtherType>& other) const {
		return (*this > other) || (*this == other);
	}
	
	template <typename OtherType>
	inline operator Complex<OtherType> () {
		return Complex<OtherType>((OtherType)this->r, (OtherType)this->i);
	}
	
	///////////////////////////////
	// Complex-Native Operations //
	///////////////////////////////
	
	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline Complex<Type> operator+(const NativeType& other) const {
		return *this + Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline Complex<Type> operator-(const NativeType& other) const {
		return *this - Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline Complex<Type> operator*(const NativeType& other) const {
		return *this * Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline Complex<Type> operator/(const NativeType& other) const {
		return *this / Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator<(const NativeType& other) const {
		return *this < Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator>(const NativeType& other) const {
		return *this > Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator==(const NativeType& other) const {
		return *this == Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator!=(const NativeType& other) const {
		return *this != Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator<=(const NativeType& other) const {
		return *this <= Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline bool operator>=(const NativeType& other) const {
		return *this >= Complex<Type>(other);
	}

	template <typename NativeType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
	inline operator NativeType () {
		return (NativeType)this->r;
	}
};


///////////////////////////////
// Native-Complex Operations //
///////////////////////////////

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline Complex<ComplexType> operator+(const NativeType& left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) + right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline Complex<ComplexType> operator-(const NativeType& left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) - right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline Complex<ComplexType> operator*(const NativeType& left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) * right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline Complex<ComplexType> operator/(const NativeType& left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) / right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator<(const NativeType left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) < right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator>(const NativeType left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) > right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator==(const NativeType left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) == right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator!=(const NativeType left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) != right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator<=(const NativeType left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) <= right;
}

template <typename NativeType, typename ComplexType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator>=(const NativeType left, const Complex<ComplexType>& right) {
	return Complex<ComplexType>(left) >= right;
}

#endif // COMPLEX_H
