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

/*
 * Project Includes
 */

#include "types.h"

/*
 * Macros
 */

#define COMPLEX_COMPLEX_OPS(OtherType) 																																	\
	inline Complex<Type> operator+(const Complex<OtherType>& other) const {																\
		return Complex<Type>(this->r + other.r, this->i + other.i);																					\
	}																																																			\
																																																				\
	inline Complex<Type> operator-(const Complex<OtherType>& other) const {																\
		return Complex<Type>(this->r - other.r, this->i - other.i);																					\
	}																																																			\
																																																				\
	inline Complex<Type> operator*(const Complex<OtherType>& other) const {																\
		return Complex<Type>(this->r * other.r - this->i * other.i, this->r * other.i - this->i * other.r);	\
	}																																																			\
																																																				\
	inline Complex<Type> operator/(const Complex<OtherType>& other) const {																\
		Type new_r, new_i;																																									\
		Type denom = this->i * this->i + other.r * other.r;																									\
																																																				\
		new_r = (this->r * other.r + this->i * other.i) / denom;																						\
		new_i = (this->r * other.i - this->i * other.r) / denom;																						\
																																																				\
		return Complex<Type>(new_r, new_i);																																	\
	}																																																			\
																																																				\
	inline bool operator<(const Complex<OtherType>& other) const {																				\
		return (this->r < other.r) || ((this->r <= other.r) && (this->i < other.i));												\
	}																																																			\
																																																				\
	inline bool operator>(const Complex<OtherType>& other) const {																				\
		return (this->r > other.r) || ((this->r >= other.r) && (this->i > other.i));												\
	}																																																			\
																																																				\
	inline bool operator==(const Complex<OtherType>& other) const {																				\
		return FP_EQUAL(this->r, other.r) && FP_EQUAL(this->i, other.i);																		\
	}																																																			\
																																																				\
	inline bool operator!=(const Complex<OtherType>& other) const {																				\
		return !(*this == other);																																						\
	}																																																			\
																																																				\
	inline bool operator<=(const Complex<OtherType>& other) const {																				\
		return (*this < other) || (*this == other);																													\
	}																																																			\
																																																				\
	inline bool operator>=(const Complex<OtherType>& other) const {																				\
		return (*this > other) || (*this == other);																													\
	}																																																			\
																																																				\
	inline operator Complex<OtherType> () {																																\
		return Complex<OtherType>((OtherType)this->r, (OtherType)this->i);																	\
	}

#define COMPLEX_NATIVE_OPS(NativeType)														\
	inline Complex<Type> operator+(const NativeType& other) const {	\
		return Complex<Type>(this->r + other, this->i);								\
	}																																\
																																	\
	inline Complex<Type> operator-(const NativeType& other) const {	\
		return Complex<Type>(this->r - other, this->i);								\
	}																																\
																																	\
	inline Complex<Type> operator*(const NativeType& other) const {	\
		return Complex<Type>(this->r * other, this->i * other);				\
	}																																\
																																	\
	inline Complex<Type> operator/(const NativeType& other) const {	\
		Type new_r, new_i;																						\
		Type denom = this->i * this->i + other * other;								\
																																	\
		new_r =  (this->r * other) / denom;														\
		new_i = -(this->i * other) / denom;														\
																																	\
		return Complex<Type>(new_r, new_i);														\
	}																																\
																																	\
	inline bool operator<(const NativeType& other) const {					\
		return this->r < other;																				\
	}																																\
																																	\
	inline bool operator>(const NativeType& other) const {					\
		return this->r > other;																				\
	}																																\
																																	\
	inline bool operator==(const NativeType& other) const {					\
		return (this->r == other) && FP_IS_ZERO(this->i);							\
	}																																\
																																	\
	inline bool operator!=(const NativeType& other) const {					\
		return !(*this == other);																			\
	}																																\
																																	\
	inline bool operator<=(const NativeType& other) const {					\
		return (*this < other) || (*this == other);										\
	}																																\
																																	\
	inline bool operator>=(const NativeType& other) const {					\
		return (*this > other) || (*this == other);										\
	}																																\
																																	\
	inline operator NativeType () {																	\
		return (NativeType)this->r;																		\
	}

#define NATIVE_COMPLEX_OPS(NativeType, ComplexType)																										\
	inline Complex<ComplexType> operator+(const NativeType& left, const Complex<ComplexType>& right) {	\
		return Complex<ComplexType>(left + right.r, right.i);																							\
	}																																																		\
																																																			\
	inline Complex<ComplexType> operator-(const NativeType& left, const Complex<ComplexType>& right) {	\
		return Complex<ComplexType>(left - right.r, right.i);																							\
	}																																																		\
																																																			\
	inline Complex<ComplexType> operator*(const NativeType& left, const Complex<ComplexType>& right) {	\
		return Complex<ComplexType>(left * right.r, left * right.i);																			\
	}																																																		\
																																																			\
	inline Complex<ComplexType> operator/(const NativeType& left, const Complex<ComplexType>& right) {	\
		ComplexType new_r, new_i;																																					\
		ComplexType denom = right.r * right.r;																														\
																																																			\
		new_r = (left * right.r) / denom;																																	\
		new_i = (left * right.i) / denom;																																	\
																																																			\
		return Complex<ComplexType>(new_r, new_i);																												\
	}																																																		\
																																																			\
	inline bool operator<(const NativeType left, const Complex<ComplexType>& right) {										\
		return left < right.r;																																						\
	}																																																		\
																																																			\
	inline bool operator>(const NativeType left, const Complex<ComplexType>& right) {										\
		return left > right.r;																																						\
	}																																																		\
																																																			\
	inline bool operator==(const NativeType left, const Complex<ComplexType>& right) {									\
		return (left == right.r) && FP_IS_ZERO(right.i);																									\
	}																																																		\
																																																			\
	inline bool operator!=(const NativeType left, const Complex<ComplexType>& right) {									\
		return !(left == right);																																					\
	}																																																		\
																																																			\
	inline bool operator<=(const NativeType left, const Complex<ComplexType>& right) {									\
		return (left < right) || (left == right);																													\
	}																																																		\
																																																			\
	inline bool operator>=(const NativeType left, const Complex<ComplexType>& right) {									\
		return (left > right) || (left == right);																													\
	}

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
		
		COMPLEX_COMPLEX_OPS(float32_t)
		COMPLEX_COMPLEX_OPS(float64_t)
		
		COMPLEX_NATIVE_OPS(uint8_t)
		COMPLEX_NATIVE_OPS(int8_t)
		COMPLEX_NATIVE_OPS(int16_t)
		COMPLEX_NATIVE_OPS(int32_t)
		COMPLEX_NATIVE_OPS(int64_t)
		COMPLEX_NATIVE_OPS(float32_t)
		COMPLEX_NATIVE_OPS(float64_t)
};

NATIVE_COMPLEX_OPS(uint8_t,		float32_t)
NATIVE_COMPLEX_OPS(int8_t,			float32_t)
NATIVE_COMPLEX_OPS(int16_t,			float32_t)
NATIVE_COMPLEX_OPS(int32_t,			float32_t)
NATIVE_COMPLEX_OPS(int64_t,			float32_t)
NATIVE_COMPLEX_OPS(float32_t,		float32_t)
NATIVE_COMPLEX_OPS(float64_t,		float32_t)

NATIVE_COMPLEX_OPS(uint8_t,		float64_t)
NATIVE_COMPLEX_OPS(int8_t,			float64_t)
NATIVE_COMPLEX_OPS(int16_t,			float64_t)
NATIVE_COMPLEX_OPS(int32_t,			float64_t)
NATIVE_COMPLEX_OPS(int64_t,			float64_t)
NATIVE_COMPLEX_OPS(float32_t,		float64_t)
NATIVE_COMPLEX_OPS(float64_t,		float64_t)

#endif
