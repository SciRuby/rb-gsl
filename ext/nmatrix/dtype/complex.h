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
		return this->r > other.r;																			\
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
	inline Complex<ComplexType>& operator+(const NativeType& left, const Complex<ComplexType>& right) {	\
		return Complex<ComplexType>(left + right.r, right.i);																							\
	}																																																		\
																																																			\
	inline Complex<ComplexType>& operator-(const NativeType& left, const Complex<ComplexType>& right) {	\
		return Complex<ComplexType>(left - right.r, right.i);																							\
	}																																																		\
																																																			\
	inline Complex<ComplexType>& operator*(const NativeType& left, const Complex<ComplexType>& right) {	\
		return Complex<ComplexType>(left * right.r, left * right.i);																			\
	}																																																		\
																																																			\
	inline Complex<ComplexType>& operator/(const NativeType& left, const Complex<ComplexType>& right) {	\
		ComplexType new_r, new_i;																																					\
		ComplexType denom = right.r * right.r;																														\
																																																			\
		new_r = (left * right.r) / denom;																																	\
		new_i = (left * right.i) / denom;																																	\
																																																			\
		return Complex<Type>(new_r, new_i);																																\
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
	}																																																		\
																																																			\
	inline operator Complex<ComplexType> (NativeType nv) {																							\
		return Complex<ComplexType>(nv, 0);																																\
	}

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

NATIVE_COMPLEX_OPS(char,					float)
NATIVE_COMPLEX_OPS(unsigned char,	float)
NATIVE_COMPLEX_OPS(short,					float)
NATIVE_COMPLEX_OPS(int,						float)
NATIVE_COMPLEX_OPS(long,					float)
NATIVE_COMPLEX_OPS(float,					float)
NATIVE_COMPLEX_OPS(double,				float)

NATIVE_COMPLEX_OPS(char,					double)
NATIVE_COMPLEX_OPS(unsigned char,	double)
NATIVE_COMPLEX_OPS(short,					double)
NATIVE_COMPLEX_OPS(int,						double)
NATIVE_COMPLEX_OPS(long,					double)
NATIVE_COMPLEX_OPS(float,					double)
NATIVE_COMPLEX_OPS(double,				double)

#endif
