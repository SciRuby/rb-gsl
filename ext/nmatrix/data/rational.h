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
// == rational.h
//
// Functions and classes for dealing with rational numbers.

#ifndef RATIONAL_H
#define RATIONAL_H

/*
 * Standard Includes
 */

#include <type_traits>
#include <ruby.h>

/*
 * Project Includes
 */

#include "types.h"

#include "util/util.h"

/*
 * Macros
 */

/*
 * Types
 */
template <typename Type> class Rational;

typedef Rational<int16_t>	Rational32;
typedef Rational<int32_t>	Rational64;
typedef Rational<int64_t>	Rational128;

/*
 * Data
 */

/*
 * Classes and Functions
 */

template <typename Type>
class Rational {
	public:
	// The numerator and denominator of the rational number.
	Type n;
	Type d;
	
	/*
	 * Default constructor.
	 */
	inline Rational(Type num = 0, Type den = 1) : n(num), d(den) {}
	
	/*
	 * Copy constructors.
	 */
	template <typename OtherType>
	inline Rational(const Rational<OtherType>& other) : n(other.n), d(other.d) {}

	template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
	inline Rational(const Complex<FloatType>& other) : n(0), d(1) {
	  rb_raise(rb_eNotImpError, "cannot convert from complex to rational");
	}

	/*
	 * Binary operator definitions for varous types.
	 */
	
	//////////////////////////////////
	// Rational-Rational Operations //
	//////////////////////////////////
	
	template <typename OtherType>
	inline Rational<Type> operator+(const Rational<OtherType>& other) const {
		long simplify;
		
		Rational<Type> result((this->n * other.d) + (other.n * this->d), this->d * other.d);
		
		simplify = gcf<Type>(result.n, result.d);
		
		result.n /= simplify;
		result.d /= simplify;
		
		return result;
	}

	template <typename OtherType>
	inline Rational<Type>& operator+=(const Rational<OtherType>& other) {
    long simplify;

    this->n = (this->n * other.d) + (other.n * this->d);
    this->d = this->d * other.d;

    simplify = gcf<Type>(this->n, this->d);

    this->n /= simplify;
    this->d /= simplify;

    return *this;
	}
	
	template <typename OtherType>
	inline Rational<Type> operator-(const Rational<OtherType>& other) const {
		long simplify;
		
		Rational<Type> result((this->n * other.d) - (other.n * this->d), this->d * other.d);
		
		simplify = gcf<Type>(result.n, result.d);
		
		result.n /= simplify;
		result.d /= simplify;
		
		return result;
	}
	
	template <typename OtherType>
	inline Rational<Type> operator*(const Rational<OtherType>& other) const {
		int g1 = gcf<Type>(this->n, other.d);
		int g2 = gcf<Type>(this->d, other.n);
		
		return Rational<Type>((this->n / g1) * (other.n / g2), (this->d / g2) * (other.d / g1));
	}


	template <typename OtherType>
	inline Rational<Type>& operator*=(const Rational<OtherType>& other) {
    int g1 = gcf<Type>(this->n, other.d);
    int g2 = gcf<Type>(this->d, other.n);

    this->n = (this->n / g1) * (other.n / g2);
    this->d = (this->d / g2) * (other.d / g1);

    return *this;
	}

	
	template <typename OtherType>
	inline Rational<Type> operator/(const Rational<OtherType>& other) const {
		return *this * Rational<OtherType>(other.d, other.n);
	}
	
	template <typename OtherType>
	inline Rational<Type> operator%(const Rational<OtherType>& other) const {
		long floor_div;
		Rational<Type> prod;
		
		floor_div = (this->n * other.n) / (this->d * other.d);
		prod			= other * Rational<long>(floor_div, 1);
		
		return Rational<long>(this->n, other.n) - prod;
	}
	
	template <typename OtherType>
	inline bool operator<(const Rational<OtherType>& other) const {
		return (this->n * other.d) < (other.n * this->d);
	}
	
	template <typename OtherType>
	inline bool operator>(const Rational<OtherType>& other) const {
		return (this->n * other.d) > (other.n * this->d);
	}
	
	template <typename OtherType>
	inline bool operator==(const Rational<OtherType>& other) const {
		return (this->n == other.n) && (this->d == other.d);
	}
	
	template <typename OtherType>
	inline bool operator!=(const Rational<OtherType>& other) const {
		return !(*this == other);
	}
	
	template <typename OtherType>
	inline bool operator<=(const Rational<OtherType>& other) const {
		return (*this < other) || (*this == other);
	}
	
	template <typename OtherType>
	inline bool operator>=(const Rational<OtherType>& other) const {
		return (*this > other) || (*this == other);
	}
	
	template <typename OtherType>
	inline operator Rational<OtherType> () const {
		return Rational<OtherType>((OtherType)this->n, (OtherType)this->d);
	}
	
	////////////////////////////////
	// Rational-Native Operations //
	////////////////////////////////
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline Rational<Type> operator+(const IntType& other) const {
		return *this + Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline Rational<Type> operator-(const IntType& other) const {
		return *this - Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline Rational<Type> operator*(const IntType& other) const {
		return *this * Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline Rational<Type> operator/(const IntType& other) const {
		return *this / Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline Rational<Type> operator%(const IntType& other) const {
		return *this % Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator<(const IntType& other) const {
		return *this < Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator>(const IntType& other) const {
		return *this > Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator==(const IntType& other) const {
		return *this == Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator!=(const IntType& other) const {
		return *this != Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator<=(const IntType& other) const {
		return *this <= Rational<Type>(other);
	}
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline bool operator>=(const IntType& other) const {
		return *this >= Rational<Type>(other);
	}
	
	template <typename NumType, typename = typename std::enable_if<std::is_arithmetic<NumType>::value>::type>
	inline operator NumType () const {
		return (NumType)this->n / (NumType)this->d;
	}
	
	/*
	 * Special casting operator for Complex numbers.
	 */
	template <typename FloatType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
	inline operator Rational<FloatType> () const {
		return Rational<FloatType>(((FloatType)this->n) / ((FloatType)this->d));
	}
};

////////////////////////////////
// Native-Rational Operations //
////////////////////////////////

/*
 * Integer Math
 */

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline Rational<RationalType> operator+(const IntType& left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) + right;
}

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline Rational<RationalType> operator-(const IntType& left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) - right;
}

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline Rational<RationalType> operator*(const IntType& left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) * right;
}

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline Rational<RationalType> operator/(const IntType& left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) / right;
}

/*
 * Floating Point Math
 */

template <typename FloatType, typename RationalType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline FloatType operator+(const FloatType& left, const Rational<RationalType>& right) {
	return left + (FloatType)right;
}

template <typename FloatType, typename RationalType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline FloatType operator-(const FloatType& left, const Rational<RationalType>& right) {
	return left - (FloatType)right;
}

template <typename FloatType, typename RationalType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline FloatType operator*(const FloatType& left, const Rational<RationalType>& right) {
	return left * (FloatType)right;
}

template <typename FloatType, typename RationalType, typename = typename std::enable_if<std::is_floating_point<FloatType>::value>::type>
inline FloatType operator/(const FloatType& left, const Rational<RationalType>& right) {
	return left / (FloatType)right;
}

/*
 * Comparisons
 */

template <typename NativeType, typename RationalType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator<(const NativeType left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) < right;
	return (left * right.d) < right.n;
}

template <typename NativeType, typename RationalType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator>(const NativeType left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) > right;
	return (left * right.d) > right.n;
}

template <typename IntType, typename RationalType>
inline bool operator==(const typename std::enable_if<std::is_integral<IntType>::value, IntType>::type left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) == right;
	return (left * right.d) == right.n;
}

template <typename FloatType, typename RationalType>
inline bool operator==(const typename std::enable_if<std::is_floating_point<FloatType>::value, FloatType>::type left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) == right;
	return FP_EQUAL(left, ((FloatType)right));
}

template <typename NativeType, typename RationalType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator!=(const NativeType left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) != right;
	return !(left == right);
}

template <typename NativeType, typename RationalType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator<=(const NativeType left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) <= right;
	return (left < right) or (left == right);
}

template <typename NativeType, typename RationalType, typename = typename std::enable_if<std::is_arithmetic<NativeType>::value>::type>
inline bool operator>=(const NativeType left, const Rational<RationalType>& right) {
	//return Rational<RationalType>(left) >= right;
	return (left > right) or (left == right); 
}

#endif // RATIONAL_H
