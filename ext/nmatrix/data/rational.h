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
	inline operator Rational<OtherType> () {
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
	
	template <typename IntType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
	inline operator IntType () {
		return (IntType)this->n / (IntType)this->d;
	}
};

////////////////////////////////
// Native-Rational Operations //
////////////////////////////////

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

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator<(const IntType left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) < right;
}

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator>(const IntType left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) > right;
}

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator==(const IntType left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) == right;
}

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator!=(const IntType left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) != right;
}

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator<=(const IntType left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) <= right;
}

template <typename IntType, typename RationalType, typename = typename std::enable_if<std::is_integral<IntType>::value>::type>
inline bool operator>=(const IntType left, const Rational<RationalType>& right) {
	return Rational<RationalType>(left) >= right;
}

#endif
