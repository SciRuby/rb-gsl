=begin
= Mathematical Functions
Contents:
(1) ((<Mathematical Constants|URL:math.html#1>))
(2) ((<Infinities and Not-a-number|URL:math.html#2>))
    (1) ((<Constants|URL:math.html#2.1>))
    (2) ((<Module functions|URL:math.html#2.2>))
(3) ((<Elementary Functions|URL:math.html#3>))
(4) ((<Small Integer Powers|URL:math.html#4>))
(5) ((<Testing the Sign of Numbers|URL:math.html#5>))
(6) ((<Testing for Odd and Even Numbers|URL:math.html#6>))
(7) ((<Maximum and Minimum functions|URL:math.html#7>))
(8) ((<Approximate Comparison of Floating Point Numbers|URL:math.html#8>))

== Mathematical Constants
--- GSL::M_E
    The base of exponentials, e
--- GSL::M_LOG2E
    The base-2 logarithm of e, log_2(e)
--- GSL::M_LOG10E
    The base-10 logarithm of e, log_10(e)
--- GSL::M_SQRT2
    The square root of two, sqrt(2)
--- GSL::M_SQRT1_2
    The square root of one-half, sqrt(1/2)
--- GSL::M_SQRT3
    The square root of three, sqrt(3)
--- GSL::M_PI
    The constant pi
--- GSL::M_PI_2
    Pi divided by two
--- GSL::M_PI_4
    Pi divided by four
--- GSL::M_SQRTPI
    The square root of pi
--- GSL::M_2_SQRTPI
    Two divided by the square root of pi
--- GSL::M_1_PI
    The reciprocal of pi, 1/pi
--- GSL::M_2_PI
    Twice the reciprocal of pi, 2/pi
--- GSL::M_LN10
    The natural logarithm of ten, ln(10)
--- GSL::M_LN2
    The natural logarithm of ten, ln(2)
--- GSL::M_LNPI
    The natural logarithm of ten, ln(pi)
--- GSL::M_EULER
    Euler's constant

== Infinities and Not-a-number

=== Constants
--- GSL::POSINF
    The IEEE representation of positive infinity, 
    computed from the expression +1.0/0.0.
--- GSL::NEGINF
    The IEEE representation of negative infinity, 
    computed from the expression -1.0/0.0.
--- GSL::NAN
    The IEEE representation of the Not-a-Number symbol,
    computed from the ratio 0.0/0.0.

=== Module functions
--- GSL::isnan(x)
    This returns 1 if ((|x|)) is not-a-number.
--- GSL::isnan?(x)
    This returns (({true})) if ((|x|))  is not-a-number, and (({false})) otherwise.
--- GSL::isinf(x)
    This returns +1 if ((|x|))  is positive infinity, 
    -1 if ((|x|))  is negative infinity and 0 otherwise.
--- GSL::isinf?(x)
    This returns (({true})) if ((|x|)) is positive or negative infinity, 
    and (({false})) otherwise.
--- GSL::finite(x)
    This returns 1 if ((|x|)) is a real number, 
    and 0 if it is infinite or not-a-number.
--- GSL::finite?(x)
    This returns (({true})) if ((|x|)) is a real number, 
    and (({false})) if it is infinite or not-a-number.

== Elementary Functions
--- GSL::log1p(x)
    This method computes the value of log(1+x) 
    in a way that is accurate for small ((|x|)). It provides an alternative 
    to the BSD math function log1p(x).
--- GSL::expm1(x)
    This method computes the value of exp(x)-1 
    in a way that is accurate for small ((|x|)). It provides an alternative 
    to the BSD math function expm1(x).
--- GSL::hypot(x, y)
    This method computes the value of sqrt{x^2 + y^2} in a way that 
    avoids overflow.
--- GSL::hypot3(x, y, z) 
    Computes the value of sqrt{x^2 + y^2 + z^2} in a way that avoids overflow. 
--- GSL::acosh(x)
    This method computes the value of arccosh(x). 
--- GSL::asinh(x)
    This method computes the value of arcsinh(x). 
--- GSL::atanh(x)
    This method computes the value of arctanh(x). 

    These methods above can take argument ((|x|)) of
    Integer, Float, Array, Vector or Matrix.

--- GSL::ldexp(x)
    This method computes the value of x * 2^e. 
--- GSL::frexp(x)
    This method splits the number ((|x|)) into its normalized fraction 
    f and exponent e, such that x = f * 2^e and 0.5 <= f < 1. 
    The method returns f and the exponent e as an array, [f, e]. 
    If ((|x|)) is zero, both f and e are set to zero. 

== Small Integer Powers
--- GSL::pow_int(x, n)
    This routine computes the power ((|x^n|)) for integer ((|n|)). 
    The power is computed efficiently -- for example, x^8 is computed as 
    ((x^2)^2)^2, requiring only 3 multiplications. 

--- GSL::pow_2(x)
--- GSL::pow_3(x)
--- GSL::pow_4(x)
--- GSL::pow_5(x)
--- GSL::pow_6(x)
--- GSL::pow_7(x)
--- GSL::pow_8(x)
--- GSL::pow_9(x)
    These methods can be used to compute small integer powers x^2, x^3, etc. 
    efficiently.

== Testing the Sign of Numbers
--- GSL::SIGN(x)
--- GSL::sign(x)
    Return the sign of ((|x|)). 
    It is defined as ((x) >= 0 ? 1 : -1). 
    Note that with this definition the sign of zero is positive 
    (regardless of its IEEE sign bit).

== Testing for Odd and Even Numbers
--- GSL::is_odd(n)
--- GSL::IS_ODD(n)
    Evaluate to 1 if ((|n|)) is odd and 0 if ((|n|)) is even. 
    The argument ((|n|)) must be of Fixnum type.
--- GSL::is_odd?(n)
--- GSL::IS_ODD?(n)
    Return (({true})) if ((|n|)) is odd and (({false})) if even.
--- GSL::is_even(n)
--- GSL::IS_EVEN(n)
    Evaluate to 1 if ((|n|)) is even and 0 if ((|n|)) is odd. 
    The argument ((|n|)) must be of Fixnum type.
--- GSL::is_even?(n)
--- GSL::IS_even?(n)
    Return (({true})) if ((|n|)) is even and (({false})) if odd.

== Maximum and Minimum functions
--- GSL::max(a, b)
--- GSL::MAX(a, b)
--- GSL::min(a, b)
--- GSL::MIN(a, b)
    
== Approximate Comparison of Floating Point Numbers
--- GSL::fcmp(a, b, epsilon = 1e-10)
    This method determines whether ((|x|)) and ((|y|)) are approximately equal to a 
    relative accuracy ((|epsilon|)).
--- GSL::equal?(a, b, epsilon = 1e-10)

== Module Constants
--- GSL::VERSION
    GSL version

--- GSL::RB_GSL_VERSION
--- GSL::RUBY_GSL_VERSION
    Ruby/GSL version

((<prev|URL:ehandling.html>))
((<next|URL:complex.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
