=begin
= Complex Numbers
Contents:
(1) ((<Class methods|URL:complex.html#1>))
(2) ((<Properties of Complex Numbers|URL:complex.html#2>))
(3) ((<Complex Arithmetic Operators|URL:complex.html#3>))
(4) ((<Elementary Complex Functions|URL:complex.html#4>))
(5) ((<Complex Trigonometric Functions|URL:complex.html#5>))
(6) ((<Inverse Complex Trigonometric Functions|URL:complex.html#6>))
(7) ((<Complex Hyperbolic Functions|URL:complex.html#7>))
(8) ((<Inverse Complex Hyperbolic Functions|URL:complex.html#8>))

== Class Methods
--- GSL::Complex.alloc(re, im)
--- GSL::Complex.rect(re, im)
--- GSL::Complex[re, im]
    These create a GSL::Complex object with real and imaginary part ((|re, im|)).

--- GSL::Complex.polar(r, theta)
    This returns a GSL::Complex object in polar representation, with the amplitude ((|r|)) and the phase (argument) ((|theta|)).

== Properties of complex numbers
--- GSL::Complex#real
--- GSL::Complex#re
--- GSL::Complex#REAL
    Returns the real part

--- GSL::Complex#imag
--- GSL::Complex#im
--- GSL::Complex#IMAG
    Returns the imaginary part

--- GSL::Complex#set(re, im)
--- GSL::Complex#set_complex(re, im)
--- GSL::Complex#SET_COMPLEX(re, im)
    Set the real and imaginary parts of the complex number.

--- GSL::Complex#set_real(re)
--- GSL::Complex#set_re(re)
--- GSL::Complex#SET_REAL(re)
--- GSL::Complex#real=(re)
--- GSL::Complex#re=(re)
--- GSL::Complex#set_imag(im)
--- GSL::Complex#set_im(im)
--- GSL::Complex#SET_IMAG(im)
--- GSL::Complex#imag=(im)
--- GSL::Complex#im=(im)
    Set the real or imaginary parts of the complex number.

--- GSL::Complex#arg
    Returns the argument

--- GSL::Complex#abs, abs2, logabs
    Returns the magnitude, squared magnitude, and the logarithm of the magnitude

== Complex arithmetic operators
--- GSL::Complex#add(b)
--- GSL::Complex#+(b)
    Return the sum of the complex numbers ((|self|)) and ((|b|)). 
--- GSL::Complex#sub(b)
--- GSL::Complex#-(b)
    Return the difference of the complex numbers ((|self|)) and ((|b|)). 
--- GSL::Complex#mul(b)
--- GSL::Complex#*(b)
    Returns the product of the complex numbers  ((|self|)) and ((|b|)).
--- GSL::Complex#div(b)
--- GSL::Complex#/(b)
    Returns the quotient of the complex numbers  ((|self|)) and ((|b|)).

--- GSL::Complex#add_real
--- GSL::Complex#sub_real
--- GSL::Complex#mul_real
--- GSL::Complex#div_real
--- GSL::Complex#add_imag
--- GSL::Complex#sub_imag
--- GSL::Complex#mul_imag
--- GSL::Complex#div_imag

--- GSL::Complex#conjugate
    Returns the complex conjugate of the complex number ((|self|)).
--- GSL::Complex#inverse
    Returns the inverse of the complex number ((|self|)).
--- GSL::Complex#negative
    Returns the negative of the complex number ((|self|)).

== Elementary Complex Functions
--- GSL::Complex#sqrt
--- GSL::Complex#pow(az)
--- GSL::Complex#pow_real(a)
--- GSL::Complex#exp
--- GSL::Complex#log
--- GSL::Complex#log10
--- GSL::Complex#log_b(b)

--- GSL::Complex.sqrt(z)
--- GSL::Complex.sqrt_real(a)
--- GSL::Complex.pow(z, za)
--- GSL::Complex.pow_real(z, a)
--- GSL::Complex.exp(z)
--- GSL::Complex.log(z)
--- GSL::Complex.log10(z)
--- GSL::Complex.log_b(z, b)

== Complex Trigonometric Functions
--- GSL::Complex#sin
--- GSL::Complex#cos
--- GSL::Complex#tan
--- GSL::Complex#sec
--- GSL::Complex#csc
--- GSL::Complex#cot

--- GSL::Complex.sin(z)
--- GSL::Complex.cos(z)
--- GSL::Complex.tan(z)
--- GSL::Complex.sec(z)
--- GSL::Complex.csc(z)
--- GSL::Complex.cot(z)

== Inverse Complex Trigonometric Functions
--- GSL::Complex#arcsin
--- GSL::Complex#arccos
--- GSL::Complex#arctan
--- GSL::Complex#arcsec
--- GSL::Complex#arccsc
--- GSL::Complex#arccot

--- GSL::Complex.arcsin(z)
--- GSL::Complex.arcsin_real(a)
--- GSL::Complex.arccos(z)
--- GSL::Complex.arccos_real(a)
--- GSL::Complex.arctan(z)
--- GSL::Complex.arcsec(z)
--- GSL::Complex.arcsec_real(a)
--- GSL::Complex.arccsc(z)
--- GSL::Complex.arccsc_real(z)
--- GSL::Complex.arccot(z)

== Complex Hyperbolic Functions
--- GSL::Complex#sinh
--- GSL::Complex#cosh
--- GSL::Complex#tanh
--- GSL::Complex#sech
--- GSL::Complex#csch
--- GSL::Complex#coth

--- GSL::Complex.sinh(z)
--- GSL::Complex.cosh(z)
--- GSL::Complex.tanh(z)
--- GSL::Complex.sech(z)
--- GSL::Complex.csch(z)
--- GSL::Complex.coth(z)

== Inverse Complex Hyperbolic Functions
--- GSL::Complex#arcsinh
--- GSL::Complex#arccosh
--- GSL::Complex#arctanh
--- GSL::Complex#arcsech
--- GSL::Complex#arccsch
--- GSL::Complex#arccoth

--- GSL::Complex#arcsinh(z)
--- GSL::Complex#arccosh(z)
--- GSL::Complex#arccosh_real(a)
--- GSL::Complex#arctanh(z)
--- GSL::Complex#arctanh_real(z)
--- GSL::Complex#arcsech(z)
--- GSL::Complex#arccsch(z)
--- GSL::Complex#arccoth(z)

((<prev|URL:math.html>))
((<next|URL:poly.html>))

((<Reference index|URL:ref.html>))
((<top|URL:index.html>))

=end
