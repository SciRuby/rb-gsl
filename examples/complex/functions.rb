#!/usr/bin/env ruby
require("gsl")
include Math

# Test GSL::Complex#sqrt
z = GSL::Complex.alloc(1, 2)
sqrtz = z.sqrt
p sqrtz
d = sqrtz*sqrtz
p z == d
sqrtz = GSL::Complex.sqrt(z)
d = sqrtz*sqrtz
p z == d

# Test GSL::Complex::sqrt_real
a = -2.0
c = GSL::Complex.sqrt_real(a)
p GSL::equal?(a, (c*c).re)

# Test GSL::Complex#exp, log
p f = z.exp
p GSL::Complex.exp(z)
p f.log
p f.log == z

p z.pow(f)  # [1.931e-02 1.752e-02], verified with Octave result

# The results below are verified with the Octave results
p GSL::Complex.sin(z)
p GSL::Complex.cos(z)
p GSL::Complex.tan(z)
p GSL::Complex.sec(z)
p GSL::Complex.csc(z)
p GSL::Complex.cot(z)
p GSL::Complex.arcsin(z)
p GSL::Complex.arcsin_real(2)
p GSL::Complex.arccos(z)
p GSL::Complex.arccos_real(2)
p GSL::Complex.arctan(z)
p GSL::Complex.arcsec(z)
p GSL::Complex.arcsec_real(2)
p GSL::Complex.arccsc(z)
p GSL::Complex.arccsc_real(2)
p GSL::Complex.arccot(z)

p GSL::Complex.sinh(z)
p GSL::Complex.cosh(z)
p GSL::Complex.tanh(z)
p GSL::Complex.sech(z)
p GSL::Complex.csch(z)
p GSL::Complex.coth(z)

p GSL::Complex.arcsinh(z)
p GSL::Complex.arccosh(z)
p GSL::Complex.arccosh_real(2)
p GSL::Complex.arctanh(z)
p GSL::Complex.arctanh_real(2)
p GSL::Complex.arcsech(z)
p GSL::Complex.arccsch(z)
p GSL::Complex.arccoth(z)

p 1/z
p z.inverse
p z.conjugate
p z.negative
p -z
z2 = GSL::Complex.alloc(2, 3)
p z*z2
p z/z2
p z - z2
p z + z2
p z + 1
p 2 + z
p 2 * z
p 2/z
p z/2

