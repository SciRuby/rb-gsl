#!/usr/bin/env ruby
require("gsl")
include GSL

puts("e")
p GSL::M_E
p Math::exp(1)


puts("log10(e)")
p GSL::M_LOG10E
p Math::log10(M_E)

puts("sqrt(2)")
p GSL::M_SQRT2
p Math::sqrt(2)

puts("sqrt(1/2)")
p GSL::M_SQRT1_2
p Math::sqrt(1.0/2)

puts("sqrt(3)")
p GSL::M_SQRT3
p Math::sqrt(3)

puts("pi")
p GSL::M_PI
p Math::PI

puts("pi/2")
p GSL::M_PI_2
p Math::PI/2.0

puts("pi/4")
p GSL::M_PI_4
p Math::PI/4.0

puts("sqrt(pi)")
p GSL::M_SQRTPI
p Math::sqrt(Math::PI)

puts("2/sqrt(pi)")
p GSL::M_2_SQRTPI
p 2.0/Math::sqrt(Math::PI)

puts("1/pi")
p GSL::M_1_PI
p 1.0/Math::PI

puts("2/pi")
p GSL::M_2_PI
p 2.0/Math::PI

puts("ln10")
p GSL::M_LN10
p Math::log(10)

puts("ln2")
p GSL::M_LN2
p Math::log(2)

puts("ln(pi)")
p GSL::M_LNPI
p Math::log(Math::PI)

puts("euler constant 0.5772156...")
p GSL::M_EULER
