#!/usr/bin/env ruby
require("gsl")

puts("gsl_log1p")
p GSL::log1p(0.1)
p Math::log(1 + 0.1)

puts("gsl_expm1")
p GSL::expm1(0.1)
p Math::exp(0.1) - 1.0

puts("gsl_hypot")
p GSL::hypot(20000, 30000)
p Math::sqrt(20000*20000 + 30000*30000)

puts("gsl_acosh")
p GSL::acosh(1.5)
p Math::acosh(1.5)

puts("gsl_asinh")
p GSL::asinh(1.5)
p Math::asinh(1.5)

puts("gsl_atanh")
p GSL::atanh(0.5)
p Math::atanh(0.5)

puts("gsl_ldexp")
p GSL::ldexp(1.5, 8)
p Math::ldexp(1.5, 8)

puts("gsl_frexp")
p GSL::frexp(100)
p Math::frexp(100)

