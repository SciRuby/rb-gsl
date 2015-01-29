#!/usr/bin/env ruby
require("gsl")

puts("sign of 5.0")
p GSL::SIGN(5.0)

puts("sign of -2.0")
p GSL::SIGN(-2.0)

puts("sign of 0: positive")
p GSL::SIGN(0)

puts("Is 1 odd?")
p GSL::IS_ODD(1)
p GSL::IS_ODD?(1)
puts("Is 1 even?")
p GSL::IS_EVEN(1)
p GSL::IS_EVEN?(1)

puts("Is 4 odd?")
p GSL::IS_ODD(4)
p GSL::IS_ODD?(4)
puts("Is 4 even?")
p GSL::IS_EVEN(4)
p GSL::IS_EVEN?(4)

puts("frcmp")

p GSL::fcmp(1.0, 1.0, 1e-10)
p GSL::fcmp(1.1, 1.0, 1e-10)
p GSL::fcmp(1.0, 1.1, 1e-10)
