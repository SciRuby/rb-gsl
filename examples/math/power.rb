#!/usr/bin/env ruby
require("gsl")

p GSL::pow_int(8, 3)
p 8*8*8.0

p GSL::pow_2(2)
p GSL::pow_3(2)
p GSL::pow_4(2)
p GSL::pow_5(2)
p GSL::pow_6(2)
p GSL::pow_7(2)
p GSL::pow_8(2)
p GSL::pow_9(2)

p GSL::pow(3.2, 4.5)
p GSL::pow([1, 2, 3, 4, 5], 3)
p GSL::pow(GSL::Vector[1, 2, 3, 4, 5], 3)
