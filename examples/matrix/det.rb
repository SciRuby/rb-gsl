#!/usr/bin/env ruby
require("gsl")
p m = GSL::Matrix[[1, 2, 3, 4, 5, 6, 7, 8 ,0], 3, 3]
puts("m.det = #{m.det}")
puts("m.trace = #{m.trace}")

p m.to_complex.det
p m.to_complex.trace

