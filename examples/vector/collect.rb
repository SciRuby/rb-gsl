#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector::Int[0..5]
p v
p v.collect { |a| a*a }
p v
v.collect! { |a| a*a }
p v

