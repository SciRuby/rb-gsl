#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector[1, 10, 100, 1000, 10000]
p v
p v.log10

