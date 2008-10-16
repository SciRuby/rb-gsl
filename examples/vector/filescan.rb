#!/usr/bin/env ruby
require("gsl")

a, b, c, d = GSL::Vector.filescan("c.dat")
p a, b, c, d

a, b, c, d = GSL::Vector::Int.filescan("c.dat")
p a, b, c, d

