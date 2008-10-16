#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector.alloc(4)
v.fread("a.dat")
p v

v2 = GSL::Vector.alloc(4)
v2.fscanf("b.dat")
p v2



