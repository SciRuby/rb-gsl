#!/usr/bin/env ruby
require("gsl")

m = GSL::Matrix.alloc(2, 2)
m.fread("a.dat")
p m


m2 = GSL::Matrix.alloc(2, 2)
m2.fscanf("b.dat")
p m2


