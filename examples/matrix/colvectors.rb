#!/usr/bin/env ruby
require("gsl")

m = GSL::Matrix.alloc(1..9, 3, 3)
n = GSL::Matrix.alloc(11..19, 3, 3)

a = m.col(1)
b = n.col(2)

c = GSL::Matrix[a, b]
p m
p n
p c
