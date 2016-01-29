#!/usr/bin/env ruby
require("gsl")
include GSL

m = NMatrix[[0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
            [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85]]

u, v, s = Linalg::SV.decomp(m)

puts "u ->"
p u

puts "v ->"
p v

puts "s ->"
p s

b = NArray[1.0, 2, 3, 4]

puts "solved ->"
p Linalg::SV.solve(u, v, s, b)
