#!/usr/bin/env ruby
require("gsl")
include GSL
include GSL::Linalg

m = GSL::Matrix[[0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85]]

b = Vector.alloc(1, 2, 3, 4)
p HH.solve(m, b)

p m.HH_solve(b)

m.HH_svx(b)
p b
