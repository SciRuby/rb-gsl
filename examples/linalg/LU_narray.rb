#!/usr/bin/env ruby
require("gsl")
include GSL
include Linalg

m = NArray[[0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
           [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85]]
b = NArray[1.0, 2, 3, 4]

lu, perm, signum = LU.decomp(m)
p m
p lu

x = LU.solve(lu, perm, b)
p x
p b

LU.svx(lu, perm, b)
p b

m2 = Matrix.alloc([1, 2, 3, 6, 5, 4, 7, 8, 1], 3, 3).to_nmatrix
lu, perm, signum = LU.decomp(m2)
p LU.det(lu, signum)

