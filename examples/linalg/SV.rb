#!/usr/bin/env ruby
require("gsl")

A = GSL::Matrix[[3, 5, 2], [6, 2, 1], [4, 7, 3]]
p A
u, v, s = A.SV_decomp
p u.class
p v.class
p s.class

# u and v are orthonormal
p u*u.trans
p v*v.trans

# Reconstruct the matrix A from u, v, and s.
p u*GSL::Matrix.diagonal(s)*v.trans
