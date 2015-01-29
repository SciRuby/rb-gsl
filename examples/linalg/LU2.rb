#!/usr/bin/env ruby
require("gsl")
include GSL
include Linalg

m = Matrix::alloc([0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                  [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85])
m.print

lu, perm = m.LU_decomp

b = [1, 2, 3, 4]
x = Vector.alloc(4)

p m.class

p LU.solve(lu, perm, b)
LU.solve(lu, perm, b, x)
p x
p LU.solve(m, b)
LU.solve(m, b, x)
p x

LU.solve(m, perm, b, x)
p x

bv = b.to_gv
LU.svx(m, perm, bv)
p bv

m.print
