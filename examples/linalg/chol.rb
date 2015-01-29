#!/usr/bin/env ruby
require("gsl")
include GSL
include Linalg

m = Matrix.alloc([4, 2], [2, 3])
c = Cholesky.decomp(m)
p c.class
p c

b = Vector[1, 2]
p Cholesky.solve(c, b)    # Expected [-0.125, 0.75]

begin
  m = Matrix::alloc([0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                  [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85])

  c = Cholesky.decomp(m)
  b = [1, 2, 3, 4].to_gv
  p Cholesky.solve(c, b)
rescue
  puts("Matrix must be positive definite.")
end

m = Matrix.pascal(6)
c = m.cholesky_decomp
a = c.lower
p (a*a.trans) == m

