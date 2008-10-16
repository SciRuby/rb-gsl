#!/usr/bin/env ruby
require("gsl")
include Math

m = GSL::Matrix.alloc(10, 10)
for i in 0...10
  for j in 0...10
    m.set(i, j, sin(i) + cos(j))
  end
end

for j in 0...10
  col = m.col(j)
  d = GSL::Blas.dnrm2(col)
  printf("matrix column %d, norm = %g\n", j, d)
end
