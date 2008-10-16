#!/usr/bin/env ruby
require("gsl")

a = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)
p a

p a.row(1)
p a.col(0)
p a.diagonal
p a.subdiagonal(1)
p a.subdiagonal(0)
p a.superdiagonal(1)


