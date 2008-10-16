#!/usr/bin/env ruby
require("gsl")

a = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)
p a

p a.max
p a.min
p a.minmax
p a.max_index
p a.min_index
p a.minmax_index

p a.isnull
p a.isnull?

a.set_zero
p a.isnull
p a.isnull?


