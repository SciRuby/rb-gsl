#!/usr/bin/env ruby
require("gsl")

a = GSL::Vector.alloc(1, 2, 3, 4, 5)

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
