#!/usr/bin/env ruby
require("gsl")

a = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)
p a

b = a.get_row(0)
c = a.get_col(2)

a.set_col(2, b)
p a

a.set_row(0, c)
p a

a.swap_rows(1, 2)
p a

a.swap_cols(1, 2)
p a

p a.transpose
p a

p a.transpose!
p a

a.transpose!
p a
