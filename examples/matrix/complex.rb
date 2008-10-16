#!/usr/bin/env ruby
require("gsl")

m = GSL::Matrix::Complex.alloc(3, 3)

m.set(1, 2, 3, 5.6)

m.print

a = m.get(1, 2)
p a
p a.class

m2 = m.submatrix(1, 1, 2, 2)
p m2
m2.print

row = m.row(1)
p row
col = m.col(2)
p col

m.each_row do |v|
  p v
end

m.each_col do |v|
  p v
end
