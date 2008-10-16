#!/usr/bin/env ruby
require("gsl") 

m = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)

p m

m.set([6,5, 6], [4, 5,7], [8, 5, 21])
m.print

m.set_col(1, [12, 3, 55].to_gv)
#m.print

m2 = m.transpose
m.print
m2.print

m.swap_rows(1, 2)
m.print 

v = m.col(0)
p v.to_a

m.diagonal.to_a

m = GSL::Matrix.alloc([1, 2, 3], [6, 5, 4], [7, 8, 1])
p m.get(1, 2)
m.print

m3 = m.LU_decomp
p m3

m3 = GSL::Matrix.alloc(5, 5)

for i in 0...5 do
  for j in 0...5 do
    val = 0.5*(i+0.4)*(j+1.2)
    m3.set(i, j, val)
  end
end

m3.print



