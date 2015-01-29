#!/usr/bin/env ruby
require("gsl")
include GSL
include Linalg

m = Matrix::alloc([0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                  [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85])
m.print

lu, perm = m.LU_decomp
m.print
lu.print

b = Vector[1, 2, 3, 4]
x = LU.solve(lu, perm, b)
x = LU.solve(lu, perm, b)
p x
p m

p m.LU_solve(b)
p m.LU_solve(perm, b)
p m.LU_solve(perm, b, x)
p lu.LU_solve(perm, b, x)

p x
LU.refine(m, lu, perm, b, x)
p x

#p b
#m.LU_svx(b)
#p b
#LU.svx(lu, perm, b)
#p b
#p m.LU_solve(perm, b)
#p m.LU_solve(perm, b, x)
#p lu.LU_solve(perm, b, x)

#exit

LU.svx(lu, perm, b)
p b

p m*x.col

b = [1, 2, 3, 4].to_gv
p LU.solve(m, b)

puts("")

m2 = Matrix.alloc([1, 2, 3, 6, 5, 4, 7, 8, 1], 3, 3)
m2.print
lu2, p2, sign = m2.LU_decomp

m3 = Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 0], 3, 3)
m3.print
lu, perm = m3.LU_decomp
LU.invert(lu, perm).print
p m3.class
p lu.class
m4 = m3.invert
LU.invert(m3).print
p m4.invert
p m3.invert.invert.invert

p m2.det
p lu2.det(p2, sign)
p m2.det(p2)

p LU.det(m2)

p LU.lndet(m2)
p LU.lndet(lu2)
p m2.lndet
p lu2.lndet

p LU.sgndet(m2)
p LU.sgndet(lu2, sign)
p m2.sgndet
p lu2.sgndet(sign)

p m
p m/b.col

__END__
