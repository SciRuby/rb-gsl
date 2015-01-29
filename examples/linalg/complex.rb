#!/usr/bin/env ruby
require("gsl")
include GSL
include Linalg

m = Matrix::alloc([0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                  [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85])
m.inv.print
p m.det

#lu, perm = m.LU_decomp
#LU.invert(lu, perm).print
#p lu.det

zm = Matrix::Complex.alloc(4, 4)

zm.set_row(0, [0.18, 0], [0.60, 0], [0.57, 0], [0.96, 0])
zm.set_row(1, [0.41, 0], [0.24, 0], [0.99, 0], [0.58, 0])
zm.set_row(2, [0.14, 0], [0.30, 0], [0.97, 0], [0.66, 0])
zm.set_row(3, [0.51, 0], [0.13, 0], [0.19, 0], [0.85, 0])

p zm.inv
p zm.det

b = Vector::Complex.alloc(4)
b.set(1.0, 2, 3, 4)

#p zm.LU_solve(b)

#lu2, perm2, signum = zm.LU_decomp
lu2, perm2, signum = LU::decomp(zm)
p lu2
p lu2.det(signum)

p x = lu2.solve(perm2, b)

p zm*x

p lu2.LU_invert(perm2)
p lu2.LU_det(signum)

#p zm.LU_svx(b)

p Linalg::Complex::LU_invert(zm)
p Linalg::Complex::LU::invert(zm)
p Linalg::Complex::LU_invert(lu2, perm2)
p Linalg::Complex::LU::invert(lu2, perm2)
p zm.invert
p lu2.invert(perm2)

p Linalg::Complex::LU_det(zm)
p Linalg::Complex::LU::det(zm)
p Linalg::Complex::LU_det(lu2, signum)
p Linalg::Complex::LU::det(lu2, signum)
p zm.det
p lu2.det(signum)

