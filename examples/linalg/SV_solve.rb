#!/usr/bin/env ruby
require("gsl")
include GSL

m = Matrix[[-2, 1, 1], [1, -2, 1], [1, 1, -2], [-2, 1, 1]]
u, v, s = m.SV_decomp
p u.class
p v.class
p s.class

m2 = Matrix[[5, 6, 8, 4], [6, 2, 4, 2], [8, 2, 3, 5], [1, 7, 2, 3]]
u2, v2, s2 = m2.SV_decomp
p u2
p v2
p s2

m = GSL::Matrix[[0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85]]

u, v, s = Linalg::SV.decomp_jacobi(m)

b = [1, 2, 3, 4].to_gv

p Linalg::SV.solve(u, v, s, b)
p Linalg::SV.solve(u, v, s, [1, 2, 3, 4])

puts "OK"
p m.SV_solve(b)
p Linalg::SV.solve(m, b)

####
A = Matrix[1..4, 2, 2]
I = Matrix.identity(2)
Ainv = A.inv

u, v, s = A.SV_decomp
sm = s.to_m_diagonal
sinv = s.collect { |elm| 1.0/elm }.to_m_diagonal
a = u*sm*v.trans
ainv = v*sinv*u.trans
p a == A
p ainv == Ainv
p (u.trans*u) == I
p (v.trans*v) == I
p A*v == u*sm
p A.trans*u == v*sm



