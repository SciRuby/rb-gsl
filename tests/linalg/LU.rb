#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test
include Linalg

m = GSL::Matrix::alloc([0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85])

A = m.clone

lu_exp = GSL::Matrix.alloc([0.51, 0.13, 0.19, 0.85], 
                    [0.352941176470588, 0.554117647058823, 0.502941176470588, 0.66],
                    [0.803921568627451, 0.244515215852796, 0.71427813163482, -0.264713375796178], [0.274509803921569, 0.476999292285916, 0.949126848480345, 0.363093705877982])

x_exp = GSL::Vector[-4.05205022957397, -12.6056113959069, 1.66091162670884, 8.69376692879523]

lu, perm, signum = m.LU_decomp
test2(m == A, "#{A.class}#LU_decomp: matrix not destroyed")
test2(lu == lu_exp, "#{A.class}#LU_decomp")
b = GSL::Vector[1, 2, 3, 4]
x = LU.solve(lu, perm, b)
test2(x == x_exp, "#{A.class}.LU_solve")

x = lu.solve(perm, b)
test2(x == x_exp, "#{lu.class}#solve")

perm, signum = m.LU_decomp!
test2(m == lu_exp, "#{A.class}#LU_decomp!")

m = A.clone

x = LU.solve(m, perm, b)
test2(x == x_exp, "#{A.class}.LU_solve")
x = m.LU_solve(perm, b)
test2(x == x_exp, "#{A.class}#LU_solve")
test2(m == A, "#{A.class}#LU_solve: matrix not destroyed")

h = GSL::Matrix.hilbert(5)
invh = GSL::Matrix.invhilbert(5)
lu, perm, _sign = h.LU_decomp
a = Linalg::LU::invert(lu, perm)
test2(a.equal?(invh, 1e-6), "#{h.class}#LU_invert, Hilbert matrix of order 5")
a =lu.invert(perm)
test2(a.equal?(invh, 1e-6), "#{h.class}#LU_invert, Hilbert matrix of order 5")
a = h.inv
test2(a.equal?(invh, 1e-6), "#{h.class}#LU_invert, Hilbert matrix of order 5")
