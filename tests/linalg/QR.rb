#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test
include Linalg

m = GSL::Matrix::alloc([0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85])

A = m.clone

x_exp = GSL::Vector[-4.05205022957397, -12.6056113959069, 1.66091162670884, 8.69376692879523]
B = GSL::Vector[1, 2, 3, 4]

qr, tau = m.QR_decomp
test2(m == A, "#{m.class}#QR_decomp: matrix not destroyed")
x = m.QR_solve(B)
test2(x == x_exp, "#{m.class}#QR_solve(b)")
x = QR::solve(m, B)
test2(x == x_exp, "GSL::Linalg::QR::solve(b)")
tau = m.QR_decomp!
test2(m != A, "#{m.class}#QR_decomp: matrix destroyed")
x = m.QR_solve(tau, B)
test2(x == x_exp, "#{m.class}#QR_solve(tau, b)")
x = qr.solve(tau, B)
test2(x == x_exp, "#{qr.class}#solve(tau, b)")

begin
  x = m.QR_solve(B)
rescue
  puts("PASS: #{m.class}#QR_solve; tau vector must be given if m is decomped")
end
begin
  x = m.solve(B)
rescue
  puts("PASS: #{m.class}#solve; tau vector must be given if m is decomped")
end
x = m.solve(tau, B)
test2(x == x_exp, "#{m.class}#solve(tau, b)")

m = A.clone
bb = B.clone
m.QR_svx(bb)
test2(bb == x_exp, "#{m.class}#QR_svx(b)")

tau = QR::decomp!(m)
bb = B.clone
m.QR_svx(tau, bb)
test2(bb == x_exp, "#{m.class}#QR_svx(tau, b)")
begin
  x = m.QR_svx(bb)
rescue
  puts("PASS: #{m.class}#QR_solve; tau vector must be given if m is decomped")
end

m = A.clone
qr, tau = m.QR_decomp
test2(m == A, "#{m.class}#QR_decomp: matrix not destroyed")
x, r = m.QR_lssolve(B)
test2(x == x_exp, "#{m.class}#QR_lssolve(b)")

r = m.QR_lssolve(B, x)
test2(x == x_exp, "#{qr.class}#QR_lssolve(b, x)")
m.QR_lssolve(B, x, r)
test2(x == x_exp, "#{qr.class}#QR_lssolve(b, x, r)")
x, r = qr.QR_lssolve(tau, B)
test2(x == x_exp, "#{qr.class}#QR_lssolve(tau, b)")
r = qr.QR_lssolve(tau, B, x)
test2(x == x_exp, "#{qr.class}#QR_lssolve(tau, b, x)")
qr.QR_lssolve(tau, B, x, r)
test2(x == x_exp, "#{qr.class}#QR_lssolve(tau, b, x, r)")
begin
  x, r = qr.QR_lssolve(bb)
rescue
  puts("PASS: #{qr.class}#QR_solve; tau vector must be given if m is decomped")
end

