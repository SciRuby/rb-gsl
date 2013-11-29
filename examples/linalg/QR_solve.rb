#!/usr/bin/env ruby
require("gsl")
include GSL

m = Matrix[[0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
           [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85]]
qr, tau = m.QR_decomp

b = Vector[1, 2, 3, 4]

qr, tau = m.QR_decomp

p Linalg::QR.solve(qr, tau, b)
p Linalg::QR.solve(m, tau, b)
p Linalg::QR.solve(m, b)
begin
  p Linalg::QR.solve(qr, b)    # Error
rescue
  puts("Error Linalg::QR.solve(qr, b)")
end

x = Vector.alloc(4)
p m.QR_solve(b)
p m.QR_solve(tau, b)
p qr.QR_solve(tau, b)
p qr.class
p qr.solve(tau, b)

m.QR_solve(b, x)
p x
m.QR_solve(tau, b, x)
p x
qr.QR_solve(tau, b, x)
p x
qr.solve(tau, b, x)
p x

bb = b.clone
qr.svx(tau, bb)
p bb
bb = b.clone
m.QR_svx(bb)
p bb

bb = b.clone
Linalg::QR.svx(m, bb)
p bb
bb = b.clone
Linalg::QR.svx(qr, tau, bb)
p bb

p m.QR_lssolve(b)
p m.QR_lssolve(tau, b)
p qr.QR_lssolve(tau, b)
p qr.lssolve(tau, b)

p Linalg::QR.lssolve(m, b)
p Linalg::QR.lssolve(qr, tau, b)

m3 = Matrix.alloc([5, 4, 1, 1], [4, 5, 1, 1], [1, 1, 4, 2], [1, 1, 2, 4])
qr, tau = m3.QR_decomp
p qr
p tau

q, r = qr.unpack(tau)
p q
p r
p q*r

q, r = Linalg::QR.unpack(qr, tau)
p q
p r
p q*r

qr, tau = m.QR_decomp
q, r = Linalg::QR.unpack(qr, tau)
p Linalg::QR.QRsolve(q, r, b)
exit
