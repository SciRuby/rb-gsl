#!/usr/bin/env ruby
require("gsl")
include GSL


m = Matrix::alloc([0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                  [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85])
qr, tau = m.QRPT_decomp
p qr.class
p tau.class

b = Vector.alloc([1, 2, 3, 4])

qr, tau, p = m.QRPT_decomp

p Linalg::QRPT.solve(qr, tau, p, b)
p Linalg::QRPT.solve(m, b)
p m.QRPT_solve(b)
p qr.QRPT_solve(tau, p, b)
p m.QRPT_solve(b)

bb = b.clone
p Linalg::QRPT.solve(qr, tau, p, bb)
bb = b.clone
p Linalg::QRPT.solve(m, bb)
bb = b.clone
p m.QRPT_solve(b)
bb = b.clone
p qr.QRPT_solve(tau, p, bb)
bb = b.clone
p m.QRPT_solve(bb)
bb = b.clone
p qr.class
p qr.solve(tau, p, bb)

q, r, tau, p = m.QRPT_decomp2
p q.class
p r.class

p Linalg::QRPT.QRsolve(q, r, p, b)

bb = b.clone
qr.svx(tau, p, bb)
p bb
bb = b.clone
m.QRPT_svx(bb)
p bb
