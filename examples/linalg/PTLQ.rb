#!/usr/bin/env ruby
require("gsl")
include GSL

m = Matrix::alloc([0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                  [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85])
m.transpose!
lq, tau = m.PTLQ_decomp
p lq.class
p tau.class

b = Vector[1, 2, 3, 4]

lq, tau, p = m.PTLQ_decomp

p Linalg::PTLQ.solve_T(lq, tau, p, b)
p Linalg::PTLQ.solve_T(m, b)
p m.PTLQ_solve_T(b)
p lq.PTLQ_solve_T(tau, p, b)
p m.PTLQ_solve_T(b)

bb = b.clone
p Linalg::PTLQ.solve_T(lq, tau, p, bb)
bb = b.clone
p Linalg::PTLQ.solve_T(m, bb)
bb = b.clone
p m.PTLQ_solve_T(b)
bb = b.clone
p lq.PTLQ_solve_T(tau, p, bb)
bb = b.clone
p m.PTLQ_solve_T(bb)
bb = b.clone
p lq.class
p lq.solve_T(tau, p, bb)

q, l, tau, p = m.PTLQ_decomp2
p q.class
p l.class

p Linalg::PTLQ.LQsolve_T(q, l, p, b)

bb = b.clone
lq.svx_T(tau, p, bb)
p bb
bb = b.clone
m.PTLQ_svx_T(bb)
p bb
