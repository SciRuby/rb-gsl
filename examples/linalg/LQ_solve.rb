#!/usr/bin/env ruby
require("gsl")
include GSL

m = GSL::Matrix::alloc([0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
        [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85])
m.transpose!
lq, tau = m.LQ_decomp

b = Vector[1, 2, 3, 4]

lq, tau = m.LQ_decomp

p Linalg::LQ.solve_T(lq, tau, b)
p Linalg::LQ.solve_T(m, tau, b)
p Linalg::LQ.solve_T(m, b)

begin
  p Linalg::LQ.solve_T(lq, b)    # Error
rescue
  puts("Error Linalg::LQ.solve_T(lq, b)")
end

x = Vector.alloc(4)
p m.LQ_solve_T(b)
p m.LQ_solve_T(tau, b)
p lq.LQ_solve_T(tau, b)
p lq.class
p lq.solve_T(tau, b)

m.LQ_solve_T(b, x)
p x
m.LQ_solve_T(tau, b, x)
p x
lq.LQ_solve_T(tau, b, x)
p x
lq.solve_T(tau, b, x)
p x

bb = b.clone
lq.svx_T(tau, bb)
p bb
bb = b.clone
m.LQ_svx_T(bb)
p bb

bb = b.clone
Linalg::LQ.svx_T(m, bb)
p bb
bb = b.clone
Linalg::LQ.svx_T(lq, tau, bb)
p bb

p m.LQ_lssolve_T(b)
p m.LQ_lssolve_T(tau, b)
p lq.LQ_lssolve_T(tau, b)
p lq.lssolve_T(tau, b)

p Linalg::LQ.lssolve_T(m, b)
p Linalg::LQ.lssolve_T(lq, tau, b)

m3 = Matrix.alloc([5, 4, 1, 1], [4, 5, 1, 1], [1, 1, 4, 2], [1, 1, 2, 4])
m3.transpose!
lq, tau = m3.LQ_decomp

q, l = lq.unpack(tau)

q, l = Linalg::LQ.unpack(lq, tau)

lq, tau = m.LQ_decomp
q, l = Linalg::LQ.unpack(lq, tau)
p Linalg::LQ.LQsolve_T(q, l, b)

