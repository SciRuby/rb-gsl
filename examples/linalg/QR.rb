#!/usr/bin/env ruby
require("gsl")
include GSL

m = Matrix[[0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
           [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85]]
p m
qr, tau = m.QR_decomp
p qr.class
p tau.class

q, r = qr.unpack(tau)
p q.class
p r.class
p q*r
p q*q.trans   # Q: orthonormal


