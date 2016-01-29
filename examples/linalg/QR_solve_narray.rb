#!/usr/bin/env ruby
require("gsl")
require 'narray'
include GSL

m = NMatrix[[0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
           [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85]]
qr, tau = Linalg::QR.decomp(m)

p qr, tau

b = NVector[1.0, 2, 3, 4]

p Linalg::QR.solve(qr, tau, b)