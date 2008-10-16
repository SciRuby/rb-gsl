#!/usr/bin/env ruby
require("gsl")

na = NArray[1.0, 2, 3, 4, 5, 6, 7 ,8 ,9]
p na
v = GSL::Vector.alloc(na)
p v

v[3] = 99.9
p v
p na

na2 = NArray[12.3, 45.6, 78.9, 1.23, 4.56, 7.89]
p na2
vref = na2.to_gv_view
p vref
vref[1] = 0.00
p na2

m = NMatrix[[1.0, 2],[3, 4]]
p m

gm = GSL::Matrix.alloc(m)
p gm

gm.set(1, 1, 99.9)
p gm
p m

gm2 = m.to_gm_view

gm2.set(1, 1, 99.9)
p gm2
p m

m = NMatrix[[0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
            [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85]]
mref = GSL::Matrix.to_gm_view(m)
mref[1][1] = 123
p m

mm = GSL::Matrix.to_gm(m)
mm[1][1] = 456
p m

