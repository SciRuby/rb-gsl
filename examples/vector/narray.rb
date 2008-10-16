#!/usr/bin/env ruby
require("gsl")

m = NMatrix[[0, 1.2, 1],[1.5, 0, 2]]
p m

gm = m.to_gv
p gm

m2 = gm.to_na
p m2
p m2.class

v = GSL::Vector.alloc(1..4)
p v
p v.class

na = v.to_na
p na
p na.class

v2 = na.to_gv
p v2

v3 = GSL::Vector.alloc(na)
p v3

v4 = na.to_gv_view
v4[2] = 123
p na
