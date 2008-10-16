#!/usr/bin/env ruby
# Test gsl/dht/test.c: test_dht_exact()
# Expected results:
#  vout[0]=0.375254649407520
#  vout[1]=(-0.133507872695560
#  vout[2]=0.044679925143840

require("gsl")
vin = GSL::Vector.alloc(1, 2, 3)
t = Dht.alloc(3, 1.0, 1.0)
vout = t.apply(vin)
p vout

vin2 = t.apply(vout)
vin2.scale!(13.323691936314223*13.323691936314223)
p vin2

