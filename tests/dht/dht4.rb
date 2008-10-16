#!/usr/bin/env ruby
# Test gsl/dht/test.c: test_dht_poly1()
# Expected results:
#  vout[0]=0.057274214
#  vout[5]=-0.000190850
#  vout[10]=0.000024342
#  vout[35]=-4.04e-07
#  vout[100]=1.0e-08

require("gsl")

N = 128
t = Dht.alloc(N, 1.0, 1.0)
vin = GSL::Vector.alloc(N)
for n in 0...N do
  x = t.x_sample(n)
  vin[n] = x*(1.0 - x*x)
end

vout = t.apply(vin)

printf("%e %e %e %e %e\n", vout[0], vout[5], vout[10], vout[35], vout[100])

