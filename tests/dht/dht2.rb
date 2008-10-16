#!/usr/bin/env ruby
# Test gsl/dht/test.c: test_dht_simple()
# Expected results:
#  vout[0]=4.00
#  vout[5]=1.84
#  vout[10]=1.27
#  vout[35]=0.352
#  vout[100]=0.0237

require("gsl")

N = 128
t = Dht.alloc(N, 0.0, 100.0)
vin = GSL::Vector.alloc(N)
for n in 0...N do
  x = t.x_sample(n)
  vin[n] = 1.0/(1.0 + x*x)
end

vout = t.apply(vin)

printf("%e %e %e %e %e\n", vout[0], vout[5], vout[10], vout[35], vout[100])

