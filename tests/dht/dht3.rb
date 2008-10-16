#!/usr/bin/env ruby
# Test gsl/dht/test.c: test_dht_exp1()
# Expected results:
#  vout[0]=0.181
#  vout[5]=0.357
#  vout[10]=0.211
#  vout[35]=0.0289
#  vout[100]=0.00221

require("gsl")

N = 128
t = Dht.alloc(N, 1.0, 20.0)
vin = GSL::Vector.alloc(N)
for n in 0...N do
  x = t.x_sample(n)
  vin[n] = Math::exp(-x)
end

vout = t.apply(vin)

printf("%e %e %e %e %e\n", vout[0], vout[5], vout[10], vout[35], vout[100])

