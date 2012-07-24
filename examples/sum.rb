#!/usr/bin/env ruby
require("gsl")
include Math

N = 20
Zeta2 = PI*PI/6.0

v = GSL::Vector.alloc(N)
sum = 0.0
for n in 0...N do
  np1 = n.to_f + 1.0
  v[n] = 1.0/(np1 * np1)
  sum += v[n]
end

#sum_accel, err, sum_plain, terms_used = v.accel_sum
#sum_accel, err, sum_plain, terms_used = v.sum_accel
sum_accel, err, sum_plain, terms_used = v.accel

printf("term-by-term sum = %.16f using %d terms\n", sum, N)
printf("term-by-term sum = %.16f using %d terms\n", sum_plain, terms_used)
printf("exact value      = %.16f\n", Zeta2)
printf("accelerated sum  = %.16f using %d terms\n", sum_accel, terms_used)
printf("Estimated error  = %.16f\n", err)
printf("Actual error     = %.16f\n", sum_accel - Zeta2)

p GSL::Sum::Levin_u.accel(v)
p GSL::Sum::Levin_utrunc.accel(v)

lu = GSL::Sum::Levin_u.alloc(N)
p lu.accel(v)

lutrunc = GSL::Sum::Levin_utrunc.alloc(N)
p lutrunc.accel(v)
