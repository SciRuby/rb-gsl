#!/usr/bin/env ruby
require("gsl")

N = 100000
k = 5

x = GSL::Vector.alloc(N)

GSL::Rng.env_setup()
T = GSL::Rng::DEFAULT
r = GSL::Rng.alloc(T)

for i in 0...N do
  x[i] = r.uniform()
end

small = x.sort_smallest(k)

printf("%d smallest values from %d\n", k, N);

for i in 0...k do
  printf("%d: %.18f\n", i, small[i]);
end
