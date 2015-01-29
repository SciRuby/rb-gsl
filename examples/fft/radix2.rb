#!/usr/bin/env ruby
require("gsl")

n = 128
data = GSL::Vector::Complex.alloc(n)

data[0] = 1.0
c = GSL::Complex[1, 0]
for i in 1..10 do
  data[i] = c
  data[n-i] = c
end

ffted = data.radix2_forward
ffted /= Math::sqrt(n)

GSL::graph(nil, data.re, ffted.re, "-T X -C -g 3 -L 'Radix-2' -x 0 #{data.size/2}")

