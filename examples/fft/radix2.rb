#!/usr/bin/env ruby
require("gsl") 

n = 128
data = GSL::FFT::Complex::PackedArray.alloc(2*n)

data.set_real(0, 1.0)
c = GSL::Complex.alloc(1, 0)
for i in 1..10 do
  data[i] = [1.0, 0.0]
  data[n-i] = [1.0, 0.0]
end

ffted = data.radix2_forward
ffted /= Math::sqrt(n)

GSL::graph(nil, data.re, ffted.re, "-T X -C -g 3 -L 'Radix-2' -x 0 #{data.size/2}")

