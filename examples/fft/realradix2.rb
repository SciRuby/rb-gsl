#!/usr/bin/env ruby
require("gsl")

n = 128
data = GSL::Vector.alloc(n)
#data = NArray.float(n)

data[0] = 1.0
for i in 1..10 do
  data[i] = 1.0
  data[n-i] = 1.0
end

#ffted = data.radix2_transform()
#ffted = data.real_radix2_transform()
# etc...
ffted = data.radix2_transform()
ffted /= Math::sqrt(n)
GSL::graph(nil, data, ffted, "-T X -C -g 3 -L 'Real Radix-2' -x 0 #{data.size}")
