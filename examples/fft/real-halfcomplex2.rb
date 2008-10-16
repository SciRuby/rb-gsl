#!/usr/bin/env ruby
require("gsl")

n = 100
data = GSL::Vector.alloc(n)
#data = NArray.float(n)

for i in (n/3)...(2*n/3) do
  data[i] = 1.0
end
org = data.clone

rtable = GSL::FFT::Real::Wavetable.alloc(n)
rwork = GSL::FFT::Real::Workspace.alloc(n)
  
data.real_transform!(rtable, rwork)

for i in 11...n do
  data[i] = 0.0
end
  
hctable = GSL::FFT::HalfComplex::Wavetable.alloc(n)
  
#data2 = ffted.halfcomplex_inverse(hctable, rwork)
#data2 = ffted.halfcomplex_inverse()
#data2 = ffted.ifft
#data2 = GSL::FFT::HalfComplex.inverse(ffted)
#data.halfcomplex_inverse!(hctable, rwork)
data.ifft!(hctable, rwork)

GSL::graph(nil, org, data, "-T X -C -g 3 -L 'Real-halfcomplex' -x 0 #{data.size}")
