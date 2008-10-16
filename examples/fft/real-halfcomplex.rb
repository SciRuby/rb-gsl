#!/usr/bin/env ruby
require("gsl")

n = 100
data = GSL::Vector.alloc(n)
#data = NArray.float(n)

for i in (n/3)...(2*n/3) do
  data[i] = 1.0
end

rtable = GSL::FFT::Real::Wavetable.alloc(n)
rwork = GSL::FFT::Real::Workspace.alloc(n)
  
#ffted = data.real_transform(1, n, rtable, rwork)
#ffted = data.real_transform(n, rtable, rwork)
#ffted = data.real_transform(rtable, rwork)
#ffted = data.real_transform(n, rtable)
#ffted = data.real_transform(n, rwork)
#ffted = data.real_transform(n)
#ffted = data.real_transform(rtable)
#ffted = data.real_transform(rwork)
#ffted = data.real_transform()
#ffted = data.fft
ffted = GSL::FFT::Real.transform(data)

for i in 11...n do
  ffted[i] = 0.0
end
  
hctable = GSL::FFT::HalfComplex::Wavetable.alloc(n)
  
#data2 = ffted.halfcomplex_inverse(hctable, rwork)
#data2 = ffted.halfcomplex_inverse()
#data2 = ffted.ifft
data2 = GSL::FFT::HalfComplex.inverse(ffted)

GSL::graph(nil, data, data2, "-T X -C -g 3 -L 'Real-halfcomplex' -x 0 #{data.size}")
