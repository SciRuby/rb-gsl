#!/usr/bin/env ruby
require("gsl")

n = 630
data = GSL::FFT::Complex::PackedArray.alloc(2*n)

table = GSL::FFT::Complex::Wavetable.alloc(n)
space = GSL::FFT::Complex::Workspace.alloc(n)

data.set_real(0, 1.0)
for i in 1..10 do
  data[i] = [1.0, 0.0]
  data[n-i] = [1.0, 0.0]
end

# Select whichever you like
#ffted = data.forward(1, n, table, space)
#ffted = data.forward(n, table, space)
#ffted = data.forward(table, space)
#ffted = data.forward(n, table)
#ffted = data.forward(n, space)
#ffted = data.forward(n)
#ffted = data.forward(table)
#ffted = data.forward(space)
#ffted = data.forward()
#ffted = data.transform(table, space, GSL::FFT::Forward)
#ffted = data.transform(GSL::FFT::Forward)
#ffted = GSL::FFT::Complex.forward(data, table, space)
ffted = GSL::FFT::Complex.forward(data)
#ffted = GSL::FFT::Complex.transform(data, table, space, GSL::FFT::Forward)
#ffted = GSL::FFT::Complex.transform(data, GSL::FFT::Forward)
ffted /= Math::sqrt(n)
GSL::graph(nil, data.re, ffted.re, "-C -g 3")
