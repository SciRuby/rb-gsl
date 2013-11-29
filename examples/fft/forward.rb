#!/usr/bin/env ruby
require("gsl")

n = 630
data = GSL::Vector::Complex::alloc(n)

table = GSL::FFT::ComplexWavetable.alloc(n)
space = GSL::FFT::ComplexWorkspace.alloc(n)

data[0] = 1.0
for i in 1..10 do
  data[i]   = GSL::Complex[1.0, 0.0]
  data[n-i] = GSL::Complex[1.0, 0.0]
end

# Select whichever you like
#ffted = data.forward(table, space)
#ffted = data.forward(table)
#ffted = data.forward(space)
#ffted = data.forward()
#ffted = data.transform(table, space, GSL::FFT::Forward)
#ffted = data.transform(GSL::FFT::Forward)
ffted = data.transform(GSL::FFT::Forward)
ffted /= Math::sqrt(n)
GSL::graph(nil, data.re, ffted.re, "-C -g 3")
