#!/usr/bin/env ruby
require("gsl")

n = 630
data = GSL::Vector::Complex.alloc(n)

table = GSL::FFT::ComplexWavetable.alloc(n)
space = GSL::FFT::ComplexWorkspace.alloc(n)

data[0] = 1.0
for i in 1..10 do
  data[i] = GSL::Complex[1.0, 0.0]
  data[n-i] = GSL::Complex[1.0, 0.0]
end
org = data.clone

# Select whichever you like
#data.forward!(table, space)
#data.forward!(table)
#data.forward!(space)
#data.forward!()
#data.transform!(table, space, GSL::FFT::Forward)
#data.transform!(GSL::FFT::Forward)
data.forward!(table, space)
data /= Math::sqrt(n)
GSL::graph(nil, org.re, data.re, "-C -g 3")
