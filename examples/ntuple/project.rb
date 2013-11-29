#!/usr/bin/env ruby
require("gsl")

sel_func = GSL::Ntuple::SelectFn.alloc { |data, params|
  x = data[0]; y = data[1]; z = data[2]
  scale = params

  e2 = x*x + y*y + z*z
  e2 > scale
}

lower = 1.5
sel_func.set_params(lower)

val_func = GSL::Ntuple::ValueFn.alloc { |data|
  x = data[0]; y = data[1]; z = data[2]
  x*x + y*y + z*z
}

v = GSL::Vector.alloc(3)
n = GSL::Ntuple.open("test.dat", v)

h = GSL::Histogram.alloc(100)
h.set_ranges_uniform(0, 10.0)

#Ntuple.project(h, n, val_func, sel_func)
n.project(h, val_func, sel_func)

h.graph("-C -X 'E2' -Y 'n' -L 'GSL::Ntuple, Select E2 > 1.5'")
File.delete("test.dat")

