#!/usr/bin/env ruby
require("gsl")

h = GSL::Histogram.alloc(5)
p h
p h.size

h.set_ranges([1, 5, 23, 45, 67, 89])
p h.range

h.set_ranges_uniform(0.0, 100)
p h.range
p h.get_range(3)
p h.max
p h.bins
p h.find(55)

#File.open("smp.dat") do |f|
#  h.fscanf(f)
#end

h.fscanf("smp.dat")

p h.max_val
p h.max_bin

