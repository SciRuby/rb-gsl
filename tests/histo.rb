#!/usr/bin/env ruby
require("gsl")

h = GSL::Histogram.alloc(10, [0, 10])
p h

p h.get_range(2)

p h.range
p h.bin


