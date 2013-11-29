#!/usr/bin/env ruby
require("gsl")

v = GSL::Vector.alloc([17.2, 18.1, 16.5, 18.3, 12.6])
p GSL::Stats.mean(v)
p GSL::Stats.variance(v)
p GSL::Stats.sd(v)
p GSL::Stats.max(v)
p GSL::Stats.min(v)
