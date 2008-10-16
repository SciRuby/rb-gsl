#!/usr/bin/env ruby
require("gsl")
p GSL::Matrix.diagonal(1..3)
p GSL::Matrix.diagonal([1, 2, 3])
p GSL::Matrix.diagonal([1, 2, 3].to_gv)
p GSL::Matrix.diagonal(1, 2, 3)
