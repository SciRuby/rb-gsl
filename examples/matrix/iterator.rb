#!/usr/bin/env ruby
require("gsl")

m = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8 ,9], 3, 3)

m.each_col do |v|
  p v
end

