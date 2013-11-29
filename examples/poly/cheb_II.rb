#!/usr/bin/env ruby
require("gsl")

for n in 0..10
  p GSL::Poly.cheb_II(n)
end
