#!/usr/bin/env ruby
require("gsl")

p GSL::Sf.hyperg_0F1(2, 3)
r = GSL::Sf.hyperg_0F1_e(2, 3)
p r.val

