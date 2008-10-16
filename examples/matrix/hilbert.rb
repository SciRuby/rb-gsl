#!/usr/bin/env ruby
require("gsl")

m = GSL::Matrix.hilbert(3)
p m

invm = m.inv
invm2 = GSL::Matrix.invhilbert(3)

p invm

p m*invm
p m*invm2
p invm == invm2
