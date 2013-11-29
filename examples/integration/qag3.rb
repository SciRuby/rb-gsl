#!/usr/bin/env ruby
require("gsl")
include Math

f = GSL::Function.alloc { |x| x*sin(1.0/x)*sqrt((1 - x).abs) }
w = GSL::Integration::Workspace.alloc(1000)

p f.qag(0, 3, w)
