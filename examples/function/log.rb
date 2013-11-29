#!/usr/bin/env ruby
require("gsl")

f_log = GSL::Function.alloc { |x| Math::log(x) }
x = GSL::Vector.linspace(1, 100, 10)
x2 = GSL::Vector.logspace2(1, 100, 10)
GSL::graph([x, f_log.eval(x)], [x2, f_log.eval(x2)], "-T X -C -g 3 -l x -S 4 -L 'log(x)'")
