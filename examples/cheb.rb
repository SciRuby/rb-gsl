#!/usr/bin/env ruby
require("gsl")

f = GSL::Function.alloc { |x|
  if x < 0.5
    0.25
  else
    0.75
  end
}

n = 1000
order = 40
cs = GSL::Cheb.alloc(order)

x = GSL::Vector.linspace(0, 1, n)
ff = f.eval(x)
cs.init(f, 0, 1)
r10 = cs.eval_n(10, x)
r40 = cs.eval(x)
GSL::graph(x, ff, r10, r40)
