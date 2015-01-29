#!/usr/bin/env ruby
require("gsl")
include Math

f = GSL::Function::alloc{ |x, params|
  a = params[0]
  b = params[1]
  c = params[2]
  (a*x + b)*x + c
}

p f.proc
p f.params
a = 1; b = 2; c = 3
f.set_params(a, b, c)
p f.params

p f.eval(2)
p f.call(4)


f.set { |x|
  x*x*x
}

p f.params

p f.eval(2)
p f[4]

f2 = GSL::Function.alloc { |x|
  sin(x) - log(x)*sqrt(x)
}

p f2.eval(2.5)
p f2.arity
