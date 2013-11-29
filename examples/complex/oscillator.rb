#!/usr/bin/env ruby
require("gsl")

m = 1.0
k = 1.0
b = 1.0
mode = GSL::Poly.complex_solve_quadratic(m, b, k)
p mode.class

f = GSL::Vector.linspace(0.01, 100, 1000)
p f.class
tf0 = 1.0/((m*mode[0]*f + b)*mode[0]*f + k)
tf1 = 1.0/((m*mode[1]*f + b)*mode[1]*f + k)
p tf0.class
GSL::graph(f, tf0.amp, tf1.amp, "-C -g 3 -l x -l y -L 'Transfer function: amplitude'")
GSL::graph(f, tf0.phase, tf1.phase, "-C -g 3 -l x -L 'Transfer function: phase'")

