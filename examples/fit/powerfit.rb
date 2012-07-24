#!/usr/bin/env ruby
require("gsl")

N = 100

y0 = 0.1
A = 2
rr = -1.5

r = GSL::Rng.alloc()
x = GSL::Vector.linspace(0.1, 1, N)
y =  y0 + A*GSL::pow(x, rr) + 0.1*r.gaussian(1, N)

coef, err, chi2, dof = GSL::MultiFit::FdfSolver.fit(x, y, "power")
y0 = coef[0]
amp = coef[1]
rr = coef[2]
p coef
p err

GSL::graph(x, y, y0+amp*GSL::pow(x, rr), "-T X -C -g 3 -l x -l y")
