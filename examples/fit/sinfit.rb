#!/usr/bin/env ruby
require("gsl")

N = 100

y0 = 1
A = 1.5
fc = 0.7
phi = 0.2

r = GSL::Rng.alloc()
x = GSL::Vector.linspace(0, 12, N)
y =  y0 + A*GSL::Sf::sin(fc*x+phi) + 0.1*r.gaussian(1, N)

coef, err, chi2, dof = GSL::MultiFit::FdfSolver.fit(x, y, "sin")
y0 = coef[0]
amp = coef[1]
fc = coef[2]
phi = coef[3]
p coef
p err
GSL::graph(x, y, y0+amp*GSL::Sf::sin(fc*x+phi))
