#!/usr/bin/env ruby
require("gsl")

N = 100

y0 = 1.0
A = 2.0
x0 = 3.0
w = 0.5

r = GSL::Rng.alloc()
x = GSL::Vector.linspace(0.01, 10, N)
sig = 1
y =  y0 + A*GSL::Sf::exp(-GSL::pow_2(GSL::Sf::log(x/x0)/w)) + 0.1*GSL::Ran::gaussian(r, sig, N)

coef, err, chi2, dof = GSL::MultiFit::FdfSolver.fit(x, y, "lognormal", [0, 3, 2, 1])
y0 = coef[0]
amp = coef[1]
x0 = coef[2]
w = coef[3]

p coef
p err

GSL::graph(x, y, y0+amp*GSL::Sf::exp(-GSL::pow_2(GSL::Sf::log(x/x0)/w)))

