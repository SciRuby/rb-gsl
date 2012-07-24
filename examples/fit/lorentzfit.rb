#!/usr/bin/env ruby
require("gsl")

N = 100

y0 = 1.0
A = 2.0
x0 = 3.0
B = 4.0

r = GSL::Rng.alloc()
x = GSL::Vector.linspace(-10, 20, N)
y =  y0 + A/(GSL::pow_2(x-x0)+B) + 0.03*GSL::Ran::gaussian(r, 1, N)

coef, err, chi2, dof = GSL::MultiFit::FdfSolver.fit(x, y, "lorentzian")
y0 = coef[0]
amp = coef[1]
x0 = coef[2]
b = coef[3]
p coef
p err
GSL::graph(x, y, y0+amp/(GSL::pow_2(x-x0)+b))
