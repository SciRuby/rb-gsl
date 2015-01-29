#!/usr/bin/env ruby
require("gsl")

N = 40

amp0 = 5.0
b0 = 0.1
y0 = 1.0

r = GSL::Rng.alloc()
x = GSL::Vector[0...N]
sigma = GSL::Vector[N]
sigma.set_all(0.1)
y = y0 + amp0*GSL::Sf::exp(-b0*x) + 0.1*r.gaussian(1, N)

coef, err, chi2, dof = GSL::MultiFit::FdfSolver.fit(x, sigma, y, "exponential")
y0 = coef[0]
amp = coef[1]
b = coef[2]
p coef
p err
GSL::graph(x, y, y0+amp*GSL::Sf::exp(-b*x))

# This will result in
# [ 1.019e+00 5.045e+00 1.040e-01 ]
# [ 3.385e-02 5.395e-02 2.826e-03 ]

# GNUPLOT results:
# y0              = 1.01925          +/- 0.03383      (3.319%)
# A               = 5.04536          +/- 0.05396      (1.069%)
# b               = 0.104049         +/- 0.002826     (2.716%)
