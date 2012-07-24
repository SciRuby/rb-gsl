#!/usr/bin/env ruby
require("gsl")

N = 100

y0 = 0.5
m = 1.0
x0 = 1
rr = 0.2

r = GSL::Rng.alloc()
x = GSL::Vector.linspace(0, 2, N)
y =  y0 + m/(1.0 + GSL::Sf::exp((x0-x)/rr)) + 0.02*r.gaussian(1, N)

coef, err, chi2, dof = GSL::MultiFit::FdfSolver.fit(x, y, "sigmoid")
y0 = coef[0]
m = coef[1]
x0 = coef[2]
rr = coef[3]
p coef
p err
GSL::graph(x, y, y0+m/(1+GSL::Sf::exp((x0-x)/rr)))

=begin
Result:
GSL::Vector
[ 4.954e-01 1.004e+00 9.916e-01 1.991e-01 ]
GSL::Vector
[ 5.653e-03 9.100e-03 5.033e-03 4.856e-03 ]

GNUPLOT result:
Final set of parameters            Asymptotic Standard Error
=======================            ==========================

y0              = 0.495446         +/- 0.005656     (1.142%)
m               = 1.00422          +/- 0.009105     (0.9067%)
x0              = 0.991588         +/- 0.005033     (0.5076%)
r               = 0.199059         +/- 0.00486      (2.442%)

=end
