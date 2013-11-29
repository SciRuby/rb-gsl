#!/usr/bin/env ruby
require("gsl")
include GSL::MultiFit

N = 100

y0 = 1
m = 2.0
xhalf = 5
rr = 4

r = GSL::Rng.alloc()
x = GSL::Vector.linspace(0.001, 10, N)
y =  y0 + (m-y0)/(1.0 + GSL::pow(xhalf/x, rr)) + 0.02*r.gaussian(1, N)

coef, err, chi2, dof = GSL::MultiFit::FdfSolver.fit(x, y, "hill")
y0 = coef[0]
m = coef[1]
xhalf = coef[2]
rr = coef[3]
p coef
p err
GSL::graph(x, y, y0+(m-y0)/(1+GSL::pow(xhalf/x, rr)))

=begin
Result:
GSL::Vector
[ 9.959e-01 1.995e+00 4.936e+00 4.035e+00 ]
GSL::Vector
[ 4.676e-03 1.035e-02 3.779e-02 1.125e-01 ]

GNUPLOT result:
Final set of parameters            Asymptotic Standard Error
=======================            ==========================

y0              = 0.995858         +/- 0.004676     (0.4695%)
m               = 1.99465          +/- 0.01034      (0.5184%)
xhalf           = 4.93605          +/- 0.0378       (0.7658%)
r               = 4.0346           +/- 0.1124       (2.787%)
=end
