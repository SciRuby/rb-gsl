#!/usr/bin/env ruby
require("gsl")

# Create data
r = GSL::Rng.alloc("knuthran")
amp1 = 5.0
x01 = 1.0
sigma1 = 1.5

amp2 = 2.0
x02 = 5.0
sigma2 = 0.5

y0 = 2.0
N = 300
x = GSL::Vector.linspace(-4, 9, N)
y = y0 + amp1*GSL::Ran::gaussian_pdf(x - x01, sigma1) + amp2*GSL::Ran::gaussian_pdf(x - x02, sigma2) + 0.05*GSL::Ran::gaussian(r, 1.0, N)

coef, err, chi2, dof =  GSL::MultiFit::FdfSolver.fit(x, y, "gaussian_2peak", [2, 4, 0.9, 1, 1, 4, 1])

p coef
y01 = coef[0]

amp1 = coef[1]*Math::sqrt(2*Math::PI)*sigma1
x01 = coef[2]
sigma1 = Math::sqrt(coef[3])

amp2 = coef[4]*Math::sqrt(2*Math::PI)*sigma2
x02 = coef[5]
sigma2 = Math::sqrt(coef[6])

y2 = y01 + amp1*GSL::Ran::gaussian_pdf(x - x01, sigma1) + amp2*GSL::Ran::gaussian_pdf(x - x02, sigma2)

GSL::graph(x, y, y2, "-C -g 3")
