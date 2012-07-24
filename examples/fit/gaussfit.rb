#!/usr/bin/env ruby
require("gsl")

# Create data
r = GSL::Rng.alloc("knuthran")
sigma = 1.5
x0 = 1.0
amp = 2.0
y0 = 3.0
N = 100
x = GSL::Vector.linspace(-4, 6, N)
y = y0 + amp*GSL::Ran::gaussian_pdf(x - x0, sigma) + 0.02*GSL::Ran::gaussian(r, 1.0, N)

coef, err, chi2, dof =  GSL::MultiFit::FdfSolver.fit(x, y, "gaussian")
sigma2 = Math::sqrt(coef[3])
x02 = coef[2]
amp2 = coef[1]*Math::sqrt(2*Math::PI)*sigma
y02 = coef[0]
y2 = y02 + amp2*GSL::Ran::gaussian_pdf(x - x02, sigma2)

GSL::graph(x, y, y2, "-C -g 3 -x -4 6")

printf("Expect:\n")
printf("sigma = #{sigma}, x0 = #{x0}, amp = #{amp}, y0 = #{y0}\n")
printf("Result:\n")
printf("sigma = %5.4e +/- %5.4e\n", sigma2, err[3])
printf("   x0 = %5.4e +/- %5.4e\n", x02, err[2])
printf("  amp = %5.4e +/- %5.4e\n", amp2, err[1])
printf("   y0 = %5.4e +/- %5.4e\n", y02, err[0])
