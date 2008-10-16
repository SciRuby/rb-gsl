#!/usr/bin/env ruby
# vim: set ts=8 sw=2 sts=2:
#require("gsl")
#require("../gsl_test2.rb")
require("./test_multifit.rb")
include GSL::Test
include Math

Maxiter = 10
N = 1000

# model: a*exp(-(x-x0)**2/2/sigma**2)
gauss_p = 3
gauss_f = Proc.new { |x, t, y, s, f|
  # x: parameters as a Vecor
  # t: observed points as a GSL::Vector
  # y: observed data as a GSL::Vector
  # s: errorbar
  # f: result
  a = x[0]
  x0 = x[1]
  sigma2 = x[2]**2
  y.size.times do |i|
    f.set(i, (a*Math::exp(-(t[i] - x0)**2/2/sigma2) - y[i])/s[i])
  end
  GSL::SUCCESS
}

gauss_df = Proc.new { |x, t, y, s, df|
  a = x[0]
  x0 = x[1]
  sigma = x[2]
  sigma2 = sigma**2
  y.size.times do |i|
    dx = t[i] - x0
    dx2 = dx**2
    f = a*Math::exp(-dx2/2/sigma2)
    df.set(i, 0, f/a/s[i])
    df.set(i, 1, f*dx/sigma2/s[i])
    df.set(i, 2, f*dx2/sigma2/sigma/s[i])
  end
  GSL::SUCCESS
}

# goal
xgoal = GSL::Vector.alloc([1, 0, 1])
parname = %w( a x0 si )

# data
t = GSL::Vector.alloc(N) # positions of data
tmin = -10.0
tmax = 10.0
t.size.times do |i|
  t[i] = tmin + (tmax - tmin)*i/(N-1)
end
stdev = xgoal[0]*0.1 
s = GSL::Vector.alloc(Array.new(t.size, stdev))  # error bar of each datum
r = Rng.alloc
e = GSL::Vector.alloc(t.size)
t.size.times do |i|
  e[i] = -r.gaussian(stdev)	# perturbation to data
end
y = GSL::Vector.alloc(t.size)
n = GSL::Vector.alloc(Array.new(t.size, 1.0))
gauss_f.call(xgoal, t, e, n, y)  # data: y = model - e

# fitting
x = GSL::Vector.alloc([0.5, 0.1, 2]) # initial guess
fdf = GSL::MultiFit::Function_fdf.alloc(gauss_f, gauss_df, gauss_p)
fdf.set_data(t, y, s)

solver = GSL::MultiFit::FdfSolver.alloc(GSL::MultiFit::FdfSolver::LMSDER, t.size, gauss_p)
solver.set(fdf, x)

solver.print_state(0)
Maxiter.times do |i|
  solver.iterate
  status = solver.test_delta(1e-6, 1e-6)
  solver.print_state(i+1)
  break unless GSL::CONTINUE == status
end

# results
covar = solver.covar(0.0)
xresult = solver.position
dof = t.size - gauss_p
chi2 = pow_2(solver.f.dnrm2)
xsigma = GSL::Vector.alloc(xresult.size)
xresult.size.times do |i|
  xsigma[i] = Math::sqrt(chi2/dof*covar[i][i]) * 2.0 
  # resulting parameters to differ two times than standard error
end
puts("a*exp(-(x-x0)**2/2/si**2), chi2/N:%.3g" % (chi2/t.size))
xresult.size.times do |i|
  test_rel(xresult[i], xgoal[i], xsigma[i], "%-2.2s" % parname[i])
  exit 1 if (xresult[i] - xgoal[i]).abs > xsigma[i]
end
