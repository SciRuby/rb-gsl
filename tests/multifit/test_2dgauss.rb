#!/usr/bin/env ruby
# vim: set ts=8 sw=2 sts=2:
#require("gsl")
#require("../gsl_test2.rb")
require("./test_multifit.rb")
include GSL::Test
include Math

Maxiter = 10
N = 33

# point
class Point
  attr_accessor :x, :y
  def initialize(x = nil, y = nil)
    @x, @y = x, y
  end
end

# model: a*exp(-((x-x0)**2+(y-y0)**2)/2/sigma**2)
gauss_f = Proc.new { |x, t, y, s, f|
  # x: parameters as a Vecor
  # t: observed points as an Array
  # y: observed data as a GSL::Vector
  # s: errorbar
  # f: result
  a = x[0]
  x0 = x[1]
  y0 = x[2]
  sigma2 = x[3]**2
  y.size.times do |i|
    f.set(i, (a*Math::exp(-((t[i].x - x0)**2+(t[i].y - y0)**2)/2/sigma2) - y[i])/s[i])
  end
  GSL::SUCCESS
}

gauss_df = Proc.new { |x, t, y, s, df|
  a = x[0]
  x0 = x[1]
  y0 = x[2]
  sigma = x[3]
  sigma2 = sigma**2
  y.size.times do |i|
    dx = t[i].x - x0; dx2 = dx**2
    dy = t[i].y - y0; dy2 = dy**2
    f = a*Math::exp(-(dx2+dy2)/2/sigma2)
    df.set(i, 0, f/a/s[i])
    df.set(i, 1, f*dx/sigma2/s[i])
    df.set(i, 2, f*dy/sigma2/s[i])
    df.set(i, 3, f*(dx2+dy2)/sigma2/sigma/s[i])
  end
  GSL::SUCCESS
}

# goal
xgoal = GSL::Vector.alloc([1, 0, 0, 1])
parname = %w( a x0 y0 si )

# data
t = Array.new
tmin = -10.0
tmax = 10.0
N.times do |j|
  N.times do |i|
    t << Point.new(tmin + (tmax - tmin)*i/(N-1), tmin + (tmax - tmin)*j/(N-1))
  end
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
x = GSL::Vector.alloc([0.5, 0.1, -0.1, 2.0]) # initial guess
fdf = GSL::MultiFit::Function_fdf.alloc(gauss_f, gauss_df, x.size)
fdf.set_data(t, y, s)

solver = GSL::MultiFit::FdfSolver.alloc(GSL::MultiFit::FdfSolver::LMSDER, t.size, x.size)
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
dof = t.size - xresult.size
chi2 = pow_2(solver.f.dnrm2)
xsigma = GSL::Vector.alloc(xresult.size)
xresult.size.times do |i|
  xsigma[i] = Math::sqrt(chi2/dof*covar[i][i]) * 2.0
  # allow resulting parameters to differ two times than standard error
end
puts "a*exp(-((x-x0)**2+(y-y0)**2)/2/si**2), chi2/N:%.3g" % (chi2/t.size)
xresult.size.times do |i|
  test_rel(xresult[i], xgoal[i], xsigma[i], "%-2.2s" % parname[i])
  if (xresult[i] - xgoal[i]).abs > xsigma[i] then
    puts "#{'%-2.2s' % parname[i]} is #{xresult[i]} +- #{xsigma[i]}"
    exit 1
  end
end
