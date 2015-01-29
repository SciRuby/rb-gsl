#!/usr/bin/env ruby
require("gsl")

FFF = GSL::Function.alloc { |t, params|
  a = params[0]; lambda = params[1]; b = params[2]
  a*Math::exp(-lambda*t) + b
}

procf = Proc.new { |x, t, y, sigma, f|
  a = x[0]; lambda = x[1]; b = x[2]
  FFF.set_params(x)
  n = t.size
  for i in 0...n do
    yi = FFF.eval(t[i])
    f[i] = (yi - y[i])/sigma[i]
  end
}

procdf = Proc.new { |x, t, y, sigma, jac|
  a = x[0]; lambda = x[1]
  n = t.size
  for i in 0...n do
    ti = t[i]
    si = sigma[i]
    ei = Math::exp(-lambda*ti)
    jac[i,0] = ei/si
    jac[i,1] = -ti*a*ei/si
    jac[i,2] = 1.0/si
  end
}

n = 20
np = 3

f = GSL::MultiFit::Function_fdf.alloc(procf, procdf, np)

r = GSL::Rng.alloc()
t = GSL::Vector.alloc(n)
y = GSL::Vector.alloc(n)
sigma = GSL::Vector.alloc(n)
File.open("expdata.dat", "w") do |fp|
  for i in 0...n do
    t[i] = i
    y[i] = 1.0 + 5*Math::exp(-0.1*t[i]) + r.gaussian(0.2)
    sigma[i] = 0.2
    fp.printf("%d %g %g\n", t[i], y[i], sigma[i])
#    fp.printf("%d %g\n", t[i], y[i])
  end
end

f.set_data(t, y, sigma)
x = GSL::Vector.alloc([1.0, 0.0, 0.0])

solver = GSL::MultiFit::FdfSolver.alloc(GSL::MultiFit::FdfSolver::LMSDER, n, np)
#solver = GSL::MultiFit::FdfSolver.alloc(GSL::MultiFit::FdfSolver::LMDER, n, np)
solver.set(f, x)

iter = 0
solver.print_state(iter)
begin
  iter += 1
  status = solver.iterate
  solver.print_state(iter)
  status = solver.test_delta(1e-4, 1e-4)
end while status == GSL::CONTINUE and iter < 500

covar = solver.covar(0.0)
position = solver.position

chi2 = GSL::pow_2(solver.f.dnrm2)
dof = n - np
printf("A      = %.5f +/- %.5f\n", position[0], Math::sqrt(chi2/dof*covar[0,0]))
printf("lambda = %.5f +/- %.5f\n", position[1], Math::sqrt(chi2/dof*covar[1,1]))
printf("b      = %.5f +/- %.5f\n", position[2], Math::sqrt(chi2/dof*covar[2,2]))

str = sprintf("%4.3f*exp(-%4.3f*x)+%4.3f", position[0], position[1], position[2])

FFF.set_params(position)
x = 0
File.open("fit.dat", "w") do |f|
  while x < 40
    f.printf("%e %e\n", x, FFF.eval(x))
    x += 0.1
  end
end

system("graph -T X -C -g 3 -L '#{str}' -X x -I e -m -1 -S 4 expdata.dat -I a -m 2 -S 0 fit.dat")
#File.delete("expdata.dat")
#File.delete("fit.dat")
