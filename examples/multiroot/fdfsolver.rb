#!/usr/bin/env ruby
require("gsl")

n = 2

procf = Proc.new { |x, params, f|
  a = params[0]
  b = params[1]
  x0 = x[0]
  x1 = x[1]
  f[0] = a*(1 - x0)
  f[1] = b*(x1 - x0*x0)
}

procdf = Proc.new { |x, params, jac|
  a = params[0]
  b = params[1]
  jac[0][0] = -a
  jac[0][1] = 0
  jac[1][0] = -2*b*x[0]
  jac[1][1] = b
}

params = [1.0, 10.0]
f = GSL::MultiRoot::Function_fdf.alloc(procf, procdf, n, params)

fdfsolver = GSL::MultiRoot::FdfSolver.alloc("gnewton", n)
#fdfsolver = GSL::MultiRoot::FdfSolver.alloc("newton", n)
#fdfsolver = GSL::MultiRoot::FdfSolver.alloc("hybridj", n)
#fdfsolver = GSL::MultiRoot::FdfSolver.alloc("hybridsj", n)
p fdfsolver.name

#x = GSL::Vector.alloc(-10.0, -5.0)
x = [-10.0, -5.0]

#p fdfsolver.x

fdfsolver.set(f, x)

iter = 0
begin
  iter += 1
  status = fdfsolver.iterate
  root = fdfsolver.root
  f = fdfsolver.f
  printf("iter = %3u x = % .3f % .3f f(x) = % .3e % .3e\n",
          iter, root[0], root[1], f[0], f[1])
  status = fdfsolver.test_residual(1e-7)
end while status == GSL::CONTINUE and iter < 1000

