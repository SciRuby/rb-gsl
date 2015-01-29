#!/usr/bin/env ruby
require("gsl")
include Math

#solver = GSL::Root::FSolver.alloc(GSL::Root::FSolver::BISECTION)
#solver = GSL::Root::FSolver.alloc("falsepos")
#solver = GSL::Root::FSolver.alloc("brent")
solver = GSL::Root::FSolver.alloc(GSL::Root::FSolver::BRENT)
puts "using #{solver.name} method"

f = GSL::Function.alloc { |x, params|
  a = params[0]
  b = params[1]
  c = params[2]
  (a*x + b)*x + c
}

f.set_params(1, 0, -5)
expected = sqrt(5.0)

printf("%5s [%9s, %9s] %9s %10s %9s\n",
          "iter", "lower", "upper", "root",
          "err", "err(est)")

solver.set(f, 0.0, 5.0)
iter = 0
status = nil
while status != GSL::SUCCESS
  iter += 1
  status = solver.iterate
  r = solver.root
  xl = solver.x_lower
  xu = solver.x_upper
  status = solver.test_interval(0, 0.001)

  if status == GSL::SUCCESS
    printf("Converged:\n")
  end
  printf("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
         iter, xl, xu, r, r - expected, xu - xl)
end

p f.fsolve(0, 5)
