#!/usr/bin/env ruby
require("gsl")
include Math

f = Proc.new { |x, params|
  a = params[0]
  b = params[1]
  c = params[2]
  (a*x + b)*x + c
}

df = Proc.new { |x, params|
  a = params[0]
  b = params[1]
  2.0*a*x + b
}

expected = sqrt(5.0)
x = 5.0

function_fdf = GSL::Function_fdf.alloc(f, df)
function_fdf.set_params([1.0, 0, -5])

solver = GSL::Root::FdfSolver.alloc(GSL::Root::FdfSolver::NEWTON)
puts "using #{solver.name} method"

solver.set(function_fdf, x)

printf("%-5s %10s %10s %10s\n",
          "iter", "root", "err", "err(est)")
iter = 0
status = nil
while status != GSL::SUCCESS
  iter += 1
  status = solver.iterate
  x0 = x
  x = solver.root
  status = GSL::Root::test_delta(x, x0, 0, 1e-3)

  if status == GSL::SUCCESS
    printf("Converged:\n")
  end
  printf("%5d %10.7f %+10.7f %10.7f\n",
         iter, x, x - expected, x - x0)

end
