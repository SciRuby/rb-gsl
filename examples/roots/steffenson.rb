#!/usr/bin/env ruby
require("gsl")
include Math

def f_to_solve(x)
  x*x - 5.0
end

def df_to_solve(x)
  2*x
end

f = Proc.new { |x|
  f_to_solve(x)
}

df = Proc.new { |x|
  df_to_solve(x)
}

expected = sqrt(5.0)
x = 5.0

function_fdf = GSL::Function_fdf.alloc(f, df)

solver = GSL::Root::FdfSolver.alloc(GSL::Root::FdfSolver::STEFFENSON)
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
