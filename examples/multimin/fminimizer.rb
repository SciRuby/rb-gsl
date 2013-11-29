#!/usr/bin/env ruby
require("gsl")
include GSL
include GSL::MultiMin

np = 2

my_f = Proc.new { |v, params|
  x = v[0]; y = v[1]
  p0 = params[0]; p1 = params[1]
  10.0*(x - p0)*(x - p0) + 20.0*(y - p1)*(y - p1) + 30.0
}

my_func = MultiMin::Function.alloc(my_f, np)
my_func.set_params([1.0, 2.0])      # parameters

x = Vector.alloc([5, 7])
ss = Vector.alloc(np)
ss.set_all(1.0)

#minimizer = FMinimizer.alloc("nmsimplex", np)
minimizer = FMinimizer.alloc(MultiMin::FMinimizer::NMSIMPLEX, np)

minimizer.set(my_func, x, ss)

iter = 0
begin
  iter += 1
  status = minimizer.iterate()
  status = minimizer.test_size(1e-2)
  if status == GSL::SUCCESS
    puts("converged to minimum at")
  end
  x = minimizer.x
  printf("%5d ", iter);
  for i in 0...np do
    printf("%10.3e ", x[i])
  end
  printf("f() = %7.3f size = %.3f\n", minimizer.fval, minimizer.size);
end while status == GSL::CONTINUE and iter < 100

