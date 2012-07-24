#!/usr/bin/env ruby
require("gsl")
include Math

# Expected: x0 ~ 0.57983, x1 ~ 2.54621
#   (by Octave)
func = GSL::MultiRoot::Function.alloc(2) { |x, f|
  x0 = x[0]
  x1 = x[1]
  f[0] = -2.0*x0*x0 + 3.0*x0*x1 + 4.0*sin(x1) - 6.0
  f[1] = 3.0*x0*x0 - 2.0*x0*x1*x1 + 3.0*cos(x0) + 4.0
}

#fsolver = GSL::MultiRoot::FSolver.alloc("dallocton", 2)
fsolver = GSL::MultiRoot::FSolver.alloc("broyden", 2)
p fsolver.name

x = GSL::Vector.alloc([1.0, 2.0])
fsolver.set(func, x)

iter = 0
begin
  iter += 1
  status = fsolver.iterate
  root = fsolver.root
  f = fsolver.f
  status = fsolver.test_residual(1e-6)
end while status == GSL::CONTINUE and iter < 1000

x0 = fsolver.root[0]
x1 = fsolver.root[1]
printf("%f %f\n", x0, x1)
