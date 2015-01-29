#!/usr/bin/env ruby
require("gsl")
include Math

# Expected: x0 ~ 0.57983, x1 ~ 2.54621
func = GSL::MultiRoot::Function.alloc(2) { |x, f|
  x0 = x[0]
  x1 = x[1]
  f[0] = -2.0*x0*x0 + 3.0*x0*x1 + 4.0*sin(x1) - 6.0
  f[1] = 3.0*x0*x0 - 2.0*x0*x1*x1 + 3.0*cos(x0) + 4.0
}

p func.solve([1.0, 2.0].to_gv, 1000, 1e-7, "hybrids")
p func.solve([1.0, 2.0].to_gv, 1000, "broyden")
p func.solve([1.0, 2.0], "hybrid")
p func.solve([1.0, 2.0], 2000, "hybrid")
p func.solve([1.0, 2.0])

#fsolver = GSL::MultiRoot::FSolver.alloc("dallocton", 2)
fsolver = GSL::MultiRoot::FSolver.alloc("broyden", 2)

x = GSL::Vector.alloc([1.0, 2.0])
fsolver.set(func, x)

p ans = fsolver.solve()
#ans = MultiRoot::FSolver.solve(fsolver)
