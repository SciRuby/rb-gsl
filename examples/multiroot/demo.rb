#!/usr/bin/env ruby
require("gsl")

params = [1.0, 10.0]
func = GSL::MultiRoot::Function.alloc(2, params) { |x, params, f|
  a = params[0]
  b = params[1]
  x0 = x[0]
  x1 = x[1]
  f[0] = a*(1 - x0)
  f[1] = b*(x1 - x0*x0)
}

fsolver = GSL::MultiRoot::FSolver.alloc("hybrids", 2)
#fsolver = GSL::MultiRoot::FSolver.alloc("hybrid", 2)
#fsolver = GSL::MultiRoot::FSolver.alloc("dallocton", 2)
#fsolver = GSL::MultiRoot::FSolver.alloc("broyden", 2)

x = GSL::Vector.alloc(-10.0, -5.0)
fsolver.set(func, x)
#p fsolver.name
#p fsolver.x

iter = 0
IO.popen("graph -T X -C -g 3 -X x -Y y -S 4", "w") do |io|
  begin
    iter += 1
    status = fsolver.iterate
    root = fsolver.root
    f = fsolver.f
    printf("iter = %3u x = % .3f % .3f f(x) = % .3e % .3e\n",
           iter, root[0], root[1], f[0], f[1])
    io.printf("%e %e\n", root[0], root[1])
    status = fsolver.test_residual(1e-7)
  end while status == GSL::CONTINUE and iter < 1000
end
