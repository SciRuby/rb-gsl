#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "multiroot/test.c"
require("gsl")
require("./gsl_test2.rb")
include GSL::Test
include Math

GC.disable

def test_fdf(desc, fdf, initpt, factor, type)
  n = fdf.n
  x = eval("#{initpt}")
  if factor != 1.0
    x = x.scale(factor)
  end
  s = GSL::MultiRoot::FdfSolver.alloc(type, n)
  s.set(fdf, x)
  iter = 0
  begin
    iter += 1
    status = s.iterate
    status = GSL::MultiRoot.test_residual(s.f, 0.0000001)
  end while status == GSL::CONTINUE and iter < 1000

  jac, _stat = GSL::MultiRoot.fdjacobian(fdf, s.x, s.f, GSL::SQRT_DBL_EPSILON)
  r = 0.0
  sum = 0.0
  for i in 0...n
    for j in 0...n
      u = jac[i,j]
      su = s.jac[i,j]
      r = (u - su).abs/(1e-6 + 1e-6 * u.abs)
      sum += r
      if (u - su).abs > (1e-6 + 1e-6 * u.abs)
        printf("broken jacobian %g\n", r)
      end
    end
    printf("avg r = %g\n", sum/(n*n))
  end

  residual = 0.0
  for i in 0...n
    residual += s.f[i].abs
  end
  GSL::Test::test(status, "#{type} on #{desc} (#{factor}), #{iter} iterations, residual = #{residual}")

end

def test_f(desc, fdf, initpt, factor, type)
  n = fdf.n
  x = eval("#{initpt}")
  x = x.scale(factor)
  function = GSL::MultiRoot::Function.alloc(fdf.f, fdf.n)
  s = GSL::MultiRoot::FSolver.alloc(type, n)
  s.set(function, x)
  iter = 0
  begin
    iter += 1
    status = s.iterate
    status = GSL::MultiRoot.test_residual(s.f,  0.0000001)
  end while status == GSL::CONTINUE and iter < 1000
  residual = 0.0
  for i in 0...n
    residual += s.f[i].abs
  end
  GSL::Test::test(status, "#{type} on #{desc} (#{factor}), #{iter} iterations, residual = #{residual}")
end

def roth_initpt
  return GSL::Vector.alloc(4.5, 3.5)
end

def wood_initpt
  return GSL::Vector.alloc(-3.0, -1.0, -3.0, -1.0)
end

def rosenbrock_initpt
  GSL::Vector.alloc(-1.2, 1.0)
end

roth_f = Proc.new { |x, f|
  u = x[0]
  v = x[1]
  f[0] = -13.0 + u + ((5.0 - v)*v - 2.0)*v;
  f[1] = -29.0 + u + ((v + 1.0)*v - 14.0)*v;
}

roth_df = Proc.new { |x, df|
  x1 = x[1]
  df.set(0, 0, 1.0)
  df.set(0, 1, -3 * x1 * x1 + 10 * x1 - 2)
  df.set(1, 0, 1.0)
  df.set(1, 1, 3 * x1 * x1 + 2 * x1 - 14)
}

rosenbrock_f = Proc.new { |x, f|
  x0 = x[0]; x1 = x[1]
  y0 = 1.0 - x0
  y1 = 10*(x1 - x0*x0)
  f[0] = y0
  f[1] = y1
  GSL::SUCCESS
}
rosenbrock_df = Proc.new { |x, df|
  x0 = x[0]
  df00 = -1.0
  df01 = 0.0
  df10 = -20*x0
  df11 = 10
  df.set(0, 0, df00)
  df.set(0, 1, df01)
  df.set(1, 0, df10)
  df.set(1, 1, df11)
  GSL::SUCCESS
}

roth = GSL::MultiRoot::Function_fdf.alloc(roth_f, roth_df, 2)
rosenbrock = GSL::MultiRoot::Function_fdf.alloc(rosenbrock_f, rosenbrock_df, 2)

fsolvers = ["dnewton", "broyden", "hybrid", "hybrids"]
fsolvers.each do |type|
  test_f("Roth", roth, "roth_initpt", 1.0, type)
  test_f("Rosenbrock", rosenbrock, "rosenbrock_initpt", 1.0, type)
end

exit
f = 1.0
fdfsolvers = ["newton", "gnewton", "hybridj", "hybridsj"]
fdfsolvers.each do |type|
  test_fdf("Roth", roth, "roth_initpt", f, type)
end
