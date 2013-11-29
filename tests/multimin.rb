#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "multimin/test.c"
require("gsl")
require("./gsl_test2.rb")
include GSL::Test
include Math

def test_fdf(desc, f, initpt, type)
  x = eval("#{initpt}")
  step_size = 0.1*Blas.dnrm2(x)
  s = GSL::MultiMin::FdfMinimizer.alloc(type, f.n)
  s.set(f, x, step_size, 0.1)
  iter = 0
  begin
    iter += 1
    status = s.iterate
    status = GSL::MultiMin.test_gradient(s.gradient, 1e-3)
  end while iter < 5000 and status == GSL::CONTINUE
  status |= s.f.abs > 1e-5 ? 1 : 0
  GSL::Test::test(status, "#{s.name}, on #{desc}: #{iter} iterations, f(x)=#{s.f}")
end

def test_f(desc, f, initpt)
  x = eval("#{initpt}")
  step_size = GSL::Vector.alloc(f.n)
  for i in 0...f.n
    step_size[i] = 1
  end
  s = GSL::MultiMin::FMinimizer.alloc("nmsimplex", f.n)
  s.set(f, x, step_size)
  iter = 0
  begin
    status = s.iterate
    status = GSL::MultiMin.test_size(s.size, 1e-3)
  end while iter < 5000 and status == GSL::CONTINUE

  status |= s.fval.abs > 1e-5 ? 1 : 0
  GSL::Test::test(status, "#{s.name}, on #{desc}: #{iter} iterations, f(x)=#{s.fval}")

  s = GSL::MultiMin::FMinimizer.alloc("nmsimplex2rand", f.n)
  s.set(f, x, step_size)
  iter = 0
  begin
    status = s.iterate
    status = GSL::MultiMin.test_size(s.size, 1e-3)
  end while iter < 5000 and status == GSL::CONTINUE

  status |= s.fval.abs > 1e-5 ? 1 : 0
  GSL::Test::test(status, "#{s.name}, on #{desc}: #{iter} iterations, f(x)=#{s.fval}")
end

def roth_initpt
  return GSL::Vector.alloc(4.5, 3.5)
end

def wood_initpt
  return GSL::Vector.alloc(-3.0, -1.0, -3.0, -1.0)
end

def rosenbrock_initpt
  return GSL::Vector.alloc(-1.2, 1.0)
end

Roth_f = Proc.new { |x|
  u = x[0]
  v = x[1]
  a = -13.0 + u + ((5.0 - v)*v - 2.0)*v;
  b = -29.0 + u + ((v + 1.0)*v - 14.0)*v;
  a * a + b * b;
}

Roth_df = Proc.new { |x, df|
  u = x[0]
  v = x[1]
  a = -13.0 + u + ((5.0 - v)*v - 2.0)*v
  b = -29.0 + u + ((v + 1.0)*v - 14.0)*v
  c = -2 + v * (10 - 3 * v)
  d = -14 + v * (2 + 3 * v)
  df[0] = 2 * a + 2 * b
  df[1] = 2 * a * c + 2 * b * d
}

Wood_f = Proc.new { |x|
  u1 = x[0]; u2 = x[1]; u3 = x[2]; u4 = x[3]
  t1 = u1*u1 - u2
  t2 = u3*u3 - u4
  100 * t1 * t1 + (1 - u1) * (1 - u1) + 90 * t2 * t2 + (1 - u3) * (1 - u3) + 10.1 * ( (1 - u2) * (1 - u2) + (1 - u4) * (1 - u4) ) + 19.8 * (1 - u2) * (1 - u4)
}

Wood_df = Proc.new { |x, df|
  u1 = x[0]; u2 = x[1]; u3 = x[2]; u4 = x[3]
  t1 = u1*u1 - u2
  t2 = u3*u3 - u4
  df[0] = 400 * u1 * t1 - 2 * (1 - u1)
  df[1] = -200 * t1 - 20.2 * (1 - u2) - 19.8 * (1 - u4)
  df[2] =  360 * u3 * t2 - 2 * (1 - u3)
  df[3] = -180 * t2 - 20.2 * (1 - u4) - 19.8 * (1 - u2)
}

Rosenbrock_f = Proc.new { |x|
  u = x[0]; v = x[1]
  a = u - 1
  b = u*u - v
  a*a + 10.0*b*b
}

Rosenbrock_df = Proc.new { |x, df|
  u = x[0]; v = x[1]
  _ = u - 1
  b = u*u - v
  df[0] = 2 * (u - 1) + 40 * u * b
  df[1] = -20 * b
}

if GSL::GSL_VERSION >= "1.8.90"
  fdfminimizers = ["steepest_descent", "conjugate_pr", "conjugate_fr",
                   "vector_bfgs", "vector_bfgs2"]
else
  fdfminimizers = ["steepest_descent", "conjugate_pr", "conjugate_fr",
                   "vector_bfgs"]
end

Rothdf = GSL::MultiMin::Function_fdf.alloc(Roth_f, Roth_df, 2)
Wooddf = GSL::MultiMin::Function_fdf.alloc(Wood_f, Wood_df, 4)
Rosenbrockdf = GSL::MultiMin::Function_fdf.alloc(Rosenbrock_f, Rosenbrock_df, 2)
fdfminimizers.each do |t|
  test_fdf("Roth", Rothdf, "roth_initpt", t)
  test_fdf("Wood", Wooddf, "wood_initpt", t)
  test_fdf("Rosenbrock", Rosenbrockdf, "rosenbrock_initpt", t)
end

rothf = GSL::MultiMin::Function.alloc(Roth_f, 2)
test_f("Roth", rothf, "roth_initpt")

woodf = GSL::MultiMin::Function.alloc(Wood_f, 4)
test_f("Wood", woodf, "wood_initpt")

rosenbrockf = GSL::MultiMin::Function.alloc(Rosenbrock_f, 2)
test_f("Rosenbrock", rosenbrockf, "rosenbrock_initpt")
