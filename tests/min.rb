#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "min/test.c"
require("gsl")
require("./gsl_test2.rb")
include GSL::Test
include Math

GSL::IEEE::env_setup()

EPSABS = 0.001
EPSREL = 0.001
MAX_ITERATIONS = 100

def within_tol(a, b, epsrel, epsabs)
  (a - b).abs < (epsrel*GSL::MIN(a.abs, b.abs) + epsabs)
end

def test_f(type, desc, f, lower, mid, upper, min)
  x_lower = lower
  x_upper = upper
  s = GSL::Min::FMinimizer.alloc(type)
  s.set(f, mid, x_lower, x_upper)
  iterations = 0
  begin
    iterations += 1
    status = s.iterate
    m = s.x_minimum
    a = s.x_lower
    b = s.x_upper
    if a > b
      desc2 = sprintf("interval is invalid (%g,%g)", a, b)
      GSL::Test::test(GSL::FAILURE, desc2)
    end
    if m < a or m > b
      desc2 = sprintf("m lies outside interval %g (%g,%g)", m, a, b)
      GSL::Test::test(GSL::FAILURE, desc2)
    end
    if status == 1; break; end
    status = GSL::Min.test_interval(a, b, EPSABS, EPSREL)
  end while status == GSL::CONTINUE and iterations < MAX_ITERATIONS
  desc2 = sprintf("%s, %s (%g obs vs %g expected) ", s.name, desc, s.x_minimum, min)
  GSL::Test::test(status, desc2)
  if !within_tol(m, min, EPSREL, EPSABS)
    desc2 = sprintf("incorrect precision (%g obs vs %g expected)", m, min)
   GSL::Test::test(GSL::FAILURE, desc2)
  end
end

def test_f_e(type, desc, f, lower, mid, upper, min)
  x_lower = lower
  x_upper = upper
  s = GSL::Min::FMinimizer.alloc(type)
  status = nil
  begin
    status = s.set(f, mid, x_lower, x_upper)
  rescue
  end
  if status != GSL::SUCCESS
    desc2 = sprintf("%s, %s", s.name, desc)
    GSL::Test::test(status == GSL::SUCCESS ? 1 : 0, desc2)
  end
  iterations = 0
  begin
    iterations += 1
    status = s.iterate
    _ = s.x_minimum
    a = s.x_lower
    b = s.x_upper
    status = GSL::Min.test_interval(a, b, EPSABS, EPSREL)
  rescue
  end while status == GSL::CONTINUE and iterations < MAX_ITERATIONS
#  desc2 = sprintf("%s, %s", s.name, desc)
#  GSL::Test::test(status == 0 ? 0 : 1, desc2)
end

F_cos = GSL::Function.alloc { |x| cos(x) }
Func1 = GSL::Function.alloc { |x| pow(x, 4.0) - 1 }
Func2 = GSL::Function.alloc { |x| sqrt(x.abs) }
Func3 = GSL::Function.alloc { |x| 
  if x < 1.0
    1
  else
    -exp(-x)
  end
}
Func4 = GSL::Function.alloc { |x| x - 30.0/(1.0 + 1e5*pow(x - 0.8, 2.0)) }

types = ["goldensection", "brent", "quad_golden"]
types.each do |t|
  test_f(t, "cos(x) [0 (3) 6]", F_cos, 0.0, 3.0, 6.0, M_PI)
  test_f(t, "x^4 - 1 [-3 (-1) 17]", Func1, -3.0, -1.0, 17.0, 0.0);
  test_f(t, "sqrt(|x|) [-2 (-1) 1.5]", Func2, -2.0, -1.0, 1.5, 0.0);
  test_f(t, "func3(x) [-2 (3) 4]", Func3, -2.0, 3.0, 4.0, 1.0);
  test_f(t, "func4(x) [0 (0.782) 1]", Func4, 0, 0.782, 1.0, 0.8);
  
  test_f_e(t, "invalid range check [4, 0]", F_cos, 4.0, 3.0, 0.0, M_PI);
  test_f_e(t, "invalid range check [1, 1]", F_cos, 1.0, 1.0, 1.0, M_PI);
  test_f_e(t, "invalid range check [-1, 1]", F_cos, -1.0, 0.0, 1.0, M_PI)
end
