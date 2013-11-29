#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "roots/test.c"
require("gsl")
require("./gsl_test2.rb")
include GSL::Test
include Math

EPSREL = 10.0*GSL::DBL_EPSILON
EPSABS = 10.0*GSL::DBL_EPSILON
MAX_ITERATIONS = 150

def WITHIN_TOL(a, b, epsrel, epsabs)
  (a - b).abs < epsrel*GSL::MIN(a.abs, b.abs) + epsabs
end

def test_f(type, desc, f, lower, upper, correct)
  s = GSL::Root::FSolver.alloc(type)
  s.set(f, lower, upper)
  iter = 0
  begin
    iter += 1
    s.iterate
    r = s.root
    a = s.x_lower
    b = s.x_upper
    if a > b
      GSL::Test::test(GSL::FAILURE, "interval is invalid (#{a},#{b})")
    end
    if r < a or r > b
      GSL::Test::test(GSL::FAILURE, "r lies outside interval #{r} (#{a},#{b})")
    end
    status = GSL::Root::test_interval(a, b, EPSABS, EPSREL)
  end while status == GSL::CONTINUE and iter < MAX_ITERATIONS
  GSL::Test::test(status, "#{s.name}, #{desc} (#{s.root} obs vs #{correct} expected)")
  if iter == MAX_ITERATIONS
    GSL::Test::test(GSL::FAILURE, "exceeded maximum number of iterations")
  end
  if !WITHIN_TOL(r, correct, EPSREL, EPSABS)
    GSL::Test::test(GSL::FAILURE, "incorrect precision (#{r} obs vs #{correct} expected)")
  end
end

def test_fdf(type, desc, fdf, root, correct)
  s = GSL::Root::FdfSolver.alloc(type)
  s.set(fdf, root)
  iter = 0
  begin
    iter += 1
    prev = s.root
    s.iterate
    status = GSL::Root::test_delta(s.root, prev, EPSABS, EPSREL)
  end while status == GSL::CONTINUE and iter < MAX_ITERATIONS
  GSL::Test::test(status, "#{s.name} #{desc} (#{s.root} obs vs #{correct} expected)")
  if iter == MAX_ITERATIONS
    GSL::Test::test(GSL::FAILURE, "exceeded maximum number of iterations")
  end
  if !WITHIN_TOL(s.root, correct, EPSREL, EPSABS)
    GSL::Test::test(GSL::FAILURE, "incorrect precision (#{r} obs vs #{correct} expected)")
  end    
end

func1 = GSL::Function.alloc { |x| pow(x, 20.0) - 1 }
func1_fdf = Proc.new { |x| 20.0*pow(x, 19.0) }

fdf1 = GSL::Function_fdf.alloc(func1.f, func1_fdf)

fsolvers = ["bisection", "brent", "falsepos"]
fsolvers.each do |type|
  test_f(type, "x^20 - 1 [0.1, 2]", func1, 0.1, 2.0, 1.0)
end

fdfsolvers = ["newton", "secant", "steffenson"]
fdfsolvers = ["newton", "steffenson"]
fdfsolvers.each do |type|
  test_fdf(type, "x^{20} - 1 {0.9}", fdf1, 0.9, 1.0)
end
