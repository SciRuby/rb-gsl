#!/usr/bin/env ruby
require("gsl")
fn1 = GSL::Function.alloc { |x| Math::cos(x) + 1.0 }
iter = 0;  max_iter = 500
m = 2.0             # initial guess
m_expected = Math::PI
a = 0.0; b = 6.0
gmf = GSL::Min::FMinimizer.alloc(GSL::Min::FMinimizer::BRENT)
gmf.set(fn1, m, a, b)
printf("Using %s method\n", gmf.name)
printf("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "min",
       "err", "err(est)")
printf("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
       iter, a, b, m, m - m_expected, b - a)
begin
  iter += 1
  status = gmf.iterate
  status = gmf.test_interval(0.001, 0.0)
  puts("Converged:") if status == GSL::SUCCESS
  a = gmf.x_lower;  b = gmf.x_upper
  m = gmf.x_minimum
  printf("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
         iter, a, b, m, m - m_expected, b - a);
end while status == GSL::CONTINUE and iter < max_iter

x = GSL::Vector.linspace(0, 6, 50)
mx = gmf.x_minimum
min = fn1.eval(mx)
GSL::graph([x, fn1.eval(x)], [GSL::Vector[mx], GSL::Vector[min]], "-C -g 3 -S 4 -m -1")
