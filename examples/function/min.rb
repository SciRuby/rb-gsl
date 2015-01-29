#!/usr/bin/env ruby
require("gsl")

humps = GSL::Function.alloc { |x|
  1.0/(GSL::pow_2(x-0.3) + 0.01) + 1.0/(GSL::pow_2(x-0.9)+0.04) - 6
}

iter = 0;  max_iter = 500
m = 0.7             # initial guess
m_expected = 0.6370
a = 0.3; b = 1.0
gmf = GSL::Min::FMinimizer.alloc(GSL::Min::FMinimizer::BRENT)
gmf.set(humps, m, a, b)
printf("Using %s method\n", gmf.name)
printf("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "min",
       "err", "err(est)")
printf("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
       iter, a, b, m, m - m_expected, b - a)
begin
  iter += 1
  status = gmf.iterate
  status = gmf.test_interval(0.0001, 0.0)
  puts("Converged:") if status == GSL::SUCCESS
  a = gmf.x_lower;  b = gmf.x_upper
  m = gmf.x_minimum
  printf("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
         iter, a, b, m, m - m_expected, b - a);
end while status == GSL::CONTINUE and iter < max_iter

x = GSL::Vector.linspace(-0.5, 1.5, 100)
GSL::graph([x, humps.eval(x)], [GSL::Vector[m], GSL::Vector[humps.eval(m)]], "-C -g 3 -S 4 -m -1")


