#!/usr/bin/env ruby
require("gsl")
include Math

def ran_power_law(slope, range, rng)
  x = rng.uniform
  xmin = GSL::pow(10.0, range*(1.0 + slope));
  x2 = (1 - xmin)*x + xmin;
  return GSL::pow(x2, 1.0/(1+slope));
end

h = GSL::Histogram.alloc(19, [1, 20])
rng = GSL::Rng.alloc
for i in 0..5000
  x = ran_power_law(-3.2, 1000, rng)
  h.increment(x)
end

result = h.fit_power
x = GSL::Vector.logspace2(1, 20, 19)
y = result[0]*GSL::pow(x, result[1])
GSL::graph(h, [x, y], "-l x -l y -x 1 20 -y 1 10000 -C -g 3")

puts("Expected power index: -3.2")
printf("Estimated power index: %5.4f +/- %5.4f\n", result[1], result[3])
