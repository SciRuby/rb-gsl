#!/usr/bin/env ruby
require("gsl")

N = 1000
GSL::Rng.env_setup()
T = GSL::Rng::DEFAULT
r = GSL::Rng.alloc(T)

x = 0.0
y = 0.0

printf("%g %g\n", x, y);

N.times do
  dx, dy = r.dir_2d()
  x += dx
  y += dy
  printf("%g %g\n", x, y);
end
