#!/usr/bin/env ruby
require("gsl")
include Math

sint2 = GSL::Function.alloc { |t| Math::sin(Math::PI/2*t*t) }
cost2 = GSL::Function.alloc { |t| Math::cos(Math::PI/2*t*t) }
w = GSL::Integration::Workspace.alloc(1000)

t = 0
STDOUT.print("Computing... ")
STDOUT.flush
IO.popen("graph -T X -C -g 3 -X 'C(t)' -Y 'S(t)' --toggle-rotate-y-label -L 'Cornu spiral'", "w") do |io|
  t = -4
  while t < 4
    c = cost2.qag([0, t], w)[0]
    s = sint2.qag([0, t], w)[0]
    io.printf("%e %e\n", c, s)
    t += 0.01
  end
  STDOUT.print("done\n")
  STDOUT.flush
  io.flush
end
