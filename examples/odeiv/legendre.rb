#!/usr/bin/env ruby
#
# Solve Legendre's differential equation
#   l = 2, m = 1
require("gsl")
dim = 2

fleg = Proc.new { |x, y, f, params|
  l = params[0]
  m = params[1]
  f[0] = y[1]
  f[1] = 2.0*x/(1-x*x)*y[1] - (l*(l + 1) - m*m/(1-x*x))/(1-x*x)*y[0]
}

solver = GSL::Odeiv::Solver.alloc(GSL::Odeiv::Step::RKF45, [1e-6, 0.0], fleg, dim)

# P21
l = 2
m = 1
solver.set_params(l, m)

x = 0.0      # initial position
xend = 0.999
hstart = 1e-8

h = hstart*1.0

# Initial conditions, at x = 0
# P21(0) = 0, P21'(0) = 3
y = GSL::Vector.alloc(0.0, 3.0)

File.open("legode.dat", "w") do |f|
  while x < xend
    x, h, status = solver.apply(x, xend, h, y)
    f.printf("%g %g\n", x, y[0])
    break if status != GSL::SUCCESS
  end
end

File.open("plm.dat", "w") do |f|
  x = 0.0
  while x < xend
    plmx = GSL::Sf::legendre_Plm(l, m, x).abs
    f.printf("%g %g\n", x, plmx)
    x += 0.01
  end
end

system("graph -T X -C -g 3 -X x -Y 'P21(x)' --toggle-rotate-y-label --title-font-size 0.05 -L 'Red: expect, Green: RKF45' -m 1 plm.dat -m -2 -S 4 legode.dat")
File.delete("legode.dat")
File.delete("plm.dat")

