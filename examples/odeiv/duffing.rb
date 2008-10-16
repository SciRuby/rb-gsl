#!/usr/bin/env ruby
#
# Duffing equation
#
require("gsl")

dim = 2

# x'' + a x' + bx^3 = f cos(omega t)
#
duffing = Proc.new { |t, x, dxdt, params|
  a = params[0]; b = params[1]
  omega = params[2]; f = params[3]
  dxdt[0] = x[1]
  dxdt[1] = -a*x[1] - b*GSL::pow_3(x[0]) + f*Math::cos(omega*t)
}

solver = GSL::Odeiv::Solver.alloc(GSL::Odeiv::Step::RK8PD, [1e-6, 0.0], duffing, dim)

a = 0.05
b = 1.0
omega = 1.0
f = 5.0
solver.set_params(a, b, omega, f)

t = 0.0      # initial position
tend = 200
hstart = 2.0*Math::PI/1000

h = hstart*1.0

# Initial conditions, at x = 0
# P21(0) = 0, P21'(0) = 3
x = GSL::Vector.alloc(0.0, 0.0)
count = 0
IO.popen("graph -T X -C -g 3", "w") do |io|
  while t < tend
    t, h, status = solver.apply(t, tend, h, x)
    if count%2 == 0
      io.printf("%g %g\n", t, x[0])
    end
    break if status != GSL::SUCCESS
    count += 1
  end
end
