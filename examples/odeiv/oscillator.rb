#!/usr/bin/env ruby
# An oscillator
#   m: mass
#   k: spring constant
#   b: resist
#   f: external force

require("gsl")

dim = 2

# x[0]: displacement, x[1]: velocity
func = Proc.new { |t, x, dxdt, params|
  m = params[0]; b = params[1]; f = params[2]; k = params[3]
  dxdt[0] = x[1]
  dxdt[1] = (f - b*x[1] - k*x[0])/m
}

gos = GSL::Odeiv::Solver.alloc(GSL::Odeiv::Step::RKF45, [1e-6, 0.0], func, dim)

m = 1.0
b = 1.0
f = 1.0
k = 10.0

gos.set_params(m, b, f, k)

t = 0.0; tend = 10.0
h = 1e-6
x = GSL::Vector.alloc([0.0, 0.0])

GSL::ieee_env_setup()

IO.popen("graph -T X -C -g 3 -S 4 -X Time -Y Amp", "w") do |io|
  while t < tend
    t, h, status = gos.apply(t, tend, h, x)
    break if status != GSL::SUCCESS
#    printf("%.5e %.5e %.5e %.5e\n", t, x[0], x[1], h)
    io.printf("%.5e %.5e\n", t, x[0])
  end
end

