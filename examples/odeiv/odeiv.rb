#!/usr/bin/env ruby
require("gsl")
dim = 2
func = Proc.new { |t, y, dydt, mu|
  dydt[0] = y[1]
  dydt[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1.0)
}
step = GSL::Odeiv::Step.alloc(GSL::Odeiv::Step::RKF45, dim)
c = GSL::Odeiv::Control.y_new(1e-6, 0.0)
evolve = GSL::Odeiv::Evolve.alloc(dim)
sys = GSL::Odeiv::System.alloc(func, dim)
mu = 10.0
sys.set_params(mu)
t = 0.0; tend = 100.0
h = 1e-6
y = GSL::Vector[1.0, 0.0]
GSL::ieee_env_setup()

N = 1500
tt = GSL::Vector[N]
yt = GSL::Vector[N]
i = 0
while t < tend and i < N
  t, h, status = evolve.apply(c, step, sys, t, tend, h, y)
  tt[i] = t
  yt[i] = y[0]
  break if status != GSL::SUCCESS
  i += 1
end
yt.subvector(i).graph(tt.subvector(i), "-T X -C -g 3 -x 0 #{tt.max} -L '#{step.name}' -S 4")


