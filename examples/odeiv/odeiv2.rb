#!/usr/bin/env ruby
require("gsl")

dim = 2

func = Proc.new { |t, y, dydt, mu|
  dydt[0] = y[1]
  dydt[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1.0)
}

jac = Proc.new { |t, y, dfdy, dfdt, mu|
  dfdy[0][0] = 0.0
  dfdy[0][1] = 1.0
  dfdy[1][0] = -2*mu*y[0]*y[1] - 1.0
  dfdy[1][1] = -mu*(y[0]*y[0] - 1.0)
  dfdt[0] = 0.0
  dfdt[1] = 0.0
}

#solver = GSL::Odeiv::Solver.alloc(GSL::Odeiv::Step::RK8PD, [1e-6, 0.0, 1, 0], func, jac, dim)
#solver = GSL::Odeiv::Solver.alloc(GSL::Odeiv::Step::RK8PD, [1e-6, 0.0, 1, 0], func, nil, dim)
#solver = GSL::Odeiv::Solver.alloc(GSL::Odeiv::Step::RK8PD, [1e-6, 0.0], func, nil, dim)
#solver = GSL::Odeiv::Solver.alloc(GSL::Odeiv::Step::RKF45, [1e-6, 0.0], func, jac, dim)
#solver = GSL::Odeiv::Solver.alloc(GSL::Odeiv::Step::RKF45, [1e-6, 0.0], func, nil, dim)
solver = GSL::Odeiv::Solver.alloc(GSL::Odeiv::Step::BSIMP, [1e-6, 0.0], func, jac, dim)
solver.set_params(10.0)

t = 0.0; tend = 100.0
h = 1e-6
y = GSL::Vector.alloc([1.0, 0.0])

GSL::ieee_env_setup()
N = 1500
tt = GSL::Vector[N]
yt = GSL::Vector[N]
i = 0
while t < tend and i < N
  t, h, status = solver.apply(t, tend, h, y)
  break if status != GSL::SUCCESS
  tt[i] = t
  yt[i] = y[0]
  i += 1
end

GSL::graph(tt.subvector(i), yt.subvector(i),  "-T X -C -g 3 -x 0 #{tt.max} -L '#{solver.step.name}' -S 4")
