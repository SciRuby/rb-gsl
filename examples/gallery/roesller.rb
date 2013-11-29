#!/usr/bin/env ruby
require("gsl")

dim = 3
roesller = Proc.new { |t, v, dvdt, params|
  a = params[0]; b = params[1]; c = params[2]
  x = v[0]; y = v[1]; z = v[2]
  dvdt[0] = - y - z
  dvdt[1] = x + a*y
  dvdt[2] = b*x - (c - x)*z
}

a = 0.344
b = 0.4
c = 4.5
solver = GSL::Odeiv::Solver.alloc(GSL::Odeiv::Step::RKF45, [1e-6, 0.0], roesller, dim)
solver.set_params(a, b, c)

t = 0.0; tend = 100.0
h = 1e-6
v = GSL::Vector.alloc(1, 0, 0)

GSL::ieee_env_setup()

IO.popen("gnuplot -persist", "w") do |io|
  io.print("set title 'Roesller equation'\n")
  io.print("set xlabel 'X'\n")
  io.print("set ylabel 'Y'\n")
  io.print("set zlabel 'Z'\n")
  io.printf("splot '-' u 2:3:4 w l\n")
  while t < tend
    t, h, status = solver.apply(t, tend, h, v)
    io.printf("%e %e %e %e\n", t, v[0], v[1], v[2])
    break if status != GSL::SUCCESS
  end
  io.printf("e\n")
  io.flush
end

