#!/usr/bin/env ruby
require("gnuplot")
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
solver = GSL::Odeiv::Solver.new(GSL::Odeiv::Step::RKF45, [1e-6, 0.0], roesller, dim)
solver.set_params(a, b, c)

t = 0.0; tend = 100.0
h = 1e-6
v = GSL::Vector[1, 0, 0]

GSL::ieee_env_setup()

N = 1000
x = GSL::Vector[N]
y = GSL::Vector[N]
z = GSL::Vector[N]

i = 0
while t < tend
  t, h, status = solver.apply(t, tend, h, v)
  x[i] = v[0]
  y[i] = v[1]
  z[i] = v[2]
  i += 1
  break if status != GSL::SUCCESS
end

Gnuplot::open do |gp|
  Gnuplot::SPlot.new( gp ) do |plot|

    plot.title "Roesller equation"

    plot.data = [
      Gnuplot::DataSet.new( [x.subvector(i), y.subvector(i), z.subvector(i)] ) { |ds|
        ds.title = ""
        ds.with = "lines"
      }
    ]
  end
end


