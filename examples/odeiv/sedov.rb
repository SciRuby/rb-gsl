#!/usr/bin/env ruby
# Self-similar solution of the Sedev-Taylor equation,
# which describes spherically expanding gas
# triggered by a supernova explosion.

require("gsl")
include Math

# Sedov-Taylor equation
#     logx: (log of) the radial coordinate in the gas shpere, x = r/R.
#           logx = 0 at the surface, and negative within the shell
#        y: GSL::GSL::Vector::View::ReadOnly
#           y[0]: gas velocity at the radius
#           y[1]: density
#           y[2]: pressure
#  dydlogx: GSL::GSL::Vector::View, derivatives
#       sh: specific heat of the gas (5/3)
Sedov = Proc.new { |logx, y, dydlogx, sh|
  a = -5.0*(6.0*y[2] - 15.0*y[2]*y[0]*sh + 2.0*y[0]*y[1] - 7.0*y[0]*y[0]*y[1] + 5.0*GSL::pow_3(y[0])*y[1])
  b = -25.0*y[2]*sh + 4.0*y[1] - 20.0*y[0]*y[1] + 25.0*y[0]*y[0]*y[1]
  dydlogx[0] = a/b

  bb = -5.0*(6.0*y[2] - 15.0*y[2]*y[0]*sh + 2.0*y[0]*y[1] - 7.0*y[0]*y[0]*y[1])

  a = -5.0*(-30.0*y[2]*y[1] + 2.0*y[0]*y[1]*y[1] - 25.0*GSL::pow_2(y[0]*y[1]) + 50.0*GSL::pow_3(y[0])*GSL::pow_2(y[1]))
  bb = (-2.0 + 5.0*y[0])*(-25.0*y[2]*sh + 4.0*y[1] - 20.0*y[0]*y[1] + 25.0*y[0]*y[0]*y[1])
  dydlogx[1] = a/bb

  a = -5.0*y[2]*(-10.0*y[2]*sh + 4.0*y[1] - 14.0*y[0]*y[1] + 10.0*y[0]*y[0]*y[1] - y[0]*sh*y[1] + 10.0*y[0]*y[0]*sh*y[1])
  dydlogx[2] = a/b
}

DIM = 3
sh = 5.0/3.0
solver = GSL::Odeiv::Solver.alloc(GSL::Odeiv::Step::RKF45, [1e-5, 1e-5], Sedov, nil, DIM)
solver.set_params(sh)

def boundary_condition(sh)
  delta = 2.0/5.0
  xx = 1
  y = GSL::Vector[DIM]
  y[0] = 2.0/(sh + 1.0)*delta
  y[1] = (sh + 1.0)/(sh - 1.0)
  y[2] = 2.0/(sh + 1.0)*delta*delta
  y00 = xx*y[0]
  y10 = y[1]
  y20 = xx*xx*y[2]
  c = adiabatic_integ(xx, y, sh);
  return y, y00, y10, y20, c
end

def adiabatic_integ(xx, y, sh)
  GSL::pow(y[1], 1.0 - sh)*y[2]*(y[0] - 2.0/5.0)*GSL::pow_5(xx)
end

# y: GSL::Vector (velocity, density, pressure)
# y00, y10, y20: initial values
# c: adiabatic integral (this must be constant throughout the computation)
y, y00, y10, y20, c = boundary_condition(sh)

# x: the radial coordinates in the shell, x = r/R.
#     x = 1 at the surface (r = R), and zero at the center (r = 0)
x = 1
logx = log(x)
logxend = log(1e-5)
h = -0.000001

N = 150
R = GSL::Vector[N]
V = GSL::Vector[N]
RHO = GSL::Vector[N]
PRESS = GSL::Vector[N]

# The values are normalized to [0, 1]
R[0] = x
V[0] = x*y[0]/y00
RHO[0] = y[1]/y10
PRESS[0] = x*x*y[2]/y20

GSL::ieee_env_setup()

n = 1
while logx > logxend
  logx, h, status = solver.apply(logx, logxend, h, y)
  break if status != GSL::SUCCESS
  break if y[0] < 0.0
  x = exp(logx)
  c = adiabatic_integ(x, y, sh);
  R[n] = x
  V[n] = x*y[0]/y00
  RHO[n] = y[1]/y10
  PRESS[n] = x*x*y[2]/y20
  n += 1
end

GSL::graph(R.subvector(n), V.subvector(n), RHO.subvector(n), PRESS.subvector(n),
      "-T X -C -g 3 -X 'x = r/R' -S 4 -L 'red: v, green: rho, blue: p'")
