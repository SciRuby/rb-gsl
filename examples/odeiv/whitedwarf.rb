#!/usr/bin/env ruby
# This example solves Lane-Emden equation to calculate white dwarf
# structure, using Chandrasekhar's model, as a polytrope gas sphere
# of self-gravitating body supported by eletron degeneration pressure.

require("gsl")
include Math
include GSL::CONST::CGSM
include GSL::Odeiv

# Degenerated electron/neutron gas
# Reference: "Black holes, white dwarfs, and neutron stars"
#            S.L. Shapiro & S.A. Teukolsky, Willey, 1983
module Degenerate
  hbar = PLANCKS_CONSTANT_HBAR    # Planck's constant
  me = MASS_ELECTRON              # Electron mass
  mn = MASS_NEUTRON               # Neutron mass
  mu = UNIFIED_ATOMIC_MASS
  c = SPEED_OF_LIGHT
  ELambda = hbar/me/c             # Compton length of electron
  NLambda = hbar/mn/c
  MeC2 = me*c*c                   # Electron rest mass energy
  MnC2 = mn*c*c
  Factor_xe = GSL::pow(3.0*PI*PI/mu, 1.0/3.0)*ELambda
  Factor_xn = GSL::pow(3.0*PI*PI/mu, 1.0/3.0)*NLambda

# See Shapiro&Teukolsky Chapter 2
  def phi(x)
    tmp = sqrt(1.0 + x*x)
    (x*tmp*(2.0*x*x/3.0 - 1.0) + log(x + tmp))/8/PI/PI
  end

  def chi(x)
    tmp = sqrt(1.0 + x*x)
    (x*tmp*(1.0 + 2*x*x) - log(x + tmp))/8/PI/PI
  end

  def xe_rho_mue(rho, mue)
    Factor_xe*GSL::pow(rho/mue, 1.0/3.0)
  end

  def xn_rho(rho)
    Factor_xn*GSL::pow(rho, 1.0/3.0)
  end

end

# Polytrope gas sphere
module Polytrope

# Lane-Emden equation
#   n: polytrope index
  EmdenEq = Proc.new { |x, y, dydx, n|
    dydx[0] = y[1]
    dydx[1] = -GSL::pow(y[0], n) - 2.0/x*y[1]
  }

  def emden_xi(n)
    dim = 2
    y = GSL::Vector[1.0, 0.0]
    dydx = GSL::Vector.alloc(dim)

    solver = Solver.alloc(Step::RKF45, [1e-6, 0], EmdenEq, dim)
    solver.set_params(n)
    solver.reset

    vx = GSL::Vector.alloc(10000)
    vy = GSL::Vector.alloc(10000)
    vdy = GSL::Vector.alloc(10000)

    x = 0.0001
    xend = 10.0
    h = 1e-6

    file = File.open("polytrope.dat", "w")
    i = 0
    while x < xend
      x, h, status = solver.apply(x, xend, h, y)
      break if GSL::isnan?(y[0])
      break if GSL::isnan?(y[1])
      vx[i] = x
      vy[i] = y[0]
      vdy[i] = y[1]
      file.printf("%e %e %e %e\n", x, y[0], y[1], GSL::pow_3(y[0]))
      i += 1
      break if status != GSL::SUCCESS
      break if y[0] <= -0.1
    end
    file.close
#    p vx.size
#    p i
    vx2 = vx.subvector(0, i)
    vy2 = vy.subvector(0, i)
    vdy2 = vdy.subvector(0, i)
    spline = GSL::Spline.alloc(GSL::Interp::AKIMA, i)
    spline.init(vy2.reverse, vx2.reverse)

# Find star surface:
# Star sufrace is defined as the zero point of density structure function
    x1 = spline.eval(0.0)
    spline.init(vx2, vdy2)
    yx2 = spline.eval(x1).abs
    return [x1, yx2*x1*x1]

  end
end

# Chandrasekhar white dwarf model:
#   * Polytrope gas sphere
#   * Support its self gravity by electron degeneration pressure
class WhiteDwarf
  include Degenerate
  include Polytrope

  G = GRAVITATIONAL_CONSTANT
  @@x1 = nil
  @@x12 = nil

  def initialize(n, rhoc, mue)
    @n = n
    @gam = 1.0 + 1.0/n
    @rhoc = rhoc
    @mue = mue
    x = xe_rho_mue(@rhoc, mue)
    phix = phi(x)
    @K = MeC2*phix/GSL::pow_3(ELambda)/GSL::pow(@rhoc, @gam)
    @a = sqrt((@n + 1)*@K*GSL::pow(@rhoc, 1.0/@n - 1.0)/4/PI/G)
    if !@@x1
      @@x1, @@x12 = emden_xi(@n)
    end
    @Pc = MeC2*phix/GSL::pow_3(ELambda)*phix   # Electron Fermi pressure
    @radius = @a*@@x1                     # White dwarf radius
    @mass = 4.0*PI*GSL::pow_3(@a)*@rhoc*@@x12  # White dwarf mass
  end

  attr_accessor :radius
  attr_accessor :mass
  attr_accessor :rhoc    # Central density
  attr_accessor :Pc      # Central pressure
  attr_accessor :gam     # Adiabatic index
  attr_accessor :n       # Polytrope index
  attr_accessor :mue     # Mean eletron number per barion
  attr_accessor :K       # Normalization constant
  attr_accessor :a
end

file = File.open("whitedwarf.dat", "w")
rho = 1e6     # central density
while rho < 1e11
  wd = WhiteDwarf.new(3, rho, 2)
  file.printf("%e %e %e %e\n", rho, wd.Pc, wd.radius/1000/100, wd.mass/SOLAR_MASS)
  rho *= 1.2
end
file.close

system("gnuplot -persist whitedwarf.gp")
File.delete("polytrope.dat")
File.delete("whitedwarf.dat")
