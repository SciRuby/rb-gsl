#!/usr/bin/env ruby
#  This solves the Saha's equation, to investigate
#  'recombination' of hot cosmic plasma.
require("gsl")
include GSL
include GSL::CONST::CGSM
include Math

KB =  BOLTZMANN
HBAR = PLANCKS_CONSTANT_HBAR

class Recombination
  def initialize(t0, bdensity, eion)
    @t = 2500
    @tCMB = t0                   # Present temperature of CMB [K]
    @n0 = bdensity/MASS_PROTON   # Barion number density [cm-3]
    @eion = eion*ELECTRON_VOLT   # Hydrogen ionization potential [eV]
    @saha_equation = Function.alloc { |x|   # x: fractional ionization
      n = @n0*pow_3(@t/@tCMB)
      tmp1 = pow(MASS_ELECTRON*KB*@t/2.0/PI/HBAR/HBAR, 1.5)
      tmp2 = exp(-@eion/KB/@t)
      x*x - (1.0 - x)*tmp1*tmp2/n
    }
  end

  def zplus1()
    @t/@tCMB
  end

  attr_accessor :t
  attr_reader :saha_equation
end

###########################
## Cosmological parameter
TCMB = 2.725        # Temperature of cosmic microwave background at present [K]
BDENSITY = 1e-29    # Barion density [g cm-3]
EIONIZE = 13.6      # Hydrogen ionization potential [eV]
unv = Recombination.new(TCMB, BDENSITY, EIONIZE)

# FSolver
solver = Root::FSolver.alloc(Root::FSolver::BRENT)

TMIN = 2500
TMAX = 6500

begin
  file = File.open("recombination.dat", "w")
  t = TMIN
  while t < TMAX
    unv.t = t  # Temperature
    x, iter, status = solver.solve(unv.saha_equation, [1e-8, 1.0], [0, 1e-4])
    file.printf("%e %e %e\n", unv.t, unv.zplus1, x)
    t *= 1.01
  end
ensure
  file.close
end
system("gnuplot -persist recombination.gp")
File.delete("recombination.dat")
#puts("Try \"gnuplot -persist recombination.gp\"")
