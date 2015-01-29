#!/usr/bin/env ruby
#      19/Apr/2004         by Yoshiki Tsunesada
#
#   This is an example to calculate the orbital evolution of
# a double neutron star (binary) system. General relativity predicts
# that the binary orbital decays by radiating gravitational waves,
# and the two stars will coalesce in time scale of 100-1000 Mega-years.
#   The values used here are of the binary system J0730-3039 discovered
# in 2003 (Burgay et al., Nature 2003). The result shows that the two
# neutron stars will merge after about 85 Mega-years. From the age of
# the system 100 Mega-year, the lifetime of the system is estimated
# about 185 Mega-years.
#
# References:
#   1. Burgay et al., Nature 426, 531 (2003)
#   2. Shapiro & Teukolsky, "Black holes, white dwarfs and neutron stars"
#        John Wiley and Sans (1983)
#

require("gsl")
include Math

class BinarySystem
  def initialize(m1, m2)
    @m1 = m1
    @m2 = m2
  end
  attr_reader :m1
  attr_reader :m2
end

GMsolarC3 = 4.925490947e-6
MegaYear = 3600*24*365*1e6

# Time evolution of the binary orbital period and the eccentricity
# due to gravitational radiation.
# The calculation is based on general relativity (See e.g. Ref.2).
#     y[0]: orbital period (pb)
#     y[1]: eccentricity (e)
#  dydt[0]: time derivative of pb
#  dydt[1]: time derivative of e

deriv = Proc.new { |t, y, dydt, binary|
  pb = y[0]            # orbital period
  e = y[1]             # eccentricity
  m1 = binary.m1       # neutron star masses
  m2 = binary.m2
  totalM = m1 + m2     # total mass
  mu = m1*m2/totalM    # reduced mass
  mm = mu*GSL::pow(totalM, 2.0/3.0)
  f_e = GSL::pow(1.0 - e*e, -3.5)*(1.0 + (73.0/24.0 + 37.0/96.0*e*e)*e*e);
  h_e = (1.0 + 121.0/304.0*e*e)*GSL::pow(1.0 - e*e, -2.5)
  tmp = GSL::pow(GMsolarC3*2.0*PI/pb, 5.0/3.0)
  dydt[0] = -192.0*PI/5.0*f_e*tmp*mm              # dP/dt
  dydt[1] = -304.0/15.0*e*h_e*tmp*(2.0*PI/pb)     # de/dt
}

# Neutron star masses in solar-mass unit.
# The values are of the binary system J0730-3039 discoverd in 2003.
# See Burgay et al., Nature 426, 531 (2003)
m1 = 1.34
m2 = 1.24
#m1 = 1.25
#m2 = 2.574 - m1
binary = BinarySystem.new(m1, m2)

# Initial data: the present values
pb = 2.45*3600        # orbital period: 2.45 hours at present
ecc = 0.088           # eccentricity
#pb = 7.67*3600
#ecc = 0.181
y = GSL::Vector[pb, ecc]

# ODE solver using RKF45 algorithm
solver = GSL::Odeiv::Solver.alloc(GSL::Odeiv::Step::RKF45, [1e-6, 0.0], deriv, 2)
solver.set_params(binary)

# the age of the binary system from birth
age = 100*MegaYear
#age = 444*MegaYear
t = 0
tend = 2500*MegaYear

# initial time step
h = 1.0*MegaYear

begin
  file = File.open("binarysystem.dat", "w")
  while t < tend
    t, h, status = solver.apply(t, tend, h, y)
    break if status != GSL::SUCCESS
    break if GSL::isnan?(y[0])
    file.printf("%e %e %e %e\n", (t+age)/MegaYear, y[0]/3600, y[1], h/MegaYear)
  end
ensure
  file.close
end

system("gnuplot -persist binarysystem.gp")
File.delete("binarysystem.dat")

__END__


