#!/usr/bin/env ruby
#  Solve Schroedinger equation
#
#  This example is taken from frei1.cpp
#  in "Numerische Physik" p201-204 (Springer),
#  which simulates the time evolution of a probability density.
#
#  Name: frei1.cpp
#  Zweck: Simuliert ein quantenmechanisches freies Teilchen
#  Gleichung: Schroedingergleichung ohne Potential verwendete
#  Bibiliothek: GSL
#
#  Reference:
#    "Numerische Physik", by Harald Wiedemann, Springer (2004)
#    ISBN: 3-540-40774-X
#    http://www.springeronline.com/sgw/cda/frontpage/0,10735,1-102-22-22345455-0,00.html

require("gsl")

#NMAX = 8192
NMAX = 256

# The wave equation:
# calculate time derivative of the wave function.
# The second spatial derivative is approximated by
#   d2_psi/dx2 ~ (psi[n+1] - 2*psi[n] + pxi[n-1])/(dx*dx)
# See "Numerische Physik" p204, Eq. (5.47).
#
# psi(x), dpsi_dt(x): Complex-valued wavefunction, expressed as
#
#   0              NMAX             2*NMAX
#  |-----------------|-----------------|
#          Real           Imaginary
#
f = Proc.new { |t, psi, dpsi_dt|
  dx2 = $dx*$dx

# Real part
  for n in 1...(NMAX-1) do
    dpsi_dt[n] = -(psi[NMAX+n+1]+psi[NMAX+n-1]-2*psi[NMAX+n])/dx2
  end
  dpsi_dt[0] = -(psi[NMAX+1]+psi[2*NMAX-1]-2*psi[NMAX])/dx2
  dpsi_dt[NMAX-1] = -(psi[NMAX]+psi[2*NMAX-2]-2*psi[2*NMAX-1])/dx2

#  Imaginary part
  for n in (NMAX+1)...(2*NMAX-1) do
    dpsi_dt[n] = +(psi[n+1-NMAX]+psi[n-1-NMAX]-2*psi[n-NMAX])/dx2
  end
  dpsi_dt[NMAX] = +(psi[1]+psi[NMAX-1]-2*psi[0])/dx2
  dpsi_dt[2*NMAX-1] = +(psi[0]+psi[NMAX-2]-2*psi[NMAX-1])/dx2
}

psi = GSL::Vector[2*NMAX]
dpsi_dt = GSL::Vector[2*NMAX]

$dx = 0.1
dt = 0.1
n_out = 20
alpha = 1
p_0 = -0.5
atol = 1e-4
rtol = 0.0

h = 1.0e-4

dx2 = $dx*$dx
sum = 0.0
for n in 0...NMAX do
  x = (n-NMAX/2) * $dx
  psi[n] = Math::exp(-GSL::pow_2(x/alpha)/2)
  sum += GSL::pow_2(psi[n])
end
sum = 1.0/Math::sqrt(sum)

for n in 0...NMAX do
  x = (n-NMAX/2) * $dx
  psi[n+NMAX] = -psi[n] * sum * Math::sin(p_0*x) # Imaginaerteil
  psi[n] = psi[n] * sum * Math::cos(p_0*x)       # Realteil
end

IO.popen("graph -T X -C -g 3", "w") do |io|
  for n1 in 0...NMAX do
    x = (n1-NMAX/2) * $dx
    io.printf("%e %e\n", x, Math::sqrt(GSL::pow_2(psi[n1]) + GSL::pow_2(psi[n1+NMAX])))
  end
  io.printf("\n")

  step = GSL::Odeiv::Step.alloc(GSL::Odeiv::Step::RKF45, 2*NMAX)
  c = GSL::Odeiv::Control.y_new(atol, rtol)
  evolve = GSL::Odeiv::Evolve.alloc(2*NMAX)
  sys = GSL::Odeiv::System.alloc(f, 2*NMAX)

  t = 0.0
  for n in 1..n_out do
    t1 = n*dt
    STDOUT.printf("t = %2.1f (%2d/%2d)\n", t1-dt, n, n_out)
    while t < t1
      t, h, status = evolve.apply(c, step, sys, t, t1, h, psi)
      break if status != GSL::SUCCESS
    end
    if n%10 == 0
      for n1 in 0...NMAX do
        x = (n1-NMAX/2) * $dx
        io.printf("%e %e\n", x, Math::sqrt(GSL::pow_2(psi[n1]) + GSL::pow_2(psi[n1+NMAX])))
      end
      io.printf("\n")
    end
  end
end
