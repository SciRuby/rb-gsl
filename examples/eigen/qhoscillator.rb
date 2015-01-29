#!/usr/bin/env ruby
#  Harmonic oscillator in quantum mechanics
#
#  This example is taken from "eigen1.cpp"
#  in "Numerische Physik" p214-216 (Springer).
#
#  Reference:
#    "Numerische Physik", by Harald Wiedemann, Springer (2004)
#    ISBN: 3-540-40774-X
#    http://www.springeronline.com/sgw/cda/frontpage/0,10735,1-102-22-22345455-0,00.html
require("gsl")

STDERR.puts("Harmonic oscillator in quantum mechanics:")

# Values used in the book:
#NMAX = 512
#dx = 0.02

# These are much faster:
NMAX = 64
dx = 0.16

# Potential
V = GSL::Vector[NMAX]
for n in 0...NMAX do
  x = (n - NMAX/2)*dx
  V[n] = 0.5*x*x
end

# Hamiltonian
H = GSL::Matrix.calloc(NMAX, NMAX)
H.set_diagonal(1.0/dx/dx + V)
tmp = -0.5/dx/dx
for n1 in 1...NMAX do
  H[n1-1,n1] = tmp
  H[n1,n1-1] = tmp
end
for n1 in 0...(NMAX-1) do
  H[n1+1,n1] = tmp
  H[n1,n1+1] = tmp
end

# Calculate eigen values and eigen vectors
STDERR.print("  Solving the eigen system of #{NMAX}X#{NMAX} dimensions...")
STDERR.flush
eval, evec = H.eigen_symmv
GSL::Eigen.symmv_sort(eval, evec, GSL::Eigen::SORT_VAL_ASC)
STDERR.puts("OK")
STDERR.flush

x2 = GSL::Vector[NMAX]
for n1 in 0...NMAX do
  x2[n1] = 0
  for n2 in 0...NMAX do
    x = (n2 - NMAX/2)*dx
    x2[n1] += GSL::pow_2(evec[n2,n1]*x)
  end
end

# Energy eigen values, see p217 "Tabelle 5.1"
# The differences with Tabelle 5.1 are from NMAX and dx.
# If we use NMAX=512 and dx=0.02, we obtain the same results (but much slower).
STDERR.puts("  Eigen values:")
STDERR.printf("  %2s Exact %5s %10s | %2s Exact %5s %10s\n",
              "n", "E", "err\(\%\)", "n", "E", "err\(\%\)")
STDERR.print("  -----------------------------------------------------\n")
for n1 in 0..6 do
  exact1 = n1 + 0.5
  exact2 = n1 + 7 + 0.5
  STDERR.printf("  %2d %4.1f %8.5f %+7.5f | %2d %4.1f %8.5f %+7.5f\n",
                n1, exact1, eval[n1], (exact1 - eval[n1])/exact1*100,
                n1+7, exact2, eval[n1+7], (exact2-eval[n1+7])/exact2*100)
end
STDERR.flush

# Eigen vectors of n = 0, 1, 2, 10. See p217 "Abb 5.3"
c = Math::sqrt(1.0/dx)
vec0 = evec.col(0).scale(c)
vec1 = evec.col(1).scale(c)
vec2 = evec.col(2).scale(c)
vec10 = evec.col(10).scale(c)
File.open("qhoscillator.dat", "w") do |fp|
  for i in 0...NMAX do
    x = (i - NMAX/2)*dx
    fp.printf("%e %e %e %e %e\n", x, -vec0[i], vec1[i], -vec2[i], -vec10[i])
  end
end
system("gnuplot -persist qhoscillator.gp")
File.delete("qhoscillator.dat")

