#!/usr/bin/env ruby
require("gsl")
include Math

# Function to be expanded
G = GSL::Function.alloc { |x| exp(-x*x) }
#G = Function.alloc { |x| sin(-x*x) }

# Sampling points
XMIN = 0.0
XMAX = 3.0
SIZE = 20

x = GSL::Vector.linspace(XMIN, XMAX, SIZE)
sample = G.eval(x)

# Discrete Hankel transform with the Bessel function J0
dht = GSL::Dht.alloc(SIZE, 0, XMAX)
g = dht.apply(sample)

num = dht.num
den = dht.den
coef = dht.coef

# Reconstruction
y = GSL::Vector[SIZE]
for n in 0...SIZE do
  val = 0.0
  for m in 0...SIZE do
    a = GSL::Sf::bessel_J0(dht.sample(n, m))
    val += (2.0/XMAX/XMAX)*a/den[m]*g[m]
#    val += (2.0/XMAX/XMAX)*num[n][m]/den[m]*g[m]
#    val += coef[n][m]*g[m]
  end
  y[n] = val
end

x0 = GSL::Vector.linspace(XMIN, XMAX, 100)
y0 = G.eval(x0)

GSL::graph([x0, y0], [x, sample], [x, y], "-T X -C -g 3 -X t -Y 'f(t)' --toggle-rotate-y-label -L 'Red: f(t), Green: sample, Blue: DHT of size #{SIZE}'")

