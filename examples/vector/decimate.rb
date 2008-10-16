#!/usr/bin/env ruby
require("gsl")

N = 1000
DECIMATE1 = 10
DECIMATE2 = 100
r = GSL::Rng.alloc
x0 = GSL::Vector.linspace(0, 20, N)

# Data: Bessel function + noise
y0 = GSL::Sf::bessel_J0(x0) + GSL::Ran::gaussian(r, 0.1, N)

y1 = y0.decimate(DECIMATE1)
y2 = y0.decimate(DECIMATE2)

x1 = GSL::Vector.linspace(0, 20, N/DECIMATE1)
x2 = GSL::Vector.linspace(0, 20, N/DECIMATE2)

# y1 and y2 are shifted vertically for visual purpose
GSL::graph([x0, y0], [x1, y1-1], [x2, y2-2])  

#graph(x2, y2-2)

#a = y2 - 2
#for i in 0...x2.size do
#  printf("%f %f\n", x2[i], a[i])
#end
	
