#!/usr/bin/env ruby
require("gsl")
n = 100
x = GSL::Vector.linspace(0.01, 4, n)
I0 = GSL::Sf::bessel_I0(x)
I1 = GSL::Sf::bessel_I1(x)
I2 = GSL::Sf::bessel_In(2, x)

y0 = GSL::Sf::bessel_K0(x)
y1 = GSL::Sf::bessel_K1(x)
y2 = GSL::Sf::bessel_Kn(2, x)
GSL::Vector.graph(x, I0, I1, I2, y0, y1, y2, "-T X -C -g 3 -y 0 6 -X x -L 'Bessel functions I0, I1, I2, K0, K1, K2'")
