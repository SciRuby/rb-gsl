#!/usr/bin/env ruby
require("gsl")
n = 100
x = GSL::Vector.linspace(0.01, 20, n)
J0 = GSL::Sf::bessel_J0(x)
J1 = GSL::Sf::bessel_J1(x)
J2 = GSL::Sf::bessel_Jn(2, x)

y0 = GSL::Sf::bessel_Y0(x)
y1 = GSL::Sf::bessel_Y1(x)
y2 = GSL::Sf::bessel_Yn(2, x)

GSL::Vector.graph(x, J0, J1, J2, y0, y1, y2, "-T X -C -g 3 -y -2 1.1 -X x -L 'Bessel functions J0, J1, J2, Y0, Y1, Y2'")
