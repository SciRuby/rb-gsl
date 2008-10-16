#!/usr/bin/env ruby
require("gsl")

y1 = GSL::Vector[0.2, 3, 2, 0.4, 0, 4, 0]
y1.graph

x = GSL::Vector.linspace(0, 10, 50)
y2 = GSL::Sf::bessel_J0(x)
y2.graph(x)

GSL::Vector.graph([nil, y1], [x, y2])
