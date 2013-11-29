#!/usr/bin/env ruby
require("gsl")
n = 50
x = GSL::Vector.linspace(0.01, 14, n)
y05 = GSL::Sf::gamma_inc_P(0.5, x)
y1 = GSL::Sf::gamma_inc_P(1, x)
y3 = GSL::Sf::gamma_inc_P(3, x)
y10 = GSL::Sf::gamma_inc_P(10, x)
GSL::Vector.graph(x, y05, y1, y3, y10, "-T X -C -g 3 -y 0 1.1 -L 'Incomplete gamma functions'")
