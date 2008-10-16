#!/usr/bin/env ruby
require("gsl")

a = 1
t = GSL::Vector.linspace(0, 3*Math::PI, 100)
cost = GSL::Sf::cos(t/3)
cost3 = cost*cost*cost
r = 4*a*cost3
x = r*GSL::Sf::cos(t)
y = r*GSL::Sf::sin(t)

GSL::graph(x, y, "-T X -C")
