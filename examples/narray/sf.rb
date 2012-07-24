#!/usr/bin/env ruby
require("gsl")

na = NArray[0.1, 0.2, 0.3, 0.4]
p GSL::Sf::legendre_Pl(2, na)

v = GSL::Vector[0.1, 0.2, 0.3, 0.4]
p GSL::Sf::legendre_Pl(2, v)

na = NArray[[1.0, 2, 3, 4], [2, 3, 4, 5]]
p GSL::Sf::sin(na)

m = GSL::Matrix[[1.0, 2, 3, 4], [2, 3, 4, 5]]
p GSL::Sf::sin(m)

n = 50
x = GSL::Vector.linspace(0.01, 1, n).to_na
y1 = GSL::Sf::beta_inc(0.5, 5.0, x)
y2 = GSL::Sf::beta_inc(1.0, 3.0, x)
y3 = GSL::Sf::beta_inc(8.0, 10.0, x)
y4 = GSL::Sf::beta_inc(5.0, 0.5, x)
p y1.class
p y2.class
p y3.class
p y4.class
GSL::Vector.graph([x, y1], [x.to_gv, y2], [x, y3.to_gv], [x, y4], "-T X -C -g 3 -y 0 1.1 -X x -L 'Incomplete beta functions'")
