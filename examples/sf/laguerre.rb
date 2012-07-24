#!/usr/bin/env ruby
require("gsl")
x = GSL::Vector.linspace(0, 5, 40)
e = GSL::Sf::exp(-x/2)

a = 0
IO.popen("graph -T X -C -g 3 -L 'Laguerre functions'", "w") do |io|
  for n in 0..5 do
    for i in 0...20 do
      y = GSL::Sf::laguerre_n(n, a, x[i])*e[i]
      io.printf("%e %e\n", x[i], y)
    end
    io.printf("\n")
  end
end




