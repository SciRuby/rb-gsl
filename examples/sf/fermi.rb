#!/usr/bin/env ruby
require("gsl")
n = 100
x = GSL::Vector.linspace(0, 10, n)
m1 = GSL::Sf::fermi_dirac_m1(x)
IO.popen("graph -T X -C -g 3 -X x -L 'Fermi-Dirac function'", "w") do |io|
  for i in 0...n do
    io.printf("%e %e\n", x[i], m1[i])
  end
end
