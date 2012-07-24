#!/usr/bin/env ruby
require("gsl")
IO.popen("graph -T X -C -g 3 -X x -Y 'psi(x)' --toggle-rotate-y-label -x -3 5 -y -5 15 -L 'Red: Digamma, Green: Trigamma'", "w") do |io|
  x = -2.999
  while x < -2
    y = GSL::Sf::psi(x)
    io.printf("%e %e\n", x, y)
    x += 0.01
  end
  x = -1.999
  while x < -1
    y = GSL::Sf::psi(x)
    io.printf("%e %e\n", x, y)
    x += 0.01
  end
  x = -0.999
  while x < 0
    y = GSL::Sf::psi(x)
    io.printf("%e %e\n", x, y)
    x += 0.01
  end
  x = 0.001
  while x < 5
    y = GSL::Sf::psi(x)
    io.printf("%e %e\n", x, y)
    x += 0.1
  end

  io.printf("\n")
  x = -2.999
  while x < -2
    y = GSL::Sf::psi_1_e(x).val
    io.printf("%e %e\n", x, y)
    x += 0.01
  end
  x = -1.999
  while x < -1
    y = GSL::Sf::psi_1_e(x).val
    io.printf("%e %e\n", x, y)
    x += 0.01
  end
  x = -0.999
  while x < 0
    y = GSL::Sf::psi_1_e(x).val
    io.printf("%e %e\n", x, y)
    x += 0.01
  end
  x = 0.001
  while x < 5
    y = GSL::Sf::psi_1_e(x).val
    io.printf("%e %e\n", x, y)
    x += 0.1
  end
end
