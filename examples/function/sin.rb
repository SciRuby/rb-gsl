#!/usr/bin/env ruby
require("gsl")

x = GSL::Vector.linspace(0, 2*Math::PI, 20)
f = GSL::Function.alloc { |x|
  GSL::Sf::sin(x)
}

f.graph(x, "-T X -g 3 -C -L 'sin(x)'")

