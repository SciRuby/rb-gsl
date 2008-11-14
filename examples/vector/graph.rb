#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create Vector to plot
y1 = GSL::Vector[0.2, 3, 2, 0.4, 0, 4, 0]

# Plot it using #graph()
y1.graph

# Create x and y2 Vectors
x = GSL::Vector.linspace(0, 10, 50)

y2 = GSL::Sf::bessel_J0(x)

# Plot x vs y using #graph(x)
y2.graph(x)

# Plot y1 and x vs y2 using GSL::Vector.graph
GSL::Vector.graph([nil, y1], [x, y2])
END
