#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create Vector of x values
x = GSL::Vector.linspace(1, 5, 5)

# Create Vector y = x**5
y = GSL::pow_5(x)

# Compute succesive differences
y1 = y.diff
y2 = y1.diff
y3 = y2.diff
y4 = y3.diff

# Show that successive differences can be computed directly from y
y.diff(2) == y2
y.diff(3) == y3
y.diff(4) == y4

# Plot x and y Vectors
GSL::graph(x, y, y1, y2, y3, "-C -g 3 -l x -l y -x 1 5 1")
END
