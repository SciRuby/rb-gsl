#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create Vector of 10 elements increasing linearly
# from 0 to 5, inclusive, using Vector.linspace
v = GSL::Vector.linspace(0, 5, 10)

# Plot v
v.graph("-C -Y v -S 4 -L 'Vector.linspace(0, 5, 10)'")

# Create Vector of 11 elements from 0 to 5, inclusive, using Vector.linspace
v = GSL::Vector.linspace(0, 5, 11)

# Plot v
v.graph("-C -Y v -S 4 -L 'Vector.linspace(0, 5, 11)'")
END
