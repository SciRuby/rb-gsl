#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create Vector of 10 elements increasing logarithmically
# from 0 to 10000, inclusive, using Vector.logspace
v = GSL::Vector.logspace(0, 4, 10)

# Plot v
v.graph("-C -Y v -l y -S 4 -L 'Vector.logspace(0, 4, 10)'")

# Create Vector of 10 elements increasing logarithmically
# from 0 to 10000, inclusive, using Vector.logspace2
v = GSL::Vector.logspace2(1, 10000, 10)

# Plot v
v.graph("-C -Y v -l y -S 4 -L 'Vector.logspace2(1, 10000, 10)'")
END
