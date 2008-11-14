#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create three test Vectors
v1 = GSL::Vector[1..2]
v2 = GSL::Vector[3..4]
v3 = GSL::Vector[5..7]

# Connect them using GSL::Vector#connect
v1.connect(v2, v3)

# Connect them using GSL::Vector.connect
GSL::Vector.connect(v1, v2, v3)
END
