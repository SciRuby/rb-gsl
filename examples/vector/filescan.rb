#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Read file "c.dat" as doubles
a, b, c, d = GSL::Vector.filescan("c.dat")

# Read file "c.dat" as ints
a, b, c, d = GSL::Vector::Int.filescan("c.dat")
END
