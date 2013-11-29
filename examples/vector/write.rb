#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create test Vector
v = GSL::Vector.alloc([1, 2, 3, 4])

# Write Vector in binary format to file a.dat
v.fwrite("a.dat")

# Create another test Vector
v2 = GSL::Vector.alloc([5, 6, 7, 8])

# Write Vector in text format to file b.dat
v2.fprintf("b.dat")
END

