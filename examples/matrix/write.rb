#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create test Matrix
m = GSL::Matrix.alloc([1, 2], [3, 4])

# Write Matrix in binary format to file a.dat
m.fwrite("a.dat")

# Create another test Matrix
m2 = GSL::Matrix.alloc([5, 6], [7, 8])

# Write Matrix in text format to file b.dat
m2.fprintf("b.dat")
END
