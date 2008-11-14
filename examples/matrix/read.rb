#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create a 2x2 test Matrix m
m = GSL::Matrix.alloc(2, 2)

# Read data into Matrix m from binary file a.dat using #fread
m.fread("a.dat")

# Show m
m

# Create a 2x2 test Matrix m2
m2 = GSL::Matrix.alloc(2, 2)

# Read data into Matrix m2 from text file b.dat using #fscanf
m2.fscanf("b.dat")

# Show m2
m2
END
