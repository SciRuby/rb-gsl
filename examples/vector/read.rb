#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create test Vector v
v = GSL::Vector.alloc(4)

# Read data into Vector v from binary file a.dat using #fread
v.fread("a.dat")

# Show v
v

# Create another test Vector
v2 = GSL::Vector.alloc(4)

# Read data into Vector v2 from text file b.dat using #fscanf
v2.fscanf("b.dat")

# Show v2
v2
END
