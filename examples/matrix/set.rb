#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create a 3x3 matrix
m = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)
# Create a 3x4 matrix
m = GSL::Matrix.indgen(3, 4)

# Set element at row 1, column 2 to 99.9
m[1,2] = 99.9

# Show matrix
m

# Set all elements to 5 using #set_all
m.set_all(5)

# Set all elements to 4.3 using #set
m.set(4.3)

# Set all elements to 2 using #[]
m[] = 2

# Show matrix
m

# Set all elements to 0
m.set_zero

# Set matrix to identity matrix
m.set_identity
END
