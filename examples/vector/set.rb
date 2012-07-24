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
v = GSL::Vector.alloc(1, 2, 3, 4, 5)

# Show third element
v[3]

# Set all elements to 9
v.set_all(9)

# Set all elements to 0
v.set_zero

# Set all elements to 3
v[] = 3

# Show v
v

# Set all elements to 0, except element 3 which is set to 1
v.set_basis(3)

# Set element 2 to 5.0
v.set(2, 5)
END
