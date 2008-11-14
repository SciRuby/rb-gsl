#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create test vector v
v = GSL::Vector.alloc(1, 2, 3, 4, 5, 6, 7, 8, 9)

# Create Vector::View of Vector v starting at element 2 and spanning 3 elements
vv = v.subvector(2, 3)

# Set element 1 of View to 9.0
vv.set(1, 9)

# Show that Vector is also modified
v

# Create Matrix::View of Vector
m = v.matrix_view(3, 3)

# Create Vector::View spanning the entire Vector
v2 = v.view

# Create Vector::View of Vector starting at element 0 and spanning 3 elements
v2 = v.view(3)

# Show size of Vector::View
v2.size
END
