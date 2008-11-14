#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create two 3x3 matrices
m = GSL::Matrix.alloc(1..9, 3, 3)
n = GSL::Matrix.alloc(11..19, 3, 3)

# Get Vector::Col::View for column 1 of Matrix a
a = m.col(1)

# Get Vector::Col::View for column 2 of Matrix b
b = n.col(2)

# Create new Matrix from Vector::Col::Views
c = GSL::Matrix[a, b]
END
