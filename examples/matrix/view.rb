#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create a 3x3 matrix m
m = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)

# Get Vector::View for row 1 of Matrix m
m.row(1)

# Get Vector::Col::View for column 0 of Matrix m
m.col(0)

# Get Vector::View for diagonal of Matrix m
m.diagonal

# Get Vector::View for subdiagonal 1 of Matrix m
m.subdiagonal(1)

# Get Vector::View for subdiagonal 0 of Matrix m
m.subdiagonal(0)

# Get Vector::View for superdiagonal 1 of Matrix m
m.superdiagonal(1)
END
