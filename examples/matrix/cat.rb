#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create a 2x2 Matrix::Int and a 2x3 Matrix::Int
a = GSL::Matrix::Int[1..4, 2, 2]
b = GSL::Matrix::Int[5..10, 2, 3]

# Concatenate them horizontally using Matrix::Int#horzcat
a.horzcat(b)

# Concatenate them horizontally using Matrix::Int.horzcat
GSL::Matrix::Int.horzcat(a, b)

# Create a 2x2 Matrix::Int and a 3x2 Matrix::Int
a = GSL::Matrix::Int[1..4, 2, 2]
b = GSL::Matrix::Int[5..10, 3, 2]

# Concatenate them vertically using Matrix::Int#vertcat
a.vertcat(b)

# Concatenate them vertically using Matrix::Int.vertcat
GSL::Matrix::Int.vertcat(a, b)
END
