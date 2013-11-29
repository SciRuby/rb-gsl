#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create matrix with diagonal given by Range
GSL::Matrix.diagonal(1..3)

# Create matrix with diagonal given by Array
GSL::Matrix.diagonal([1, 2, 3])

# Create matrix with diagonal given by GSL::Vector
GSL::Matrix.diagonal(GSL::Vector.indgen(3,1))

# Create matrix with diagonal given by individual elements
GSL::Matrix.diagonal(1, 2, 3)
END
