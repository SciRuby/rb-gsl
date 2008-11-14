#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create three Vectors: x, y, z
x = GSL::Vector::Int[1, 2, 3]
y = GSL::Vector::Int[1, 2, 5]
z = GSL::Vector::Int[0, 2, 9]

# Test element-wise "==" method
x.eq(y)

# Test element-wise "!=" method
x.ne(y)

# Test element-wise ">=" method
x.ge(y)

# Test element-wise "<" method
x.lt(z)
END
