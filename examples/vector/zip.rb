#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create three test vectors of varying lengths
a = GSL::Vector[0..4]
b = GSL::Vector[2, 3, 4]
c = GSL::Vector[5, 7, 4, 8, 9, 2]

# Zip them together using Vector#zip
a.zip(b, c)

# Zip them together using GSL::Vector.zip
GSL::Vector.zip(a, b, c)

# Convert test Vectors to Vector::Complex
aa = a.to_complex
bb = b.to_complex
cc = c.to_complex

# Zip Vector::Complex objects together using Vector::Complex#zip
aa.zip(bb, cc)


# Zip Vector::Complex objects together using GSL::Vector::Complex.zip
GSL::Vector::Complex.zip(aa, bb, cc)
END
