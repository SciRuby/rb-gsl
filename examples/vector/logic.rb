#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create two test Vectors
a = GSL::Vector[1, 0, 3, 0]
b = GSL::Vector[3, 4, 0, 0]

# Show result of (a AND b)
a.and(b)
b.and(a)

# Show result of (a OR b)
a.or(b)
b.or(a)

# Show result of (a XOR b)
a.xor(b)
b.xor(a)

# Show result of (NOT a)
p a.not

# Show result of (NOT b)
p b.not
END
