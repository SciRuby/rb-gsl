#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create 3x3 Hilbert matrix m
m = GSL::Matrix.hilbert(3)

# Compute inverse of m
invm = m.inv

# Create inverse of 3x3 Hilbert matrix directly
invm2 = GSL::Matrix.invhilbert(3)

# Show that both inverse matrices are inverses of m
m*invm
m*invm2

# Show that the two inverse matrices are equal
# to absolute accuracy eps = 1e-10
invm == invm2

# Show that they may not be exactly equal
invm - invm2
END
