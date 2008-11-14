#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Initialize random number generator with fixed seed for repeatable sequence
r = GSL::Rng.alloc(GSL::Rng::MT19937, ?r-?b+?g-?s-?l)

# Create 3x3 matrix initialized with numbers from uniform distribution
u = GSL::Matrix.rand(3, 3, r)

# Create 3x3 matrix initialized with numbers from normal distribution
n = GSL::Matrix.randn(3, 3, r)
END
