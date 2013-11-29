#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create 10,000 random numbers from Gaussian distribution
N = 10000
r = GSL::Rng.alloc
v = r.gaussian(1.0, N)    # Generate N random numbers

# Bin them into 50 bins spanning -4 to +4
h = v.histogram(50, [-4, 4])

# Plot them using graph
h.graph("-T X -C -g 3")
END
