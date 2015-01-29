#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Setup constants
N = 1000
DECIMATE1 = 10
DECIMATE2 = 100

# Setup random number generator
r = GSL::Rng.alloc

# Create Vector of x values
x0 = GSL::Vector.linspace(0, 20, N)

# Data: Bessel function + noise
y0 = GSL::Sf::bessel_J0(x0) + GSL::Ran::gaussian(r, 0.1, N)

# Decimate y0 by DECIMATE1
y1 = y0.decimate(DECIMATE1)

# Decimate y0 by DECIMATE2
y2 = y0.decimate(DECIMATE2)

# Create Vectors of decimated x values
x1 = GSL::Vector.linspace(0, 20, N/DECIMATE1)
x2 = GSL::Vector.linspace(0, 20, N/DECIMATE2)

# y1 and y2 are shifted vertically for visual purpose
GSL::graph([x0, y0], [x1, y1-1], [x2, y2-2])
END
