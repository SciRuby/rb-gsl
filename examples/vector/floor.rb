#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
#Create test Vector
v = GSL::Vector[1.1, 2.7, 3.5, 5.8, -1.2, -2.8, -3.5]

# Show result of #floor
v.floor

# Show result of #ceil
v.ceil

# Show result of #round
v.round
END
