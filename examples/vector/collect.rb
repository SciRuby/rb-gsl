#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create Vector v
v = GSL::Vector::Int[0..5]

# Create new Vector whose elements are squares of v's elements
v.collect { |a| a*a }

# Show that v us unmodified
v

# Square elements of v in-place
v.collect! { |a| a*a }

# Show that v is modified
v
END
