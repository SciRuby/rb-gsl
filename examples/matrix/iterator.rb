#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create test matrix
m = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8 ,9], 3, 3)

# Iterate through columns
m.each_col do |v|
  p v
end
END
