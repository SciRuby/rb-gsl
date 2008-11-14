#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create test Vector::Int
z = GSL::Vector::Int[1, 0, 0, 5, 2, 0, 9]


# #where without a block indices of non-zero elements.
z.where

# #where with a block returns indices of elements for which block returns true.
z.where {|e| e >= 2}

# #where2 without a block returns ["indices of non-zero elemens", "indices of
# zero elements"]
z.where2

# #where2 with a block returns ["indices of elements for which block returned
# true", "indices of elements for which block returned false"]
z.where2 {|e| e >= 2}
END
