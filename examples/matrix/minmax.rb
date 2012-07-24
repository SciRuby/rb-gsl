#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create test Matrix
a = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)

# Find maximum
a.max

# Find minimum
a.min

# Find minimum and maximum
a.minmax

# Find index of maximimum
a.max_index

# Find index of minimimum
a.min_index

# Find indices of minimum and maximum
minmax_idx = a.minmax_index

# Use minmax_idx to get minimnum
a[minmax_idx[0]]

# Use minmax_idx to get maximnum
a[minmax_idx[1]]

# Show that #isnull returns 0 for non-null Matrix
a.isnull

# Show that #isnull? returns false for non-null Matrix
a.isnull?

# Set all elements to zero
a.set_zero

# Show that #isnull returns 1 for null Matrix
a.isnull

# Show that #isnull? returns true for null Matrix
a.isnull?
END
