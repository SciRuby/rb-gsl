#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create test Vector
a = GSL::Vector.alloc(1, 2, 3, 4, 5)

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
a.minmax_index

# Show that #isnull returns 0 for non-null Vector
a.isnull

# Show that #isnull? returns false for non-null Vector
a.isnull?

# Set all elements to zero
a.set_zero

# Show that #isnull returns 1 for null Vector
a.isnull

# Show that #isnull? returns true for null Vector
a.isnull?
END
