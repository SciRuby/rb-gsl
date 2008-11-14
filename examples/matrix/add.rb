#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create to 3x3 matrices
a = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)
b = GSL::Matrix.alloc([6, 7, 8], [2, 3, 4], [3, 4, 5])

# Add together to create new matrix
a + b

# Add b to a, modifying a
a += b

# Show that a is modified
a

# Subtract b from a, modifying a
a -= b

# Add 2 to all elements of a, creating new matrix
a + 2

# Another way to add 2 to all elements of a, creating new matrix
2 + a

# Add 2 to all elements of a, modifying a
a += 2

# Show that a is modified
a

# Subtract 2 from all elements of a, modifying a
a -= 2

# Show that a is modified
a
END
