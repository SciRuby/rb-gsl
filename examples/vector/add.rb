#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create two Vectors
a = GSL::Vector.alloc(1, 2, 3, 4, 5)
b = GSL::Vector.alloc(6, 7, 8, 9, 10)

# Show the sum of a and b
a + b

# Add b to a
a += b

# Show a
a

# Subtract b from a
a -= b

# Show a
a

# Create new Vector containing 2 added to all elements of a
a + 2

# Show that a remains unmodified
a

# Add 2 to all elements of a in-place
a += 2

# Show the modified a
a

# Subtract 2 from all elements of a in-place
a -= 2

# Show that a remains unmodified
a

# Show coersion
2 + a

5 - a

# Show that a remains unmodified
a
END
