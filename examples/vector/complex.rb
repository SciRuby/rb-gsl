#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create a Vector::Complex of 5 elements, all 0+0i
v = GSL::Vector::Complex.alloc(5)

# Set element 2 to 3+4i
v[2] = [3, 4]

# Show vector
v

# Show element 2
v[2]

# Use #map! to modify each element of vector in-place
i = 0
v.map! do |elm|
  i += 1
  elm += i
end

# Show vector
v

# Show element 3
v[3]

# Set all elements to 2+4.7i
v.set_all([2, 4.7])

# Create subvector starting at element 1 and spanning 3 elements
v2 = v.subvector(1, 3)

# Show size of subvector
v2.size

# Get a Vector::View of the real components of Vector::Complex v
p v.real

# Convert v to an Array
v.to_a
END
