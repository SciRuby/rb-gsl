#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Crete test matrix m
m = GSL::Matrix[[1, 2, 3, 4, 5, 6, 7, 8 ,0], 3, 3]

# Calculate determinant of m
m.det

# Calculate trace of m (sum of diagonal elements)
m.trace

# Convert to Matrix::Complex mz
mz = m.to_complex

# Calulate determinant of mz
mz.det

# Calculate trace of mz (sum of diagonal elements)
mz.trace
END
