#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create two 3x3 test matrices
a = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)
b = GSL::Matrix.alloc([6, 7, 8], [2, 3, 4], [3, 4, 5])

# Multiply elements of Matrix a by 2 (Matrix a remains unmodifided)
a * 2

# Multiply elements of Matrix a by 2 (Matrix a remains unmodifided)
2 * a

# Multiply elements of Matrix a by 2, modifying Matrix a
a *= 2

# Show a
a

# Divide elements of Matrix a by 2, modifying Matrix a
a /= 2

# Show a
a

# Matrix-multiply Matrix a and Matrix b
a*b

# Do element-wise multiplication of a and b
a.mul_elements(b)
END
