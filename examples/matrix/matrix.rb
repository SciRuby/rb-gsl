#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create test matrix m
m = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)

# Set rows of m from Arrays
m.set([6, 5, 6], [4, 5, 7], [8, 5, 21])

# Set column 1 of m from GSL::Vector
m.set_col(1, GSL::Vector[12, 3, 55])

# Create transpose of m
m2 = m.transpose

# Swap rows 1 and 2 of m
m.swap_rows(1, 2)

# Create Vector::Col::View for column 0 of m
v = m.col(0)

# Create Vector::View of diagonal of m
m.diagonal

# Create Array containing diagonal elements of m
m.diagonal.to_a

# Create another test matrix m
m = GSL::Matrix.alloc([1, 2, 3], [6, 5, 4], [7, 8, 1])

# Get element at row 1, column 2
m.get(1, 2)

# Perform LU decomposition of m
lu, perm, sign = m.LU_decomp

# Create 5x5 zero matrix m5
m5 = GSL::Matrix.alloc(5, 5)

# Initialize elements of m5
for i in 0...5 do
  for j in 0...5 do
    m5[i, j] = 0.5*(i+0.4)*(j+1.2)
  end
end

# Show m5
m5
END
