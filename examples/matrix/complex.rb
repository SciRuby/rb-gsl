#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create 3x3 Matrix::Complex mz
mz = GSL::Matrix::Complex.alloc(3, 3)

# Set element at row 1, column 2 to 3+5.6i
mz.set(1, 2, GSL::Complex[3, 5.6])

# Get element at row 1, column 2
a = mz.get(1, 2)

# Create Matrix::Complex::View of mz
# starting at row 1, column 1 and
# spanning 2 rows and 2 columns
mzv = mz.submatrix(1, 1, 2, 2)

# Create a Vector::Complex::View for row 1 of mz
row = mz.row(1)

# Create a Vector::Complex::Col::View for column 2 of mz
col = mz.col(2)

# Iterate through rows of mz
mz.each_row do |v|
  p v
end

# Iterate through columns of mz
mz.each_col do |v|
  p v
end
END
