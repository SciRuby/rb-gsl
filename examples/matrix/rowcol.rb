#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create 3x3 test Matrix a
a = GSL::Matrix.alloc([1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)

# Create Vector initialized from row 0 of Matrix a
b = a.get_row(0)

# Create Vector::Col for column 2 of Matrix a
c = a.get_col(2)

# Set column 2 of Matrix a from Vector from row 0
a.set_col(2, b)

# Set row 0 of Matrix a from Vector::Col from column 2
a.set_row(0, c)

# Create new Matrix from a with rows 1 and 2 swapped
a.swap_rows(1, 2)

# Show that Matrix a remains unmodified
a

# Swap columns 1 and 2 of Matrix a in-place
a.swap_cols!(1, 2)

# Show that Matrix a is modified
a

# Create new Matrix that is transpose of Matrix a
atrans = a.transpose

# Transpose Matrix a in-place
a.transpose!

# Transpose Matrix a in-place again
a.transpose!
END
