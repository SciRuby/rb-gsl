#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# These examples show all(?) the ways that Vector#get or its alias Vector#[]
# can be invoked.  For a single Fixnum argument, a single element is returned.
# For a single Array or GSL::Permutation argument, a new GSL::Vector containing
# the specified elements, in the specified order, is returned.  For all other
# cases, Vector#get is essentially an alias for Vector#subvector.  See
# examples/vector/view_all.rb for more examples.

# Create test vector v
v = GSL::Vector.indgen(9)

# Vector#[] with zero args returns a Vector::View of entire Vector
v[]

# Vector#[] with one Fixnum argument, i, returns the i'th element if i is
# positive or the (size+i)'th element if i is negative.
v[3]
v[-3]

# Vector#[] with single Array argument.  Notice the inner pair of brackets!
v[[1,4,-9]]

# Vector#[] with a single GSL::Permutation argument.
p = GSL::Permutation.calloc(4).reverse
v[p]

# Vector#[] with one Range argument returns a Vector::View of the specified
# elements.  If the begin value is greater than the end value, the View will
# have the elements in reverse order.  If begin and/or end value is negative,
# the value is taken to be "size-n".
v[1..4]
v[4..1]
v[1...4]
v[4...1]

v[4..-2]
v[-2..4]
v[4...-2]
v[-2...4]

v[-4..8]
v[8..-4]
v[-4...8]
v[8...-4]

v[-5..-2]
v[-2..-5]
v[-5...-2]
v[-2...-5]

# Vector#[] with a Range argument and a Fixnum argument is like a single
# Range argument, but with a stride given by the Fixnum argument.
v[1..7, 3]
v[7..1, 3]
v[1...7, 3]
v[7...1, 3]

# Vector#[] with two Fixnum arguments is offset, length.  If offset is
# negative, it means size+offset.  If length is negative, it means step is -1.
v[2, 4]
v[4, 2]
v[-4, 2]
v[-2, -4]

# Vector#[] with three Fixnum arguments is offset, stride, length.  If offset
# is negative, it means size+offset.  If length is negative, the sign of both
# stride and length is inverted.
v[1, 2, 3]
v[1, -2, -3]
v[-1, -2, 3]
v[-1, 2, -3]
END
