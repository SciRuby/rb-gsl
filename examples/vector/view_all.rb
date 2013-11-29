#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# These examples show all(?) the ways that a Vector::View can be created using
# Vector#subvector or its alias Vector#view.  Note that Vector#get or, more
# commonly, its alias Vector#[] can also be used to create a Vector::View.  See
# examples/vector/get_all.rb for more examples.

# Create test vector v
v = GSL::Vector.indgen(9)

# Vector#view with zero args returns a Vector::View of entire Vector
v.view

# Vector#view with one Fixnum argument, i, returns a Vector::View of the first
# i'th elements if i is positive or the last i'th elements if i is negative.
v.view(3)
v.view(-3)

# Vector#view with one Range argument returns a Vector::View of the specified
# elements.  If the begin value is greater than the end value, the View will
# have the elements in reverse order.  If begin and/or end value is negative,
# the value is taken to be "size-n".
v.view(1..4)
v.view(4..1)
v.view(1...4)
v.view(4...1)

v.view(4..-2)
v.view(-2..4)
v.view(4...-2)
v.view(-2...4)

v.view(-4..8)
v.view(8..-4)
v.view(-4...8)
v.view(8...-4)

v.view(-5..-2)
v.view(-2..-5)
v.view(-5...-2)
v.view(-2...-5)

# Vector#view with a Range argument and a Fixnum argument is like a single
# Range argument, but with a stride given by the Fixnum argument.
v.view(1..7, 3)
v.view(7..1, 3)
v.view(1...7, 3)
v.view(7...1, 3)

# Vector#view with two Fixnum arguments is offset, length.  If offset is
# negative, it means size+offset.  If length is negative, it means step is -1.
v.view(2, 4)
v.view(4, 2)
v.view(-4, 2)
v.view(-2, -4)

# Vector#view with three Fixnum arguments is offset, stride, length.  If offset
# is negative, it means size+offset.  If length is negative, the sign of both
# stride and length is inverted.
v.view(1, 2, 3)
v.view(1, -2, -3)
v.view(-1, -2, 3)
v.view(-1, 2, -3)
END
