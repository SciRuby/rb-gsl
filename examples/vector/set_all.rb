#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# These examples show some of the ways that Vector#set or its alias Vector#[]=
# can be invoked.  For a single argument, this is equivalent to Vector#set_all.
# For two arguments with the first a Fixnum i, this sets the i'th element (or
# the (size-i)'th element to the value of the second argument.  All other forms
# treat all but the last argument as with Vector#subvector and set the
# specified elements based on the last argument, which can be a Vector (or
# Vector::View), an Array, a Range, or a Numeric.  Vector, Array, and Range
# rvalues must have the same number of elements as the specified subvector.
# For a Numeric rvalue, all elements of the subvector are set to that value.
#
# Note the different return values of Vector#set and Vector#[]=.  Vector#set
# return self, but Vector[]= return the value to the right of the = sign.  This
# must be standard Ruby behavior since the underlying code returns the same
# value to Ruby regardless of whether it is invoked as #set or #[]=.
#
# Also be careful is setting part of a Vector from another part of the same
# vector.  The GSL method that performs this operation uses memcpy, which does
# not handle overlapping memory regions in a well defined way.  See the last
# two examples.
#
# See examples/vector/view_all.rb for additional examples of how to specify
# subvectors.

# Create test vector v
v = GSL::Vector.indgen(9)

# Vector#set and Vector#[]= with one arg sets all elements
v.set(1.2)
v[] = 3.4
v

# Vector#[i]= Numeric sets the i'th element if i is
# positive or the (size+i)'th element if i is negative.
v[3] = 5.6
v[-8] = 7.8
v

# Specifying subvector using Range with various rvalue types
v[1..4] = GSL::Vector[2, 3, 5, 7] # rvalue is Vector
v

v[1..4] = [11, 13, 17, 19] # rvalue is Array
v

v[1..4] = 24..27 # rvalue is Range
v

v[1..4] = 1.0 # rvalue is Numeric
v

# Specifying subvector using Range and stride with various rvalue types
v[0..4, 2] = GSL::Vector[2, 3, 5] # rvalue is Vector
v

v[0..4, 2] = [7, 11, 13] # rvalue is Array
v

v[0..4, 2] = 8..10 # rvalue is Range
v

v[0..4, 2] = 1.0 # rvalue is Numeric
v

# Specifying subvector using two Fixnums (offset, length) with various rvalue
# types
v[2, 4] = GSL::Vector[2, 3, 5, 7] # rvalue is Vector
v

v[2, 4] = [11, 13, 17, 19] # rvalue is Array
v

v[2, 4] = 24..27 # rvalue is Range
v

v[2, 4] = 1.0 # rvalue is Numeric
v

# Specifying subvector using three Fixnum arguments (offset, stride, length)
# with various rvalue types
v[1, 2, 3] = GSL::Vector[2, 3, 5] # rvalue is Vector
v

v[1, 2, 3] = [7, 11, 13] # rvalue is Array
v

v[1, 2, 3] = 8..10 # rvalue is Range
v

v[1, 2, 3] = 1.0 # rvalue is Numeric
v

# Copying part of a Vector to another part of the same Vector can potentially
# be problematic if the regions overlap.
v.indgen!
v[0, 3] = v[2, 3]
v

v.indgen!
v[2, 3] = v[0, 3]
v

# But it's OK if the regions do not overlap
v.indgen!
v[0, 3] = v[3, 3]
v

v.indgen!
v[3, 3] = v[0, 3]
v
END
