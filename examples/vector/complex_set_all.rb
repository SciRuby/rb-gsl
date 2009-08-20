#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# These examples show some of the ways that Vector::Complex#set or its alias
# Vector::Complex#[]= can be invoked.  For a single argument, this is
# equivalent to Vector::Complex#set_all.  For two arguments with the first
# being a Fixnum i, this sets the i'th element (or the (size-i)'th element to
# the complex value derived from the second argument.  The second argument may
# be a Fixnum, Bignum, Float or two element Array.  For the Array case, the
# first element is the real component and the second element is the imaginary
# component.  A nil component leave that component unchanged.  All other forms
# treat all but the last argument as with Vector::Complex#subvector and set the
# specified elements based on the last argument, which can be a Vector::Complex
# (or Vector::Complex::View), Array, Range, Fixnum, Bignum, or Float.
# Vector::Complex, Array, and Range rvalues must have the same number of
# elements as the specified subvector.  For a Fixnum, Bignum, or Float rvalue,
# all elements of the subvector are set to that value.
#
# Note the different return values of Vector::Complex#set and
# Vector::Complex#[]=.  Vector::Complex#set return self, but Vector::Complex[]=
# return the value to the right of the = sign.  This must be standard Ruby
# behavior since the underlying code returns the same value to Ruby regardless
# of whether it is invoked as #set or #[]=.
#
# Also be careful is setting part of a Vector::Complex from another part of the
# same vector.  The GSL method that performs this operation uses memcpy, which
# does not handle overlapping memory regions in a well defined way.  See the
# last two examples.
#
# See examples/vector/complex_view_all.rb for additional examples of how to
# specify subvectors.

# Create test vector v
v = GSL::Vector::Complex.indgen(9)

# Vector::Complex#set and Vector::Complex#[]= with one arg sets all elements
v.set(1.2)
v[] = 3.4
v

# Vector::Complex#[i]= Numeric sets the i'th element if i is
# positive or the (size+i)'th element if i is negative.
v[3] = 5.6
v[-8] = 7.8
v
v[-8] = [nil, 1.0] # Set imaginary component only
v[-8]
v[-8] = [nil, 0.0] # Set imaginary component only
v[-8]

# Specifying subvector using Range with various rvalue types
v[1..4] = GSL::Vector::Complex[[2,3],[5,7],[11,13],[17,19]]
v

v[1..4] = [11, 13, 17, 19] # rvalue is Array
v

v[1..4] = 24..27 # rvalue is Range
v

v[1..4] = 1.0 # rvalue is Float
v

# Specifying subvector using Range and stride with various rvalue types
v[0..4, 2] = GSL::Vector::Complex[[2,3],[5,7],[11,13]]
v

v[0..4, 2] = [7, 11, 13] # rvalue is Array
v

v[0..4, 2] = 8..10 # rvalue is Range
v

v[0..4, 2] = 1.0 # rvalue is Float
v

# Specifying subvector using two Fixnums (offset, length) with various rvalue
# types
v[2, 4] = GSL::Vector::Complex[[2,3],[5,7],[11,13],[17,19]]
v

v[2, 4] = [11, 13, 17, 19] # rvalue is Array
v

v[2, 4] = 24..27 # rvalue is Range
v

v[2, 4] = 1.0 # rvalue is Float
v

# Specifying subvector using three Fixnum arguments (offset, stride, length)
# with various rvalue types
v[1, 2, 3] = GSL::Vector::Complex[[2,3],[5,7],[11,13]]
v

v[1, 2, 3] = [7, 11, 13] # rvalue is Array
v

v[1, 2, 3] = 8..10 # rvalue is Range
v

v[1, 2, 3] = 1.0 # rvalue is Float
v

# Copying part of a Vector::Complex to another part of the same Vector::Complex can potentially
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
