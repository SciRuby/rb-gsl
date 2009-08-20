#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# These examples show some of the ways that Matrix#set or its alias Matrix#[]=
# can be invoked. First, create test matrix m...

m = GSL::Matrix[3,4]

# For a single Array argument, i.e. m.set([row0,row1,...]) or
# m[]=[row0,row1,...], the Array's elements are taken as row contents.  Each
# given row must have exactly the same number of elements as the Matrix has
# columns, but the number of rows given need not match the Matrix's row count.
# Extra given rows are ignored, while Matrix rows beyond those given are not
# affected.

m[] = [[1,2,3,4],[5,6,7,8],[9,8,7,6]]

m

# Note the different return values of Matrix#set and Matrix#[]=.  Matrix#set
# return self (see below), but Matrix[]= returns the value to the right of the
# = sign (see above).  This must be standard Ruby behavior since the underlying
# code returns the same value to Ruby regardless of whether it is invoked as
# #set or #[]=.

m.set([[9,8,7,6],[5,4,3,2],[1,0,1,2]])

# For a single non-Array argument, Matrix#set and Matrix#[] are equivalent to
# Matrix#set_all (other than the difference in the return value of Matrix#[] as
# noted above).

m.set(1.2)             # could also use:  m[] = 1.2

# For two arguments with the first being an Array and the second a non-Array,
# i.e. m.set([i,j], x) or m[[i,j]]=x (note the double square brackets), the
# first two elements of the Array must be Fixnums which specify the row and
# column of the element that will be set to the value of the second (non-Array)
# argument.  This special case exists to allow values returned by
# Matrix#max_index and Matrix#min_index to be used as indexes.

m.indgen!

m[m.max_index] = 100

m[m.min_index] = -100

m

# For three arguments with the first two being Fixnums i and j, this sets
# element (i,j) to the value of the last argument.

m[1,2] = 50; m[-2,-3] = -50; m

# For multiple arguments with the first two being Arrays, i.e.
# m.set(row0,row1,...), this behaves as if the rows were given in a single
# Array (see the first case above).

m.set([1,2,3,4], [5,6,7,8], [9,8,7,6])

# All other forms treat all but the last argument as with Matrix#submatrix and
# set the specified elements based on the last argument, which can be a Matrix
# (or Matrix::View), an Array (of Numerics or Arrays of Numerics), a Range, or
# a Numeric.  Matrix, Array, and Range rvalues must have the same number of
# elements as the specified submatrix.  For a Numeric rvalue, all elements of
# the submatrix are set to that value.
#
# See examples/matrix/view_all.rb for additional examples of how to specify
# submatrices.

m[nil,1] = 0; m

m[1,nil] = 1; m

m[1..2,1..3] = 1..6; m

# Also be careful when setting part of a Matrix from another part of the same
# Matrix.  The GSL method that performs this operation uses memcpy, which does
# not handle overlapping memory regions in a well defined way.

m.indgen!

# This is faster but has problems with overlap
m[1..2,1..2] = m[0..1,0..1]; m

n = GSL::Matrix[3,4].indgen!

# Converting right hand side to Array avoids the problem, but is slower
n[1..2,1..2] = n[0..1,0..1].to_a; n

# See the difference at element [2,2]
n-m
END
