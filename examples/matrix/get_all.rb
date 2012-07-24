#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# These examples show all(?) the ways that Matrix#get or its alias Matrix#[]
# can be invoked.  For one or two Fixnum arguments or a single two-element
# Array argument, a single element is returned.  For all other cases,
# Matrix#get is essentially an alias for Matrix#submatrix.  See
# examples/matrix/view_all.rb for more examples.

# Create 4x4 test matrix m
m = GSL::Matrix.indgen(4, 4)

# Matrix#[] with zero args returns a Matrix::View of entire Matrix
m[]

# Matrix#[] with one Fixnum argument, i, treats the Matrix as a Vector and
# returns a Matrix::View of the i'th element if i is positive or the
# (i+size1*size2)'th element if i is negative.
m[3]
m[-3]

# Matrix#[Fixnum, Fixnum] returns a Matrix::View of the single element at
# the specified row and column.
m[2, 3]
m[-1, -3]

# Matrix#[[Fixnum, Fixnum]] (note the double square brackets) is treated the
# same as Matrix#[Fixnum, Fixnum].  This special case exists to allow values
# returned by Matrix#max_index and Matrix#min_index to be used as indexes.
m[[2, 3]]
m[[-1, -3]]
m.max_index
m[m.max_index]
m.min_index
m[m.min_index]

# When Matrix#[] is called with two arguments, the first specifies which
# row(s) the view will cover and the second specifies which column(s) the view
# will cover.  The arguments may be nil (indicating all rows or columns),
# a Fixnum (indicating a single row or column), or a Range (indicating a range
# of rows or columns).  The return type is Matrix::View unless exactly one
# argument is a Fixnum in which case a Vector::View or Vector::Col::View is
# returned.
#
# NOTE: GSL 1.11 (and maybe earlier versions) has a bug that can prevent the
# exactly-one-Fixnum case from working properly if the Matrix does not have an
# equal number of rows and columns!!!

# Matrix#[nil, nil] returns a Matrix::View of entire Matrix
m[nil, nil]

# Matrix#[Range, Range] returns a Matrix::View of the rows and columns
# specified by the two Ranges.
m[0...2, 1..-2]
m[-3..-1, -4...4]

# Matrix#[Fixnum, nil] returns a Vector::View of the entire row specified
# by the Fixnum argument.  A negative value is treated as counting backwards
# from the end of the corresponding dimension.  NOTE: GSL 1.11 (and maybe
# earlier versions) has a bug that prevents this from working if the Matrix has
# more columns than rows!!!
m[1, nil]
m[-2, nil]

# Matrix#[nil, Fixnum] returns a Vector::Col::View of the entire column
# specified by the Fixnum argument.  A negative value is treated as counting
# backwards from the end of the corresponding dimension.  NOTE: GSL 1.11 (and
# maybe earlier versions) has a bug that prevents this from working if the
# Matrix has more rows than columns!!!
m[nil, 1]
m[nil, -2]

# Matrix#[Range,nil] returns a Matrix::View of all columns and the rows
# specified by the Range argument.  Note that negative begin and/or end values
# are treated as counting backwards from the end of corresponding dimension.
m[1...3, nil]
m[0..-2, nil]
m[-2..3, nil]
m[-2..-1, nil]

# Matrix#[nil, Range] returns a Matrix::View of all rows and the columns
# specified by the Range argument.  Note that negative begin and/or end values
# are treated as counting backwards from the end of corresponding dimension.
m[nil, 1...3]
m[nil, 0..-2]
m[nil, -2..3]
m[nil, -2..-1]

# Matrix#[Range, Fixnum] returns a Vector::Col::View of the rows specified
# by the Range argument of the column specified by the Fixnum argument.  A
# negative value is treated as counting backwards from the end of the
# corresponding dimension.  NOTE: GSL 1.11 (and maybe earlier versions) has a
# bug that prevents this from working if the Matrix has more rows than
# columns!!!
m[1...3, 0]
m[0..-2, 1]
m[-2..3, -2]
m[-2..-1, -1]

# Matrix#[Fixnum, Range] returns a Vector::View of the columns specified
# by the Range argument of the row specified by the Fixnum argument.  A
# negative value is treated as counting backwards from the end of the
# corresponding dimension.  NOTE: GSL 1.11 (and maybe earlier versions) has a
# bug that prevents this from working if the Matrix has more rows than
# columns!!!
m[0, 1...3]
m[1, 0..-2]
m[-2, -2..3]
m[-1, -2..-1]

# When Matrix#[] is called with three arguments, the first or last argument
# must be nil or a Range and the other two arguments must be Fixnums.  The two
# Fixnums indicate a span whose offset is given by the first Fixnum and whose
# length is given by the second Fixnum.  If they are the first two arguments,
# they indicate which rows the returned view will cover.  If they are the last
# two arguments, they indicate which columns the returned view will cover.  The
# nil or Range argument indicate what portion of the other dimension will be
# covered by the returned view (nil means all rows or columns).

# Matrix#[nil, Fixnum, Fixnum] returns a Matrix::View covering all rows of
# the column span specified by the two Fixnums.
m[nil, 1, 2]
m[nil, -2, 2]

# Matrix#[Fixnum, Fixnum, nil] returns a Matrix::View covering all columns
# of the row span specified by the two Fixnums.
m[nil, 1, 2]
m[nil, -2, 2]

# Matrix#[Range, Fixnum, Fixnum] returns a Matrix::View covering Range rows
# of the column span specified by the two Fixnums.
m[0...2, -3, 2]
m[1..-2, 1, 2]
m[-3..-1, 3, 1]
m[-4...4, -4, 2]

# Matrix#[Fixnum, Fixnum, Range] returns a Matrix::View covering Range
# columns of the row span specified by the two Fixnums.
m[-3, 2, 0...2]
m[1, 2, 1..-2]
m[3, 1, -3..-1]
m[-4, 2, -4...4]

# When Matrix#[] is called with four arguments, all four arguments must be
# Fixnums.  The first two Fixnums specify the Matrix element that will be the
# upper left corner of the view (negative values are treated as counting
# backwrds from the end of the corresponding dimension).  The last two Fixnums
# specify the number of rows and columns that the view will have.
m[0, 1, 2, 3]
m[-2, -3, 2, 1]
END
