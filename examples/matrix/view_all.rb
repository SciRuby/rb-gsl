#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# These examples show all(?) the ways that a Matrix::View can be created using
# Matrix#submatrix or its alias Matrix#view.  Note that Matrix#get or, more
# commonly, its alias Matrix#[] can also be used to create a Matrix::View.  See
# examples/matrix/get_all.rb for more examples.

# Create 4x4 test matrix m
m = GSL::Matrix.indgen(4, 4)

# Matrix#view with zero args returns a Matrix::View of entire Matrix
m.view

# Matrix#view with one Fixnum argument, i, treats the Matrix as a Vector and
# returns a Matrix::View of the i'th element if i is positive or the
# (i+size1*size2)'th element if i is negative.
m.view(3)
m.view(-3)

# When Matrix#view is called with two arguments, the first specifies which
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

# Matrix#view(nil, nil) returns a Matrix::View of entire Matrix
m.view(nil, nil)

# Matrix#view(Fixnum, Fixnum) returns a Matrix::View of the single element at
# the specified row and column.
m.view(2, 3)
m.view(-1, -3)

# Matrix#view(Range, Range) returns a Matrix::View of the rows and columns
# specified by the two Ranges.
m.view(0...2, 1..-2)
m.view(-3..-1, -4...4)

# Matrix#view(Fixnum, nil) returns a Vector::View of the entire row specified
# by the Fixnum argument.  A negative value is treated as counting backwards
# from the end of the corresponding dimension.  NOTE: GSL 1.11 (and maybe
# earlier versions) has a bug that prevents this from working if the Matrix has
# more columns than rows!!!
m.view(1, nil)
m.view(-2, nil)

# Matrix#view(nil, Fixnum) returns a Vector::Col::View of the entire column
# specified by the Fixnum argument.  A negative value is treated as counting
# backwards from the end of the corresponding dimension.  NOTE: GSL 1.11 (and
# maybe earlier versions) has a bug that prevents this from working if the
# Matrix has more rows than columns!!!
m.view(nil, 1)
m.view(nil, -2)

# Matrix#view(Range,nil) returns a Matrix::View of all columns and the rows
# specified by the Range argument.  Note that negative begin and/or end values
# are treated as counting backwards from the end of corresponding dimension.
m.view(1...3, nil)
m.view(0..-2, nil)
m.view(-2..3, nil)
m.view(-2..-1, nil)

# Matrix#view(nil, Range) returns a Matrix::View of all rows and the columns
# specified by the Range argument.  Note that negative begin and/or end values
# are treated as counting backwards from the end of corresponding dimension.
m.view(nil, 1...3)
m.view(nil, 0..-2)
m.view(nil, -2..3)
m.view(nil, -2..-1)

# Matrix#view(Range, Fixnum) returns a Vector::Col::View of the rows specified
# by the Range argument of the column specified by the Fixnum argument.  A
# negative value is treated as counting backwards from the end of the
# corresponding dimension.  NOTE: GSL 1.11 (and maybe earlier versions) has a
# bug that prevents this from working if the Matrix has more rows than
# columns!!!
m.view(1...3, 0)
m.view(0..-2, 1)
m.view(-2..3, -2)
m.view(-2..-1, -1)

# Matrix#view(Fixnum, Range) returns a Vector::View of the columns specified
# by the Range argument of the row specified by the Fixnum argument.  A
# negative value is treated as counting backwards from the end of the
# corresponding dimension.  NOTE: GSL 1.11 (and maybe earlier versions) has a
# bug that prevents this from working if the Matrix has more rows than
# columns!!!
m.view(0, 1...3)
m.view(1, 0..-2)
m.view(-2, -2..3)
m.view(-1, -2..-1)

# When Matrix#view is called with three arguments, the first or last argument
# must be nil or a Range and the other two arguments must be Fixnums.  The two
# Fixnums indicate a span whose offset is given by the first Fixnum and whose
# length is given by the second Fixnum.  If they are the first two arguments,
# they indicate which rows the returned view will cover.  If they are the last
# two arguments, they indicate which columns the returned view will cover.  The
# nil or Range argument indicate what portion of the other dimension will be
# covered by the returned view (nil means all rows or columns).

# Matrix#view(nil, Fixnum, Fixnum) returns a Matrix::View covering all rows of
# the column span specified by the two Fixnums.
m.view(nil, 1, 2)
m.view(nil, -2, 2)

# Matrix#view(Fixnum, Fixnum, nil) returns a Matrix::View covering all columns
# of the row span specified by the two Fixnums.
m.view(nil, 1, 2)
m.view(nil, -2, 2)

# Matrix#view(Range, Fixnum, Fixnum) returns a Matrix::View covering Range rows
# of the column span specified by the two Fixnums.
m.view(0...2, -3, 2)
m.view(1..-2, 1, 2)
m.view(-3..-1, 3, 1)
m.view(-4...4, -4, 2)

# Matrix#view(Fixnum, Fixnum, Range) returns a Matrix::View covering Range
# columns of the row span specified by the two Fixnums.
m.view(-3, 2, 0...2)
m.view(1, 2, 1..-2)
m.view(3, 1, -3..-1)
m.view(-4, 2, -4...4)

# When Matrix#view is called with four arguments, all four arguments must be
# Fixnums.  The first two Fixnums specify the Matrix element that will be the
# upper left corner of the view (negative values are treated as counting
# backwrds from the end of the corresponding dimension).  The last two Fixnums
# specify the number of rows and columns that the view will have.
m.view(0, 1, 2, 3)
m.view(-2, -3, 2, 1)
END
