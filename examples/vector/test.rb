#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create one Vector and two Vector::Ints
x = GSL::Vector[1, 2, 3]
y = GSL::Vector::Int[1, 0, 5]
a = GSL::Vector::Int[0, 0, 0]

# Call #any? on them
x.any?
y.any?
a.any?

# Call #all? on them
x.all?
y.all?
a.all?

# Call #none? on them
x.none?
y.none?
a.none?

# Call x.any? with blocks
x.any? { |val| val > 5 }
x.any? { |val| val > 2 }

# Call x.all? with blocks
x.all? { |val| val >= 1 }
x.all? { |val| val >= 2 }

# Call x.none? with blocks
x.none? { |val| val == 1 }
x.none? { |val| val == 5 }
END
