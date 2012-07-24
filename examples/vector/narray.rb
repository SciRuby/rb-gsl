#!/usr/bin/env ruby
# Turn on warnings
$-w = true

require 'irb/xmp'
require 'gsl'

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create NMatrix m
m = NMatrix[[0, 1.2, 1],[1.5, 0, 2]]

# Convert NMatrix m to Vector gv
gv = m.to_gv

# Convert Vector gv to NArray m2
m2 = gv.to_na
m2.class

# Create GSL::Vector v
v = GSL::Vector.alloc(1..4)

# Convert v to NArray na
na = v.to_na
na.class

# Convert na back to Vector
v2 = na.to_gv

# Create new Vector copy of na
v3 = GSL::Vector.alloc(na)

# Create Vector::View of na
v4 = na.to_gv_view

# Set element of Vector::View
v4[2] = 123

# Show that na was modified
na

# Show that v3 was not modified
v3
END
