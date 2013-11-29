#!/usr/bin/env ruby
# Turn on warnings
$-w = true

# This script is written by Cameron McBride.
# Thank you Cameron!

#gem 'gsl', '=1.10.3'
require 'irb/xmp'
require 'gsl'
puts "Using rb-gsl #{GSL::RB_GSL_VERSION}"

# Seed random number generator for repeatable results
srand(?r-?b+?g-?s-?l)

# The creation of the test Vector below can be replaced with any of these lines
# to show the behavior with GSL::Block, GSL::Vector::Int, and
# GSL::Vector::Int::Block.  Note that the block objects should NOT be created
# using "v = GSL::Vector[0..9].block" because the Vector object will be subject
# to garbage collection which will free the memory that the Block object points
# to.  See the BUGS file for more information.
#
# vv = GSL::Vector[0..9]; v = vv.block
# v = GSL::Vector::Int[0..9]
# vv = GSL::Vector::Int[0..9]; v = vv.block

# Apparently, IRB::Frame has a bug that prevents the defaults from working, so
# an XMP instance must be created explicitly this way instead of using the
# otherwise convenient xmp method.
XMP.new(IRB::Frame.top(-1)).puts <<END
# Create test Vector v
v = GSL::Vector[0..9]

# Create mask
mask = (v > 2)

# Show mask.where with no block
mask.where

# Show mask.where with "true" block
mask.where { true }

# Show mask.where with "false" block
# NB: Returns nil, not empty object!
mask.where { false }

# Show mask.where with no block
mask.where2

# Show mask.where with "true" block
mask.where2 { true  }

# Show mask.where with "false" block
mask.where2 { false }

# Show mask.where with random block
mask.where2 { rand > 0.5 }
END
