#!/usr/bin/env ruby
# This script is written by Cameron McBride.
# Thank you Cameron!

require 'gsl'

v = GSL::Vector[0..9]
#v = GSL::Vector[0..9].block
#v = GSL::Vector::Int[0..9]
#v = GSL::Vector::Int[0..9].block

puts "Vector:"
p v

print "\n" + '='*50 + "\n"

puts "mask = Vector > 2"
p mask = (v > 2)

print "\n" + '='*50 + "\n"

puts "mask where (and with block { true } / { false })"
p mask.where
p mask.where { true }
p mask.where { false}

print "\n" + '='*50 + "\n"

puts "mask where2 (and with block { true } / { false } / { rand > 0.5 })"
p mask.where2
p mask.where2 { true  }
p mask.where2 { false }
p mask.where2 { rand > 0.5 }


