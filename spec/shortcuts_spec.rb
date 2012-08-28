# = NMatrix
#
# A linear algebra library for scientific computation in Ruby.
# NMatrix is part of SciRuby.
#
# NMatrix was originally inspired by and derived from NArray, by
# Masahiro Tanaka: http://narray.rubyforge.org
#
# == Copyright Information
#
# SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
# NMatrix is Copyright (c) 2012, Ruby Science Foundation
#
# Please see LICENSE.txt for additional copyright notices.
#
# == Contributing
#
# By contributing source code to SciRuby, you agree to be bound by
# our Contributor Agreement:
#
# * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
#
# == shortcuts_spec.rb
#
# Basic tests for the shortcuts used in NMatrix and NVector.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")

describe "NMatrix shortcuts" do
    
  it "creates a matrix of zeros" do
  end
  
  it "creates a matrix of ones" do
  end
  
  it "creates an identity matrix" do
  end
  
  it "creates a matrix of random numbers" do
    # How to test this method?
  end
  
  it "creates a matrix of integers, sequentially" do
  end
end

describe "NVector shortcuts" do
      
  it "creates a vector of zeros" do
  end
  
  it "creates a vector of ones" do
  end
  
  it "creates a vector of random numbers" do
    # How to test this method?
  end
  
  it "creates a vector of integers, sequentially" do
  end
  
  it "creates a vector with n values equally spaced between a and b" do
  end
end

describe "Inline constructor" do
  
  it "creates a NMatrix with the given values" do
  end
end