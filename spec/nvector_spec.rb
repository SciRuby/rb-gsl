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
# == nvector_spec.rb
#
# Basic tests for NVector.
#

require "./lib/nmatrix"

describe NVector do
  it "correctly initializes" do
    v = NVector.new 5, :float64
    v.shape[0].should == 5
    v.shape[1].should == 1
  end

  it "permits setting and getting contents" do
    v = NVector.new 5, :float64
    v[0] = 1.555
    v[0].should == 1.555
  end

end