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
# == lapack_spec.rb
#
# Tests for properly exposed LAPACK functions.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")

describe NMatrix::LAPACK do
  [:rational32, :rational64, :rational128, :float32, :float64, :complex64, :complex128].each do |dtype|
    context dtype do
      it "exposes clapack getrf" do
        a = NMatrix.new(:dense, 3, [4,9,2,3,5,7,8,1,6], dtype)
        NMatrix::LAPACK::clapack_getrf(:row, 3, 3, a, 3)
        a[0,0].should == 8
        a[0,1].should == 1
        a[0,2].should == 6
        a[1,0].should == 1.quo(2)
        a[1,1].should == 17.quo(2)
        a[1,2].should == -1
        a[2,0].should == 3.quo(8)
        # FIXME: these are rounded, == won't work
        #a[2,1].should == 0.544118
        #a[2,2].should == 5.294118
      end

      it "exposes clapack getrs" do
        a     = NMatrix.new(:dense, 3, [-2,4,-3,3,-2,1,0,-4,3], dtype)
        ipiv  = NMatrix::LAPACK::clapack_getrf(:row, 3, 3, a, 3)
        b     = NVector.new(3, [-1, 17, -9], dtype)

        NMatrix::LAPACK::clapack_getrs(:row, false, 3, 1, a, 3, ipiv, b, 1)

        b[0].should == 5
        b[1].should == -15.quo(2)
        b[2].should == -13
      end

      it "exposes cblas trsm, with B as a matrix"
    end
  end
end