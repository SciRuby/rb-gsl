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
# == blas_spec.rb
#
# Tests for properly exposed BLAS functions.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")

describe NMatrix::BLAS do
  [:rational32, :rational64, :rational128, :float32, :float64, :complex64, :complex128].each do |dtype|
    context dtype do
      it "exposes cblas trsm" do
        a     = NMatrix.new(:dense, 3, [4,-1.quo(2), -3.quo(4), -2, 2, -1.quo(4), -4, -2, -1.quo(2)], dtype)
        b     = NVector.new(3, [-1, 17, -9], dtype)
        NMatrix::BLAS::cblas_trsm(:row, :right, :lower, :transpose, :nonunit, 1, 3, 1.0, a, 3, b, 3)

        # These test results all come from actually running a matrix through BLAS. We use them to ensure that NMatrix's
        # version of these functions (for rationals) give similar results.

        b[0].should == -1.quo(4)
        b[1].should == 33.quo(4)
        b[2].should == -13

        NMatrix::BLAS::cblas_trsm(:row, :right, :upper, :transpose, :unit, 1, 3, 1.0, a, 3, b, 3)

        b[0].should == -15.quo(2)
        b[1].should == 5
        b[2].should == -13
      end
    end
  end
end