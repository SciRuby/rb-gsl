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
      it "exposes cblas trsm, with B as a vector" do
        a = NMatrix.new(:dense, 3, [5,4,-1,0,10,-3,0,0,1], dtype)
        b = NVector.new(3, [0,11,3], dtype)
        alpha = 1
        NMatrix::BLAS::cblas_trsm(:row, :left, :upper, :no_transpose, :nonunit, b.shape[0], b.shape[1], alpha, a, 3, b, 1)
        b[0].should == -1
        b[1].should == 2
        b[2].should == 3
      end

      it "exposes cblas trsm, with B as a matrix"
    end
  end
end