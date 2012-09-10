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
# == math_spec.rb
#
# Tests for non-BLAS and non-LAPACK math functions, or for simplified
# versions of unfriendly BLAS and LAPACK functions.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")

describe "math" do

  [:float32, :float64, :complex64, :complex128, :rational32, :rational64, :rational128].each do |dtype|
    context dtype do
      it "should correctly factorize a matrix" do
        m = NMatrix.new(:dense, 3, [4,9,2,3,5,7,8,1,6], dtype)
        a = m.factorize_lu
        a[0,0].should == 8
        a[0,1].should == 1
        a[0,2].should == 6
        a[1,0].should == 0.5
        a[1,1].should == 8.5
        a[1,2].should == -1
        a[2,0].should == 0.375
      end
    end
  end

  [:float32, :float64, :complex64, :complex128].each do |dtype|
    context dtype do

      # Note: this exposes gemm, not cblas_gemm (which is the unfriendly CBLAS no-error-checking version)
      it "exposes gemm" do
        #STDERR.puts "dtype=#{dtype.to_s}"
        #STDERR.puts "1"
        n = NMatrix.new([4,3], dtype)
        n[0,0] = 14.0
        n[0,1] = 9.0
        n[0,2] = 3.0
        n[1,0] = 2.0
        n[1,1] = 11.0
        n[1,2] = 15.0
        n[2,0] = 0.0
        n[2,1] = 12.0
        n[2,2] = 17.0
        n[3,0] = 5.0
        n[3,1] = 2.0
        n[3,2] = 3.0

        m = NMatrix.new([3,2], dtype)

        m[0,0] = 12.0
        m[0,1] = 25.0
        m[1,0] = 9.0
        m[1,1] = 10.0
        m[2,0] = 8.0
        m[2,1] = 5.0

        #c = NMatrix.new([4,2], dtype)
        r = NMatrix::BLAS.gemm(n, m) #, c)
        #c.should equal(r) # check that both are same memory address

        r[0,0].should == 273.0
        r[0,1].should == 455.0
        r[1,0].should == 243.0
        r[1,1].should == 235.0
        r[2,0].should == 244.0
        r[2,1].should == 205.0
        r[3,0].should == 102.0
        r[3,1].should == 160.0
      end

    end
  end
end