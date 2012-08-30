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
# == nmatrix_spec.rb
#
# Element-wise operation tests.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")

describe NMatrix do

  context "list" do
    before :each do
      @n = NMatrix.new(:list, 2, 0, :int64)
      @m = NMatrix.new(:list, 2, 0, :int64)
      @n[0,0] = 52
      @m[0,0] = -48
      @n[1,1] = 40
    end

    it "adds" do
      r = NMatrix.new(:dense, 2, [4, 0, 0, 40], :int64).cast(:list, :int64)
      (@n+@m).should == r
    end

    it "subtracts" do
      r = NMatrix.new(:dense, 2, [100, 0, 0, 40], :int64).cast(:list, :int64)
      (@n-@m).should == r
    end

    it "multiplies" do
      r = NMatrix.new(:dense, 2, [-2496, 0, 0, 0], :int64).cast(:list, :int64)
      (@n*@m).should == r
    end

    it "divides" do
      m = @m.clone
      m[1,1] = 2
      r = NMatrix.new(:dense, 2, [-2, 0, 0, 1], :int64).cast(:list, :int64)
      (@n/@m).should == r
    end

    it "modulo"

    it "equals" do
      (@n =~ @m).should == NMatrix.new(:list, 2, 0, :byte)
    end

    it "is not equal to" do
      (@n !~ @m).should == NMatrix.new(:list, 2, 1, :byte)
    end

    it "is less than" do
      (@n < @m).should == NMatrix.new(:list, 2, 0, :byte)
    end

    it "is greater than" do
      (@n > @m).should == NMatrix.new(:dense, 2, [1, 0, 0, 1], :byte).cast(:list, :byte)
    end

    it "is greater than or equal to" do
      (@n >= @m).should == NMatrix.new(:dense, 2, [1, 1, 1, 1], :byte).cast(:list, :byte)
    end

    it "is less than or equal to" do
      (@n <= @m).should == NMatrix.new(:dense, 2, [0, 1, 1, 0], :byte).cast(:list, :byte)
    end
  end

  context "dense" do
    context "elementwise arithmetic" do
      before :each do
        @n = NMatrix.new(:dense, 2, [1,2,3,4], :int64)
        @m = NMatrix.new(:dense, 2, [-4,-1,0,66], :int64)
      end

      it "adds" do
        r = @n+@m
        r.should == NMatrix.new(:dense, [2,2], [-3, 1, 3, 70], :int64)
      end

      it "subtracts" do
        r = @n-@m
        r.should == NMatrix.new(:dense, [2,2], [5, 3, 3, -62], :int64)
      end

      it "multiplies" do
        r = @n*@m
        r.should == NMatrix.new(:dense, [2,2], [-4, -2, 0, 264], :int64)
      end

      it "divides in the Ruby way" do
        m = @m.clone
        m[1,0] = 3
        r = @n/m
        r.should == NMatrix.new(:dense, [2,2], [-1, -2, 1, 0], :int64)
      end

      it "modulo" do
        r = @n % @m
        r.should == NMatrix.new(:dense, [2,2], [0, 1, 2, 3], :int64)
      end
    end

    context "elementwise comparisons" do
      before :each do
        @n = NMatrix.new(:dense, 2, [1,2,3,4], :int64)
        @m = NMatrix.new(:dense, 2, [-4,-1,3,2], :int64)
      end

      it "equals" do
        r = @n =~ @m
        r.should == NMatrix.new(:dense, [2,2], [0, 0, 1, 0], :byte)
      end

      it "is not equal" do
        r = @n !~ @m
        r.should == NMatrix.new(:dense, [2,2], [1, 1, 0, 1], :byte)
      end

      it "is less than" do
        r = @n < @m
        r.should == NMatrix.new(:dense, [2,2], 0, :byte)
      end

      it "is greater than" do
        r = @n > @m
        r.should == NMatrix.new(:dense, [2,2], [1, 1, 0, 1], :byte)
      end

      it "is less than or equal to" do
        r = @n <= @m
        r.should == NMatrix.new(:dense, [2,2], [0, 0, 1, 0], :byte)
      end

      it "is greater than or equal to" do
        n = NMatrix.new(:dense, [2,2], [ 1,  2, 2, 4], :int64)
        r = n >= @m
        r.should == NMatrix.new(:dense, [2,2], [1, 1, 0, 1], :byte)
      end
    end
  end
end