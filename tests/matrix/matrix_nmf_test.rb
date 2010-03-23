#!/usr/bin/env ruby
# Written by Roman Shterenzon
#
#$:.unshift(*['lib','ext'].collect{|d| File.join(File.dirname(__FILE__),'..','..',d)})
require 'rubygems'
require 'test/unit'
require 'gsl'

class MatrixNmfTest < Test::Unit::TestCase

  def setup
    @m1 = GSL::Matrix.alloc([6, 7, 8], [2, 3, 4], [3, 4, 5])
    @m2 = GSL::Matrix.alloc([6, 7, 8], [2, 3, 4], [3, 4, 7])
  end

  def test_difcost
    assert_equal(0, GSL::Matrix::NMF.difcost(@m1, @m1))
    assert_equal(4, GSL::Matrix::NMF.difcost(@m1, @m2))
  end

  def test_nmf
    [2, 3, 4, 5].each do |cols|
      res = GSL::Matrix::NMF.nmf(@m1, cols)
      assert_equal([3,cols], res[0].size)
      assert_equal([cols,3], res[1].size)
      cost = GSL::Matrix::NMF.difcost(@m1, res[0]*res[1])
      assert(cost <= 0.000001, "Cols: #{cols}, Delta: #{cost}")
    end
  end
  def test_matrix_nmf
    [2, 3, 4, 5].each do |cols|
      res = @m1.nmf(cols)
      assert_equal([3,cols], res[0].size)
      assert_equal([cols,3], res[1].size)
      cost = GSL::Matrix::NMF.difcost(@m1, res[0]*res[1])
      assert(cost <= 0.000001, "Cols: #{cols}, Delta: #{cost}")
    end
  end
end
