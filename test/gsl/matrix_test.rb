require 'test_helper'

class MatrixTest < GSL::TestCase

  def test_ispos_neg
    m = GSL::Matrix::Int.alloc([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)
    assert_equal 0,     m.ispos
    assert_equal false, m.ispos?
    assert_equal 0,     m.isneg
    assert_equal false, m.isneg?

    m += 1
    assert_equal 1,     m.ispos
    assert_equal true,  m.ispos?
    assert_equal 0,     m.isneg
    assert_equal false, m.isneg?

    m -= 100
    assert_equal 0,     m.ispos
    assert_equal false, m.ispos?
    assert_equal 1,     m.isneg
    assert_equal true,  m.isneg?
  end

  def test_isnonneg
    m = GSL::Matrix::Int.alloc([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 3, 3)
    assert_equal 1,     m.isnonneg
    assert_equal true,  m.isnonneg?
    assert_equal 0,     m.isneg
    assert_equal false, m.isneg?

    m -= 100
    assert_equal 0,     m.isnonneg
    assert_equal false, m.isnonneg?
    assert_equal 1,     m.isneg
    assert_equal true,  m.isneg?

    m += 200
    assert_equal 1,     m.isnonneg
    assert_equal true,  m.isnonneg?
    assert_equal 1,     m.ispos
    assert_equal true,  m.ispos?
  end

  def test_eye
    z = GSL::Complex[1, 0]
    m = GSL::Matrix::Complex.eye(2, z)

    assert_equal z, m[0, 0]
    assert_equal GSL::Complex[0, 0], m[0, 1]
    assert_equal GSL::Complex[0, 0], m[1, 0]
    assert_equal z, m[1, 1]
  end

  def test_set_row
    z0 = GSL::Complex[1, 0]
    z1 = GSL::Complex[2, 0]

    m = GSL::Matrix::Complex[2, 2]
    m.set_row(0, z0, z1)

    assert_equal z0, m[0, 0]
    assert_equal z1, m[0, 1]
  end

  def test_set_col
    z0 = GSL::Complex[1, 0]
    z1 = GSL::Complex[2, 0]

    m = GSL::Matrix::Complex[2, 2]
    m.set_col(0, z0, z1)

    assert_equal z0, m[0, 0]
    assert_equal z1, m[1, 0]
  end

end
