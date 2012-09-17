require 'minitest/autorun'
require 'minitest/unit'

require 'gsl'

class MatrixComplex < MiniTest::Unit::TestCase

  def test_eye
    z = GSL::Complex[1,0]
    m = GSL::Matrix::Complex.eye(2, z)
    assert_equal(z, m[0,0])
    assert_equal(GSL::Complex[0,0], m[0,1])
    assert_equal(GSL::Complex[0,0], m[1,0])
    assert_equal(z, m[1,1])
  end

  def test_set_row
    z0 = GSL::Complex[1,0]
    z1 = GSL::Complex[2,0]
    m = GSL::Matrix::Complex[2,2]
    m.set_row(0,z0,z1)
    assert_equal(z0, m[0,0])
    assert_equal(z1, m[0,1])
  end

  def test_set_col
    z0 = GSL::Complex[1,0]
    z1 = GSL::Complex[2,0]
    m = GSL::Matrix::Complex[2,2]
    m.set_col(0,z0,z1)
    assert_equal(z0, m[0,0])
    assert_equal(z1, m[1,0])
  end

end

