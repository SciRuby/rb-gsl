require 'minitest/autorun'
require 'minitest/unit'

require 'gsl'

class Sf < MiniTest::Unit::TestCase

  def test_mode
    # TODO!!!
    z = GSL::Complex[1,0]
    m = GSL::Matrix::Complex.eye(2, z)
    assert_equal(z, m[0,0])
    assert_equal(GSL::Complex[0,0], m[0,1])
    assert_equal(GSL::Complex[0,0], m[1,0])
    assert_equal(z, m[1,1])
  end

end

