require 'test_helper'

class OperTest < GSL::TestCase

  def test_multiplication_matrix
    matrix = GSL::Matrix.ones(1)

    mul_int   = 2 * matrix
    mul_float = 0.2 * matrix

    assert_equal 2,   mul_int[0, 0]
    assert_equal 0.2, mul_float[0, 0]
  end

  def test_multiplication_matrix_int
    matrix = GSL::Matrix::Int.ones(1)

    mul_int   = 2 * matrix
    mul_float = 0.2 * matrix

    assert_equal 2, mul_int[0, 0]
    assert_equal 0, mul_float[0, 0]
  end

  def test_multiplication_matrix_complex
    matrix = GSL::Matrix::Complex.eye(1)

    result = 0.2 * matrix

    assert_equal 0.2, result[0][0]
  end

  def test_multiplication_vector
    vector = GSL::Vector[1, 2]

    mul_int   = 2 * vector
    mul_float = 0.2 * vector

    assert_equal 2,   mul_int[0]
    assert_equal 0.2, mul_float[0]
  end

  def test_multiplication_vector_int
    vector = GSL::Vector::Int[1, 2]

    mul_int   = 2 * vector
    mul_float = 0.2 * vector

    assert_equal 2, mul_int[0]
    assert_equal 0, mul_float[0]
  end

  def test_multiplication_vector_complex
    re = GSL::Vector[1..4]
    im = GSL::Vector[5..8]

    vector = GSL::Vector::Complex[re, im]

    mul_int   = 2 * vector
    mul_float = 0.2 * vector

    assert_equal 10,  mul_int[0][1]
    assert_equal 1.0, mul_float[0][1]
  end

  def test_division_poly
    poly = GSL::Poly.alloc([2])

    a = GSL::Poly[1]; a[0] = 2

    result   = 2 / poly
    expected = GSL::Rational.new(a, poly)

    assert_equal expected.num, result.num
    assert_equal expected.den, result.den
  end

  def test_division_vector_col
    vector = GSL::Vector[1, 2].col

    result1 = 2 / vector
    result2 = 2 / result1

    assert_in_epsilon 0.4, result1[0]
    assert_equal result2, vector
  end

  def test_division_vector_int_col
    vector = GSL::Vector::Int[1, 2].col

    result1 = 2 / vector
    result2 = 2 / result1

    assert_in_epsilon 0.4, result1[0]
    assert_equal result2.to_a.map(&:to_i), vector.to_a
  end

end
