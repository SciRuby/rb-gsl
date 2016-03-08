require 'test_helper'

class NMatrixGslTest < GSL::TestCase
  def setup
    @gsl_vector = GSL::Vector.alloc(2.354, 4.443, 6.776)
    @nm_vector  = NMatrix.new([3], [2.354, 4.443, 6.776], dtype: :float64)

    @gsl_int_vector = GSL::Vector::Int[1,2,3,4,5]
    @nm_int_vector  = NMatrix.new([5], [1,2,3,4,5], dtype: :int32)

    @gsl_complex_vector = GSL::Vector::Complex.alloc([[1,0], [2,0], [3,0]])
    @nm_complex_vector  = NMatrix.new([3], [1,2,3], dtype: :complex128)

    @gsl_matrix = GSL::Matrix.alloc(
      [1,2,3],
      [4,5,6],
      [7,8,9]
    )
    @nm_matrix  = NMatrix.new([3,3], [1,2,3,4,5,6,7,8,9], dtype: :float64)

    @gsl_int_matrix = GSL::Matrix::Int.alloc(
      [1,2,3],
      [4,5,6],
      [7,8,9] 
    )
    @nm_int_matrix = NMatrix.new([3,3], [1,2,3,4,5,6,7,8,9], dtype: :int32)

    @gsl_complex_matrix = GSL::Matrix::Complex.alloc(2,2)
    @gsl_complex_matrix.set(0,0, [1.1,1.1])
    @gsl_complex_matrix.set(0,1, [2.2,2.2])
    @gsl_complex_matrix.set(1,0, [3.3,3.3])
    @gsl_complex_matrix.set(1,1, [4.4,4.4])

    @nm_complex_matrix = NMatrix.new([2,2], 
      [Complex(1.1,1.1), Complex(2.2,2.2), Complex(3.3,3.3), Complex(4.4,4.4)], dtype: :complex128)
  end

  # GSL::Vector to 1D NMatrix
  def test_gsl_vector_to_nmatrix
    assert_equal @nm_vector        , @gsl_vector.to_nm        , 'floating point GSL::Vector to NMatrix'
    assert_equal @nm_int_vector    , @gsl_int_vector.to_nm    , 'integer GSL::Vector to NMatrix'
    assert_equal @nm_complex_vector, @gsl_complex_vector.to_nm, 'complex GSL::Vector to NMatrix'
  end

  # GSL::Matrix to NMatrix
  def test_gsl_matrix_to_nmatrix
    assert_equal @nm_matrix, @gsl_matrix.to_nm, 'floating point GSL::Matrix to 2D NMatrix'
    assert_equal @nm_int_matrix, @gsl_int_matrix.to_nm, 'GSL::Matrix::Int to 2D NMatrix'
    assert_equal @nm_complex_matrix, @gsl_complex_matrix.to_nm, 'GSL::Matrix::Complex to 2D NMatrix'
  end

  # NMatrix to GSL::Vector
  def test_nmatrix_to_gsl_vector
    assert_equal @gsl_vector        , @nm_vector.to_gslv        , 'floating point NMatrix to GSL::Vector'
    assert_equal @gsl_int_vector    , @nm_int_vector.to_gslv    , 'int NMatrix to GSL::Vector::Int'
    assert_equal @gsl_complex_vector, @nm_complex_vector.to_gslv, 'complex NMatrix to GSL::Vector::Complex'
  end

  # NMatrix to GSL::Matrix
  def test_nmatrix_to_gsl_matrix
    assert_equal @gsl_matrix        , @nm_matrix.to_gslm        , 'floating NMatrix to GSL::Matrix'
    assert_equal @gsl_int_matrix    , @nm_int_matrix.to_gslm    , 'int NMatrix to GSL::Matrix::Int'
    assert_equal @gsl_complex_matrix, @nm_complex_matrix.to_gslm, 'complex NMatrix to GSL::Matrix::Complex'
  end
end