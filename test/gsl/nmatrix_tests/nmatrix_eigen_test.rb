require 'test_helper'

class NMatrixEigenTest < GSL::TestCase
  def setup
    @nmatrix = NMatrix.new([4,4],
      [1.0, 1/2.0, 1/3.0, 1/4.0, 
       1/2.0, 1/3.0, 1/4.0, 1/5.0,
       1/3.0, 1/4.0, 1/5.0, 1/6.0,
       1/4.0, 1/5.0, 1/6.0, 1/7.0], dtype: :float64)
  end

  def test_symm_symmv
    eigen_values = NMatrix.new([4], 
      [1.50021, 0.169141, 0.00673827, 9.67023e-05], dtype: :float64)
    eigen_vectors = NMatrix.new([4,4],
      [0.792608, 0.582076,-0.179186,-0.0291933, 
       0.451923,-0.370502, 0.741918, 0.328712 , 
       0.322416,-0.509579,-0.100228,-0.791411 , 
       0.252161,-0.514048,-0.638283, 0.514553], dtype: :float64)
    
    assert_enum_abs GSL::Eigen.symm(@nmatrix), eigen_values, 0.001, "GSL::Eigen.symm(nmatrix)"

    # val, vec = GSL::Eigen.symmv(@nmatrix)

    # assert_enum_abs val, eigen_values , 0.001, "GSL::Eigen.symmv(nmatrix)"
    # assert_enum_abs vec, eigen_vectors, 0.001, "GSL::Eigen.symmv(nmatrix)"
  end
end
