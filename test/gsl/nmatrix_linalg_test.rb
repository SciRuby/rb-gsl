require 'test_helper'

class NMatrixGslTest < GSL::TestCase
  def test_lu
    nm = NMatrix.new([4,4], 
      [
        [0.18, 0.60, 0.57, 0.96], 
        [0.41, 0.24, 0.99, 0.58],
        [0.14, 0.30, 0.97, 0.66], 
        [0.51, 0.13, 0.19, 0.85]
      ], dtype: :float64)
    b = NMatrix.new([4], [1,2,3,4], dtype: :float64)
    lu, perm, signum = GSL::Linalg::LU.decomp(nm)

    lu_exp = NMatrix.new([4,4],
     [0.51,              0.13,              0.19,              0.85,
      0.352941176470588, 0.554117647058823, 0.502941176470588, 0.66,
      0.803921568627451, 0.244515215852796, 0.71427813163482, -0.264713375796178,
      0.274509803921569, 0.476999292285916, 0.949126848480345, 0.363093705877982],
      dtype: :float64)
    x_exp = NMatrix.new([4], [-4.05205022957397, -12.6056113959069, 1.66091162670884, 8.69376692879523], dtype: :flaot64)

    assert lu == lu_exp, "GSL::Linalg::LU.decomp(A) with NMatrix"
    assert GSL::Linalg::LU.solve(lu, b) == x_exp, "GSL::Linalg::LU.solve(lu, b) with NMatrix"

    ##########################################################################

    nmatrix = NMatrix.new([2,2], [1,1,14,1], dtype: :float64)
    lu, perm, sign = GSL::Linalg::LU.decomp(nmatrix)
    inverted  = NMatrix.new([2,2], 
      [-0.076923, 0.076923, 1.076923, -0.076923], dtype: :float64)

    assert GSL::Linalg::LU.invert(lu, perm) == inverted, "GSL::Linalg::LU.invert(lu, perm) with NMatrix"
    assert GSL::Linalg::LU.det(lu, sign)    == -13, "GSL::Linalg::LU.det(lu, sign) with NMatrix"

    # TODO: Write tests for GSL::Linalg::LU.svx(lu, bx) and lndet(lu)
  end

  def test_qr
    
  end
end