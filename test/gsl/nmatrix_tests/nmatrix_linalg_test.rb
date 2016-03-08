require 'test_helper'

class NMatrixLinalgTest < GSL::TestCase
  def setup
    @nm = NMatrix.new([4,4], 
      [
        0.18, 0.60, 0.57, 0.96, 
        0.41, 0.24, 0.99, 0.58,
        0.14, 0.30, 0.97, 0.66, 
        0.51, 0.13, 0.19, 0.85
      ], dtype: :float64)
    @b = NMatrix.new([4], [1,2,3,4], dtype: :float64)
    @x_exp = NMatrix.new([4], [-4.05205022957397, -12.6056113959069, 1.66091162670884, 8.69376692879523], dtype: :float64)   
  end

  def test_lu
    lu, perm, signum = GSL::Linalg::LU.decomp(@nm)

    lu_exp = NMatrix.new([4,4],
     [0.51,              0.13,              0.19,              0.85,
      0.352941176470588, 0.554117647058823, 0.502941176470588, 0.66,
      0.803921568627451, 0.244515215852796, 0.71427813163482, -0.264713375796178,
      0.274509803921569, 0.476999292285916, 0.949126848480345, 0.363093705877982],
      dtype: :float64)

    assert_enum_abs lu, lu_exp, 0.0001, "GSL::Linalg::LU.decomp(A) with NMatrix"
    assert_enum_abs GSL::Linalg::LU.solve(lu, perm, @b), @x_exp, 0.0001, "GSL::Linalg::LU.solve(lu, perm, b) with NMatrix"

    ##########################################################################

    nm = NMatrix.new([2,2], [1,1,14,1], dtype: :float64)
    lu, perm, sign = GSL::Linalg::LU.decomp(nm)
    inverted  = NMatrix.new([2,2], 
      [
        -0.076923,  0.076923,
         1.076923, -0.076923
      ], dtype: :float64)

    assert_enum_abs GSL::Linalg::LU.invert(lu, perm), inverted, 0.001, "GSL::Linalg::LU.invert(lu, perm) with NMatrix"
    assert GSL::Linalg::LU.det(lu, sign) == -13, "GSL::Linalg::LU.det(lu, sign) with NMatrix"
  end

  def test_qr
    qr_answer = NMatrix.new([4,4], 
      [-0.692965, -0.454136, -1.06961 , -1.35144, 
       0.469664 ,  0.564146,  0.72597 , 0.726696, 
       0.160373 , -0.159838, -0.781606, 0.063932, 
       0.584216 ,  0.593044, -0.332286, 0.239865], dtype: :float64)
    tau_answer = NMatrix.new([4], [1.25975, 1.45217, 1.80113, 0.0], dtype: :float64)
    qr, tau = GSL::Linalg::QR.decomp(@nm)
    
    assert_enum_abs qr_answer , qr , 0.001, "GSL::Linalg::QR.decomp(nmatrix)"
    assert_enum_abs tau_answer, tau, 0.001, "GSL::Linalg::QR.decomp(nmatrix)"

    assert_enum_abs GSL::Linalg::QR.solve(qr, tau, @b), @x_exp, 0.001, 
      "GSL::Linalg::QR.solve(qr, tau, b)"
  end

  def test_sv
    u_answer = NMatrix.new([4,4],
      [-0.545591,-0.312561, 0.709796,-0.317529, 
       -0.5298  , 0.418583,-0.475268,-0.564111, 
       -0.524621, 0.436573, 0.112083, 0.722229, 
       -0.382642,-0.73246 ,-0.507688, 0.243598 ], dtype: :float64)
    v_answer = NMatrix.new([4,4],
      [-0.260024,-0.288729 ,-0.732277,-0.559279 , 
       -0.294582,-0.0751933, 0.659393,-0.687581 , 
       -0.630928, 0.762631 ,-0.12665 , 0.0654517, 
       -0.668983,-0.573912 , 0.113712, 0.458427 ], dtype: :float64)
    s_answer = NMatrix.new([4], [2.24602,0.682566,0.423782,0.112813], dtype: :float64)

    u, v, s = GSL::Linalg::SV.decomp(@nm)

    assert_enum_abs u, u_answer, 0.001, "GSL::Linalg::SV.decomp(nmatrix) -> u"
    assert_enum_abs v, v_answer, 0.001, "GSL::Linalg::SV.decomp(nmatrix) -> v"
    assert_enum_abs s, s_answer, 0.001, "GSL::Linalg::SV.decomp(nmatrix) -> s"

    assert_enum_abs GSL::Linalg::SV.solve(u, v, s, @b), @x_exp, 0.001, 
      "GSL::Linalg::SV.solve(u,v,s,b)"
  end

  def test_cholesky
    m = NMatrix.new([2,2], [4.0, 2, 2, 3], dtype: :float64)
    b = NMatrix.new([2], [1,2], dtype: :float64)

    cholesky = NMatrix.new([2,2], [2.0, 1.0, 1.0, 1.41421], dtype: :float64)
    x_exp = NMatrix.new([2], [-0.125, 0.75], dtype: :float64)

    c = GSL::Linalg::Cholesky
    assert_enum_abs c.decomp(m)         , cholesky, 0.001, "GSL::Linalg::Cholesky.decomp"
    assert_enum_abs c.solve(cholesky, b), x_exp   , 0.001, "GSL::Linalg::Cholesky.solve"
    assert_enum_abs c.svx(cholesky, b)  , x_exp   , 0.001, "GSL::Linalg::Cholesky.svx"
  end

  def test_hh
    hh = GSL::Linalg::HH
    assert_enum_abs hh.solve(@nm, @b), @x_exp, 0.001, "GSL::Linalg::HH.solve(m, b)"
    assert_enum_abs hh.svx(@nm, @b)  , @x_exp, 0.001, "GSL::Linalg::HH.svx(m, b)"
  end
end