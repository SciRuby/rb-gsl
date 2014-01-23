require 'test_helper'

class LinalgTest < GSL::TestCase

  def _check(x, actual, eps, desc)
    if x == actual
      assert true
    elsif actual.zero?
      refute x.abs > eps, desc
    else
      refute((x - actual).abs / actual.abs > eps, desc)
    end
  end

  def _create_general_matrix(size1, size2)
    m = GSL::Matrix.alloc(size1, size2)

    size1.times { |i|
      size2.times { |j|
        m.set(i, j, 1.0 / (i + j + 1.0))
      }
    }

    m
  end

  def _create_hilbert_matrix(size)
    _create_general_matrix(size, size)
  end

  def _create_vandermonde_matrix(size)
    m = GSL::Matrix.alloc(size, size)

    size.times { |i|
      size.times { |j|
        m.set(i, j, GSL.pow(i + 1.0, size - j - 1.0))
      }
    }

    m
  end

  def setup
    @hilb2    = _create_hilbert_matrix(2)
    @hilb3    = _create_hilbert_matrix(3)
    @hilb4    = _create_hilbert_matrix(4)
    @hilb12   = _create_hilbert_matrix(12)
    @vander2  = _create_vandermonde_matrix(2)
    @vander3  = _create_vandermonde_matrix(3)
    @vander4  = _create_vandermonde_matrix(4)
    @vander12 = _create_vandermonde_matrix(12)

    @hilb2_solution  = GSL::Vector.alloc(-8.0, 18.0)
    @hilb3_solution  = GSL::Vector.alloc(27.0, -192.0, 210.0)
    @hilb4_solution  = GSL::Vector.alloc(-64.0, 900.0, -2520.0, 1820.0)
    @hilb12_solution = GSL::Vector.alloc(
      -1728.0, 245388.0, -8528520.0,
       127026900.0, -1009008000.0, 4768571808.0,
      -14202796608.0, 27336497760.0, -33921201600.0,
       26189163000.0, -11437874448.0, 2157916488.0
    )

    @vander2_solution  = GSL::Vector.alloc(1.0, 0.0)
    @vander3_solution  = GSL::Vector.alloc(0.0, 1.0, 0.0)
    @vander4_solution  = GSL::Vector.alloc(0.0, 0.0, 1.0, 0.0)
    @vander12_solution = GSL::Vector.alloc(0.0, 0.0, 0.0, 0.0,
                                           0.0, 0.0, 0.0, 0.0,
                                           0.0, 0.0, 1.0, 0.0)
  end

  def test_matmult
    a = GSL::Matrix.alloc([10.0, 5.0, 1.0, 20.0], 2, 2)
    b = GSL::Matrix.alloc([10.0, 5.0, 2.0, 1.0, 3.0, 2.0], 2, 3)
    c = a * b

    refute((c[0, 0] - 105.0).abs > GSL::DBL_EPSILON)
    refute((c[0, 1] -  65.0).abs > GSL::DBL_EPSILON)
    refute((c[0, 2] -  30.0).abs > GSL::DBL_EPSILON)
    refute((c[1, 0] -  30.0).abs > GSL::DBL_EPSILON)
    refute((c[1, 1] -  65.0).abs > GSL::DBL_EPSILON)
    refute((c[1, 2] -  42.0).abs > GSL::DBL_EPSILON)
  end

  def _test_bidiag_decomp_dim(m, eps, desc)
    eps *= 2 * GSL::DBL_EPSILON

    mm = m.size1
    nn = m.size2

    a = m.duplicate
    b = GSL::Matrix.calloc(nn, nn)

    u, v, d, sd = GSL::Linalg::Bidiag.unpack(*GSL::Linalg::Bidiag.decomp(a))

    b.set_diagonal(d)
    (nn - 1).times { |i| b[i, i + 1] = sd[i] }

    a = u * b * v.trans

    mm.times { |i|
      nn.times { |j|
        _check(aij = a[i, j], mij = m[i, j], eps,
          '%s: (%d,%d)[%d,%d]: %22.18g %22.18g' % [desc, mm, nn, i, j, aij, mij])
      }
    }
  end

  def test_bidiag_decomp
    m53 = _create_general_matrix(5, 3)
    m97 = _create_general_matrix(9, 7)

    _test_bidiag_decomp_dim(m53,       64.0, 'bidiag_decomp m(5,3)')
    _test_bidiag_decomp_dim(m97,       64.0, 'bidiag_decomp m(9,7)')
    _test_bidiag_decomp_dim(@hilb2,     8.0, 'bidiag_decomp hilbert(2)')
    _test_bidiag_decomp_dim(@hilb3,    64.0, 'bidiag_decomp hilbert(3)')
    _test_bidiag_decomp_dim(@hilb4,  1024.0, 'bidiag_decomp hilbert(4)')
    _test_bidiag_decomp_dim(@hilb12, 1024.0, 'bidiag_decomp hilbert(12)')
  end

  def test_cholesky
    m = GSL::Matrix.pascal(6)

    c_exp = GSL::Matrix[[1, 0, 0, 0, 0, 0],
                   [1, 1, 0, 0, 0, 0],
                   [1, 2, 1, 0, 0, 0],
                   [1, 3, 3, 1, 0, 0],
                   [1, 4, 6, 4, 1, 0],
                   [1, 5, 10, 10, 5, 1]]

    c = m.cholesky_decomp
    a = c.lower

    assert a == c_exp, "#{m.class}#cholesky_decomp"
    assert a * a.trans == m, "#{m.class}#cholesky_decomp"
  end

  def _test_HH_solve_dim(m, actual, eps, desc)
    eps *= GSL::DBL_EPSILON if eps > 1

    dim = m.size1

    x = GSL::Vector.indgen(dim) + 1
    GSL::Linalg::HH.svx(m.duplicate, x)

    dim.times { |i|
      _check(si = x[i], ai = actual[i], eps,
        '%s: %d[%d]: %22.18g %22.18g' % [desc, dim, i, si, ai])
    }
  end

  def test_HH_solve
    _test_HH_solve_dim(@hilb2,    @hilb2_solution,       8.0,  'HH_solve Hilbert(2)')
    _test_HH_solve_dim(@hilb3,    @hilb3_solution,     128.0,  'HH_solve Hilbert(3)')
    _test_HH_solve_dim(@hilb4,    @hilb4_solution,    2048.0,  'HH_solve Hilbert(4)')
    _test_HH_solve_dim(@hilb12,   @hilb12_solution,      0.5,  'HH_solve Hilbert(12)')
    _test_HH_solve_dim(@vander2,  @vander2_solution,     8.0,  'HH_solve Vander(2)')
    _test_HH_solve_dim(@vander3,  @vander3_solution,    64.0,  'HH_solve Vander(3)')
    _test_HH_solve_dim(@vander4,  @vander4_solution,  1024.0,  'HH_solve Vander(4)')
    _test_HH_solve_dim(@vander12, @vander12_solution,    0.05, 'HH_solve Vander(12)')
  end

  def test_LU
    m = GSL::Matrix.alloc([0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                          [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85])

    a = m.clone
    assert m == a, "#{a.class}#LU_decomp: matrix not destroyed"

    lu_exp = GSL::Matrix.alloc([0.51,              0.13,              0.19,              0.85],
                               [0.352941176470588, 0.554117647058823, 0.502941176470588, 0.66],
                               [0.803921568627451, 0.244515215852796, 0.71427813163482, -0.264713375796178],
                               [0.274509803921569, 0.476999292285916, 0.949126848480345, 0.363093705877982])

    x_exp = GSL::Vector[-4.05205022957397, -12.6056113959069, 1.66091162670884, 8.69376692879523]

    lu, perm, _sign = m.LU_decomp
    assert lu == lu_exp, "#{a.class}#LU_decomp"

    b = GSL::Vector[1, 2, 3, 4]
    x = GSL::Linalg::LU.solve(lu, perm, b)
    assert x == x_exp, "#{a.class}.LU_solve"

    x = lu.solve(perm, b)
    assert x == x_exp, "#{lu.class}#solve"

    perm, _sign = m.LU_decomp!
    assert m == lu_exp, "#{a.class}#LU_decomp!"

    m = a.clone

    x = GSL::Linalg::LU.solve(m, perm, b)
    assert x == x_exp, "#{a.class}.LU_solve"

    x = m.LU_solve(perm, b)
    assert x == x_exp, "#{a.class}#LU_solve"
    assert m == a, "#{a.class}#LU_solve: matrix not destroyed"

    h    = GSL::Matrix.hilbert(5)
    invh = GSL::Matrix.invhilbert(5)
    lu, perm, _sign = h.LU_decomp

    a = GSL::Linalg::LU.invert(lu, perm)
    assert a.equal?(invh, 1e-6), "#{h.class}#LU_invert, Hilbert matrix of order 5"

    a = lu.invert(perm)
    assert a.equal?(invh, 1e-6), "#{h.class}#LU_invert, Hilbert matrix of order 5"

    a = h.inv
    assert a.equal?(invh, 1e-6), "#{h.class}#LU_invert, Hilbert matrix of order 5"
  end

  def test_QR
    m = GSL::Matrix.alloc([0.18, 0.60, 0.57, 0.96], [0.41, 0.24, 0.99, 0.58],
                          [0.14, 0.30, 0.97, 0.66], [0.51, 0.13, 0.19, 0.85])

    a = m.clone
    assert m == a, "#{m.class}#QR_decomp: matrix not destroyed"

    x_exp = GSL::Vector[-4.05205022957397, -12.6056113959069, 1.66091162670884, 8.69376692879523]
    b = GSL::Vector[1, 2, 3, 4]

    qr, tau = m.QR_decomp

    x = m.QR_solve(b)
    assert x == x_exp, "#{m.class}#QR_solve(b)"

    x = GSL::Linalg::QR.solve(m, b)
    assert x == x_exp, 'GSL::Linalg::QR::solve(b)'

    tau = m.QR_decomp!
    assert m != a, "#{m.class}#QR_decomp: matrix destroyed"

    x = m.QR_solve(tau, b)
    assert x == x_exp, "#{m.class}#QR_solve(tau, b)"

    x = qr.solve(tau, b)
    assert x == x_exp, "#{qr.class}#solve(tau, b)"

    assert_raises(ArgumentError) { m.QR_solve(b) }
    assert_raises(ArgumentError) { m.solve(b) }

    x = m.solve(tau, b)
    assert x == x_exp, "#{m.class}#solve(tau, b)"

    m = a.clone
    bb = b.clone
    m.QR_svx(bb)
    assert bb == x_exp, "#{m.class}#QR_svx(b)"

    tau = GSL::Linalg::QR.decomp!(m)
    bb = b.clone
    m.QR_svx(tau, bb)
    assert bb == x_exp, "#{m.class}#QR_svx(tau, b)"
    assert_raises(ArgumentError) { m.QR_svx(bb) }

    m = a.clone
    qr, tau = m.QR_decomp
    assert m == a, "#{m.class}#QR_decomp: matrix not destroyed"

    x, r = m.QR_lssolve(b)
    assert x == x_exp, "#{m.class}#QR_lssolve(b)"

    r = m.QR_lssolve(b, x)
    assert x == x_exp, "#{qr.class}#QR_lssolve(b, x)"

    m.QR_lssolve(b, x, r)
    assert x == x_exp, "#{qr.class}#QR_lssolve(b, x, r)"

    x, r = qr.QR_lssolve(tau, b)
    assert x == x_exp, "#{qr.class}#QR_lssolve(tau, b)"

    r = qr.QR_lssolve(tau, b, x)
    assert x == x_exp, "#{qr.class}#QR_lssolve(tau, b, x)"

    qr.QR_lssolve(tau, b, x, r)
    assert x == x_exp, "#{qr.class}#QR_lssolve(tau, b, x, r)"
    assert_raises(ArgumentError) { qr.QR_lssolve(bb) }
  end

  def test_SV
    a = GSL::Matrix.alloc([1, 2, 3, 4], 2, 2)
    i = GSL::Matrix.identity(2)
    ainv = a.inv

    u, v, s = a.SV_decomp
    sm = s.to_m_diagonal
    sinv = s.map { |x| 1.0 / x }.to_m_diagonal

    assert u * sm * v.trans == a, "#{a.class}#SV_decomp"
    assert v * sinv * u.trans == ainv, "#{a.class}#SV_decomp"

    assert u.trans * u == i, "#{a.class}#SV_decomp"
    assert v.trans * v == i, "#{a.class}#SV_decomp"

    assert a * v == u * sm, "#{a.class}#SV_decomp"
    assert a.trans * u == v * sm, "#{a.class}#SV_decomp"
  end

  def _test_TDN_solve_dim(dim, d, a, b, actual, eps, desc)
    eps *= GSL::DBL_EPSILON

    abovediag = GSL::Vector.alloc(dim - 1)
    belowdiag = GSL::Vector.alloc(dim - 1)

    diag = GSL::Vector.alloc(dim)
    diag.set_all(d)

    rhs = GSL::Vector.indgen(dim) + 1

    abovediag.set_all(a)
    belowdiag.set_all(b)

    x = GSL::Linalg.solve_tridiag(diag, abovediag, belowdiag, rhs)

    dim.times { |i|
      _check(si = x[i], ai = actual[i], eps,
        '%s: %d[%d]: %22.18g %22.18g' % [desc, dim, i, si, ai])
    }
  end

  def test_TDN_solve
    actual = GSL::Vector.alloc(5)

    actual[0] = -7.0 / 3.0
    actual[1] =  5.0 / 3.0
    actual[2] =  4.0 / 3.0
    _test_TDN_solve_dim(3, 1.0, 2.0, 1.0, actual, 2.0, 'solve_TDN dim=2 A')

    actual[0] = 0.75
    actual[1] = 0.75
    actual[2] = 2.625
    _test_TDN_solve_dim(3, 1.0, 1.0 / 3.0, 1.0 / 2.0, actual, 2.0, 'solve_TDN dim=2 B')

    actual[0] =   99.0 / 140.0
    actual[1] =   41.0 /  35.0
    actual[2] =   19.0 /  10.0
    actual[3] =   72.0 /  35.0
    actual[4] =  139.0 /  35.0
    _test_TDN_solve_dim(5, 1.0, 1.0 / 4.0, 1.0 / 2.0, actual, 35.0 / 8.0, 'solve_TDN dim=5')
  end

  def _test_TDN_cyc_solve_dim(dim, d, a, b, actual, eps, desc)
    eps *= GSL::DBL_EPSILON

    abovediag = GSL::Vector.alloc(dim)
    belowdiag = GSL::Vector.alloc(dim)

    diag = GSL::Vector.alloc(dim)
    rhs = GSL::Vector.indgen(dim) + 1

    abovediag.set_all(a)
    belowdiag.set_all(b)

    diag.set_all(d)

    x = GSL::Linalg.solve_cyc_tridiag(diag, abovediag, belowdiag, rhs)
    dim.times { |i|
      _check(si = x[i], ai = actual[i], eps,
        '%s: %d[%d]: %22.18g %22.18g' % [desc, dim, i, si, ai])
    }
  end

  def test_TDN_cyc_solve
    actual = GSL::Vector.alloc(5)

    actual[0] =  3.0 / 2.0
    actual[1] = -1.0 / 2.0
    actual[2] =  1.0 / 2.0
    _test_TDN_cyc_solve_dim(3, 1.0, 2.0, 1.0, actual, 32.0, 'solve_TDN_cyc dim=2')

    actual[0] =  -5.0 / 22.0
    actual[1] =  -3.0 / 22.0
    actual[2] =  29.0 / 22.0
    actual[3] =  -9.0 / 22.0
    actual[4] =  43.0 / 22.0
    _test_TDN_cyc_solve_dim(5, 3.0, 2.0, 1.0, actual, 66.0, 'solve_TDN_cyc dim=5')
  end

  def _test_TDS_solve_dim(dim, d, od, actual, eps, desc)
    eps *= GSL::DBL_EPSILON

    diag = GSL::Vector.alloc(dim)
    diag.set_all(d)

    rhs = GSL::Vector.indgen(dim) + 1

    offdiag = GSL::Vector.alloc(dim - 1)
    offdiag.set_all(od)

    x = GSL::Linalg.solve_symm_tridiag(diag, offdiag, rhs)
    dim.times { |i|
      _check(si = x[i], ai = actual[i], eps,
        '%s: %d[%d]: %22.18g %22.18g' % [desc, dim, i, si, ai])
    }
  end

  def test_TDS_solve
    actual = GSL::Vector[0.0, 2.0]
    _test_TDS_solve_dim(2, 1.0, 0.5, actual, 8.0, 'solve_TDS dim=2 A')

    actual = GSL::Vector[3.0 / 8.0, 15.0 / 8.0]
    _test_TDS_solve_dim(2, 1.0, 1.0 / 3.0, actual, 8.0, 'solve_TDS dim=2 B')

    actual = GSL::Vector[5.0 / 8.0, 9.0 / 8.0, 2.0, 15.0 / 8.0, 35.0 / 8.0]
    _test_TDS_solve_dim(5, 1.0, 1.0 / 3.0, actual, 8.0, 'solve_TDS dim=5')
  end

  def _test_TDS_cyc_solve_one(dim, d, od, r, actual, eps, desc)
    eps *= GSL::DBL_EPSILON

    diag = d.duplicate
    offdiag = od.duplicate
    rhs = r.duplicate

    x = GSL::Linalg.solve_symm_cyc_tridiag(diag, offdiag, rhs)
    dim.times { |i|
      _check(si = x[i], ai = actual[i], eps,
        '%s: %d[%d]: %22.18g %22.18g' % [desc, dim, i, si, ai])
    }
  end

  def test_TDS_cyc_solve
    diag = GSL::Vector.alloc(1)
    diag[0] = 2

    offdiag = GSL::Vector.alloc(1)
    offdiag[0] = 3

    rhs = GSL::Vector.alloc(1)
    rhs[0] = 7

    actual = GSL::Vector.alloc(1)
    actual[0] = 3.5

    # XXX GSL::ERROR::EBADLEN: Ruby/GSL error code 19, size of cyclic system must be
    # 3 or more (file tridiag.c, line 531), matrix/vector sizes are not conformant
    #_test_TDS_cyc_solve_one(1, diag, offdiag, rhs, actual, 28.0, 'solve_TDS_cyc dim=1')

    diag = GSL::Vector[1, 2]
    offdiag = GSL::Vector[3, 4]
    rhs = GSL::Vector[7, -7]
    actual = GSL::Vector[-5, 4]

    # XXX GSL::ERROR::EBADLEN: Ruby/GSL error code 19, size of cyclic system must be
    # 3 or more (file tridiag.c, line 531), matrix/vector sizes are not conformant
    #_test_TDS_cyc_solve_one(2, diag, offdiag, rhs, actual, 28.0, 'solve_TDS_cyc dim=2')

    diag = GSL::Vector[1, 1, 1]
    offdiag = GSL::Vector[3, 3, 3]
    rhs = GSL::Vector[7, -7, 7]
    actual = GSL::Vector[-2, 5, -2]

    _test_TDS_cyc_solve_one(3, diag, offdiag, rhs, actual, 28.0, 'solve_TDS_cyc dim=3')

    diag = GSL::Vector[4, 2, 1, 2, 4]
    offdiag = GSL::Vector[1, 1, 1, 1, 1]
    rhs = GSL::Vector[30, -24, 3, 21, -30]
    actual = GSL::Vector[12, 3, -42, 42, -21]

    _test_TDS_cyc_solve_one(5, diag, offdiag, rhs, actual, 35.0, 'solve_TDS_cyc dim=5')
  end

end
