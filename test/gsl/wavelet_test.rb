require 'test_helper'

class WaveletTest < GSL::TestCase

  MEMBERS = [309, 307, 305, 303, 301, 208, 206, 204, 202, 105, 103]

  def _urand
    x = 1
    x = (1103515245 * x + 12345) & 0x7fffffff
    x / 2147483648.0
  end

  def _test_1d(n, stride, type, member)
    nn = n * stride
    data = GSL::Vector.alloc(nn)
    nn.times { |i| data[i] = 12345.0 + i }

    v1 = data.view_with_stride(0, stride, n)
    n.times { |i| v1[i] = _urand }

    v2 = GSL::Vector.alloc(n)
    GSL::Vector.memcpy(v2, v1)

    vdelta = GSL::Vector.alloc(n)

    work = GSL::Wavelet::Workspace.alloc(n)
    w = GSL::Wavelet.alloc(type, member)
    w.transform_forward(v2, work)
    w.transform_inverse(v2, work)

    n.times { |i| vdelta[i] = (v1[i] - v2[i]).abs }

    i = vdelta.max_index
    x1, x2 = v1[i], v2[i]

    refute((x2 - x1).abs > n * 1e-15,
      "#{w.name}(#{member}), n = #{n}, stride = #{stride}, maxerr = #{(x2 - x1).abs}")

    assert((0...nn).all? { |j| j % stride == 0 || data[j] == 12345.0 + j },
      "#{w.name}(#{member}) other data untouched, n = #{n}, stride = #{stride}") if stride > 1
  end

  def _test_2d(n, tda, t, member, type)
    nn = n * tda
    data = GSL::Vector.alloc(nn)
    nn.times { |i| data[i] = 12345.0 + i }

    m1 = data.matrix_view_with_tda(n, n, tda)
    n.times { |i| n.times { |j| m1.set(i, j, _urand) } }

    m2 = GSL::Matrix.alloc(n, n)
    GSL::Matrix.memcpy(m2, m1)

    mdelta = GSL::Matrix.alloc(n, n)

    work = GSL::Wavelet::Workspace.alloc(n)
    w = GSL::Wavelet.alloc(t, member)

    typename = case type
      when 1
        GSL::Wavelet2d.transform_matrix_forward(w, m2, work)
        GSL::Wavelet2d.transform_matrix_inverse(w, m2, work)
        'standard'
      when 2
        GSL::Wavelet2d.nstransform_matrix_forward(w, m2, work)
        GSL::Wavelet2d.nstransform_matrix_inverse(w, m2, work)
        'nonstd'
    end

    n.times { |i| n.times { |j| mdelta.set(i, j, (m1[i, j] - m2[i, j]).abs) } }

    i, j = mdelta.max_index
    x1, x2 = m1[i, j], m1[i, j]

    refute((x2 - x1).abs > n * 1e-15,
      "#{w.name}(#{member})-2d #{typename}, n = #{n}, tda = #{tda}, maxerr = #{(x2 - x1).abs}")

    assert((0...n).to_a.product((n...tda).to_a).all? { |k, l| data[k * tda + l] == 12345.0 + k * tda + l },
      "#{w.name}(#{member})-2d #{typename} other data untouched, n = #{n}, tda = #{tda}") if tda > n
  end

  def _each_pow(n = 14)
    n.times { |i| yield 2 ** i }
  end

  def _each_n(n = 9, m = 4)
    n.times { |i| yield m + 2 * i }
  end

  def test_1d_bspline
    _each_pow { |n| MEMBERS.each { |m| _test_1d(n, 1, 'bspline', m) } }
  end

  def test_1d_bspline_centered
    _each_pow { |n| MEMBERS.each { |m| _test_1d(n, 1, 'bspline_centered', m) } }
  end

  def test_1d_daubechies
    _each_pow { |n| _each_n { |i| _test_1d(n, 1, 'daubechies', i) } }
  end

  def test_1d_daubechies_centered
    _each_pow { |n| _each_n { |i|_test_1d(n, 1, 'daubechies_centered', i) } }
  end

  def test_1d_haar
    _each_pow { |n| _test_1d(n, 1, 'haar', 2) }
  end

  def test_1d_haar_centered
    _each_pow { |n| _test_1d(n, 1, 'haar_centered', 2) }
  end

  def test_2d_bspline_standard
    _each_pow(6) { |n| MEMBERS.each { |m| _test_2d(n, n, 'bspline', m, 1) } }
  end

  def test_2d_bspline_centered_standard
    _each_pow(6) { |n| MEMBERS.each { |m| _test_2d(n, n, 'bspline_centered', m, 1) } }
  end

  def test_2d_bspline_nonstd
    _each_pow(6) { |n| MEMBERS.each { |m| _test_2d(n, n, 'bspline', m, 2) } }
  end

  def test_2d_bspline_centered_nonstd
    _each_pow(6) { |n| MEMBERS.each { |m| _test_2d(n, n, 'bspline_centered', m, 2) } }
  end

end
