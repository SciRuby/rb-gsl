#!/usr/bin/env ruby
require("gsl")
require("./gsl_test2.rb")
include GSL::Test
include Math

N_BS = 11

def urand()
  x = 1
  x = (1103515245 * x + 12345) & 0x7fffffff
  return x / 2147483648.0
end

def test_1d(n, stride, type, member)
  nn = n*stride
  data = GSL::Vector.alloc(nn)
  for i in 0...nn
    data[i] = 12345.0 + i
  end
  v1 = data.view_with_stride(0, stride, n)
  for i in 0...n
    v1[i] = urand()
  end

  v2 = GSL::Vector.alloc(n)
  GSL::Vector.memcpy(v2, v1)

  vdelta = GSL::Vector.alloc(n)

  work = GSL::Wavelet::Workspace.alloc(n)
  w = GSL::Wavelet.alloc(type, member)
  w.transform_forward(v2, work)
  w.transform_inverse(v2, work)
  for i in 0...n
    x1 = v1[i]
    x2 = v2[i]
    vdelta[i] = (x1 - x2).abs
  end
  i = vdelta.max_index
  x1 = v1[i]
  x2 = v2[i]
  GSL::Test::test((x2 - x1).abs > n*1e-15 ? 1 : 0, "#{w.name}(#{member}), n = #{n}, stride = #{stride}, maxerr = #{(x2 - x1).abs}")

  if stride > 1
    status = 0
    for i in 0...nn
      next if i%stride == 0
      status |= (data[i] != 12345.0 + i) ? 1 : 0
    end
    GSL::Test::test(status, "#{w.name}(#{member}) other data untouched, n = #{n}, stride = #{stride}")
  end
end

def test_2d(n, tda, t, member, type)
  nn = n*tda
  data = GSL::Vector.alloc(nn)
  typename = (type == 1) ? "standard" : "nonstd"
  for i in 0...nn
    data[i] = 12345.0 + i
  end
  m1 = data.matrix_view_with_tda(n, n, tda)
  for i in 0...n
    for j in 0...n
      m1.set(i, j, urand())
    end
  end
  m2 = GSL::Matrix.alloc(n, n)
  GSL::Matrix.memcpy(m2, m1)
  mdelta = GSL::Matrix.alloc(n, n)
  work = GSL::Wavelet::Workspace.alloc(n)
  w = GSL::Wavelet.alloc(t, member)
  case type
  when 1
    GSL::Wavelet2d::transform_matrix_forward(w, m2, work)
    GSL::Wavelet2d::transform_matrix_inverse(w, m2, work)
  when 2
    GSL::Wavelet2d::nstransform_matrix_forward(w, m2, work)
    GSL::Wavelet2d::nstransform_matrix_inverse(w, m2, work)
  end

  for i in 0...n
    for j in 0...n
      x1 = m1[i][j]
      x2 = m2[i][j]
      mdelta.set(i, j, (x1 - x2).abs)
    end
  end
  i, j = mdelta.max_index
  x1 = m1[i][j]
  x2 = m1[i][j]
  GSL::Test::test((x2 - x1).abs > n*1e-15 ? 1 : 0, "#{w.name}(#{member})-2d #{typename}, n = #{n}, tda = #{tda}, maxerr = #{(x2 - x1).abs}")
  if tda > n
    status = 0
    for i in 0...n
      for j in n...tda
        status |= (data[i*tda+j] != (12345.0 + (i*tda+j))) ? 1 : 0
      end
    end
    GSL::Test::test(status, "#{w.name}(#{member})-2d #{typename} other data untouched, n = #{n}, tda = #{tda}")
  end
end

Member = [309, 307, 305, 303, 301, 208, 206, 204, 202, 105, 103]
S = 1
NS = 2

n = 1
while n <= 16384
  for stride in 1..1
    for i in 0...N_BS
      test_1d(n, stride, "bspline", Member[i]);
      test_1d(n, stride, "bspline_centered", Member[i]);
    end
    i = 4
    while i <= 20
      test_1d(n, stride, "daubechies", i)
      test_1d(n, stride, "daubechies_centered", i)
      i += 2
    end
      test_1d(n, stride, "haar", 2)
      test_1d(n, stride, "haar_centered", 2)
  end
  n *= 2
end

n = 1
while n <= 64
  for tda in n..n
    for i in 0...N_BS
      test_2d(n, tda, "bspline", Member[i], S);
      test_2d(n, tda, "bspline_centered", Member[i], S);
      
      test_2d(n, tda, "bspline", Member[i], NS);
      test_2d(n, tda, "bspline_centered", Member[i], NS);
    end
  end
  n *= 2
end
          


