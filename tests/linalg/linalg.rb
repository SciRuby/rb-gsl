#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
include GSL::Test

def check(x, actual, eps)
  if x == actual
    return 0
  elsif actual.zero?
    return x.abs > eps ? 1 : 0
  else
    return ((x - actual).abs/actual.abs > eps) ? 1 : 0
  end
end

def create_hilbert_matrix(size)
  m = GSL::Matrix.alloc(size, size)
  for i in 0...size
    for j in 0...size
      m.set(i, j, 1.0/(i+j+1.0))
    end
  end
  return m
end

def create_general_matrix(size1, size2)
  m = GSL::Matrix.alloc(size1, size2)
  for i in 0...size1
    for j in 0...size2
      m.set(i, j, 1.0/(i+j+1.0))
    end
  end
  return m
end

def create_singular_matrix(size1, size2)
  m = create_general_matrix(size1, size2)
  for j in 0...m.size2
    m.set(0, j, 0.0)
  end
  return m
end

def create_vandermonde_matrix(size)
  m = GSL::Matrix.alloc(size, size)
  for i in 0...size
    for j in 0...size
      m.set(i, j, pow(i + 1.0, size - j - 1.0))
    end
  end
  return m
end

def create_moler_matrix(size)
  m = GSL::Matrix.alloc(size, size)
  for i in 0...size
    for j in 0...size
      m.set(i, j, GSL::MIN(i+1, j+1) - 2.0)
    end
  end
  return m
end

def create_complex_matrix(size)
  m = GSL::Matrix::Complex.alloc(size, size)
  for i in 0...size
    for j in 0...size
      z = GSL::Complex.alloc(1.0/(i+j+1.0), 1/(i*i+j*j+0.5))
      m.set(i, j, z)
    end
  end
  return m
end

def create_row_matrix(size1, size2)
  m = GSL::Matrix.alloc(size1, size2)
  for i in 0...size1
    m.set(i, 0, 1.0/(i + 1.0))
  end
  return m
end

def create_2x2_matrix(a11, a12, a21, a22)
  return GSL::Matrix.alloc(a11, a12, a21, a22)
end

def create_diagonal_matrix(a, size)
  m = GSL::Matrix.alloc(size, size)
  for i in 0...size
    m.set(i, i, a[i])
  end
  return m
end

Hilb2 = create_hilbert_matrix(2)
Hilb3 = create_hilbert_matrix(3)
Hilb4 = create_hilbert_matrix(4)
Hilb12 = create_hilbert_matrix(12)
Vander2  = create_vandermonde_matrix(2)
Vander3  = create_vandermonde_matrix(3)
Vander4  = create_vandermonde_matrix(4)
Vander12  = create_vandermonde_matrix(12)

_inf5_data = GSL::Vector.alloc(1.0, 0.0, -3.0, 0.0, -5.0)
_m53_lssolution = GSL::Vector.alloc(52.5992295702070, -337.7263113752073, 351.8823436427604)
Hilb2_solution = GSL::Vector.alloc(-8.0, 18.0)
Hilb3_solution = GSL::Vector.alloc(27.0, -192.0, 210.0)
Hilb4_solution = GSL::Vector.alloc(-64.0, 900.0, -2520.0, 1820.0)
Hilb12_solution = GSL::Vector.alloc(-1728.0, 245388.0, -8528520.0, 
                             127026900.0, -1009008000.0, 4768571808.0, 
                             -14202796608.0, 27336497760.0, -33921201600.0,
                             26189163000.0, -11437874448.0, 2157916488.0)

_c7_solution = GSL::Vector.alloc(2.40717272023734e+01, -9.84612797621247e+00,
                         -2.69338853034031e+02, 8.75455232472528e+01,
                         2.96661356736296e+03, -1.02624473923993e+03,
                         -1.82073812124749e+04, 5.67384473042410e+03,
                         5.57693879019068e+04, -1.61540963210502e+04,
                         -7.88941207561151e+04, 1.95053812987858e+04,
                         3.95548551241728e+04, -7.76593696255317e+03)

Vander2_solution = GSL::Vector.alloc(1.0, 0.0) 
Vander3_solution = GSL::Vector.alloc(0.0, 1.0, 0.0) 
Vander4_solution = GSL::Vector.alloc(0.0, 0.0, 1.0, 0.0) 
Vander12_solution = GSL::Vector.alloc(0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0, 
                            0.0, 0.0, 1.0, 0.0) 

def test_matmult()
  s = 0
  a = GSL::Matrix.alloc([10.0, 5.0, 1.0, 20.0], 2, 2)
  b = GSL::Matrix.alloc([10.0, 5.0, 2.0, 1.0, 3.0, 2.0], 2, 3)
  c = a**b
  s += ((c[0,0] - 105.0).abs > GSL::DBL_EPSILON) ? 1 : 0
  s += ((c[0,1] -  65.0).abs > GSL::DBL_EPSILON) ? 1 : 0
  s += ((c[0,2] -  30.0).abs > GSL::DBL_EPSILON) ? 1 : 0
  s += ((c[1,0] -  30.0).abs > GSL::DBL_EPSILON) ? 1 : 0
  s += ((c[1,1] -  65.0).abs > GSL::DBL_EPSILON) ? 1 : 0
  s += ((c[1,2] -  42.0).abs > GSL::DBL_EPSILON) ? 1 : 0
  return s
end

def test_matmult_mod()
  GSL::Matrix[[10.0, 5.0, 1.0], [1.0, 20.0, 5.0], [1.0, 3.0, 7.0]]
  GSL::Matrix[[10.0, 5.0, 2.0], [1, 3, 2], [1, 3, 2]]
  GSL::Matrix[[10, 5, 1], [1, 20, 5]]
  
end

def test_LU_solve_dim(m, actual, eps)
  dim = m.size1
  _perm = Permutation.alloc(dim)
  _rhs = GSL::Vector[1..dim]
  _lu = m.clone

end
