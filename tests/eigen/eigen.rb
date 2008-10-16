#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "eigen/test.c"
require("gsl")
require("../gsl_test.rb")
include GSL::Test

def create_random_symm_matrix(size1, size2, rng, lower, upper)
	m = GSL::Matrix.alloc(size1, size2)
  for i in 0...size1 do
  	for j in i...size2 do
  		x = rng.uniform()*(upper-lower) + lower
  		m[i][j] = x
  		m[j][i] = x
		end
	end
	m
end

def create_random_herm_matrix(size1, size2, rng, lower, upper)
	m = GSL::Matrix::Complex.alloc(size1, size2)
	for i in 0...size1 do
		for j in i...size2 do
			re = rng.uniform()* (upper - lower) + lower
      if i == j
      	im = 0.0
      else
				im = rng.uniform()* (upper - lower) + lower
			end
			z = GSL::Complex.alloc(re, im)
			m[i][j] = z
			m[j][i] = z.conjugate
		end
	end
	m
end

def create_random_posdef_matrix(size1, size2, rng)
	m = GSL::Matrix.alloc(size1, size2)
	x = rng.uniform()
  for i in 0...size1 do
  	for j in i...size2 do
        a = pow(x, (j - i).to_f)
				m[i][j] = a
				m[j][i] = a
		end
	end
	m
end

def create_random_complex_posdef_matrix(size1, size2, rng)
	m = GSL::Matrix::Complex.calloc(size1, size2)
	n = size1

	for i in 0...n do
		x = rng.uniform()
		z = GSL::Complex.alloc(x, 0.0)
		m[i][i] = z
	end
	
	work = GSL::Vector::Complex.alloc(n)
	for i in 0...n do
		for j in 0...n do
			x = 2.0*rng.uniform() - 1.0
			y = 2.0*rng.uniform() - 1.0			
			z = GSL::Complex.alloc(x, y)
			work[j] 
		end
		tau = GSL::Linalg::Complex::householder_transform(work)
		GSL::Linalg::Complex::householder_hm(tau, work, m)
		GSL::Linalg::Complex::householder_mh(tau.conjugate, work, m)		
	end
	m
end

def create_random_nonsymm_matrix(size1, size2, rng, lower, upper)
	m = GSL::Matrix.alloc(size1, size2)
  for i in 0...size1 do
  	for j in 0...size2 do
  		m[i][j] = rng.uniform()*(upper - lower) + lower
  	end
  end
  m
end

def test_eigen_results(n, m, eval, evec, desc, desc2)
  x = GSL::Vector.alloc(n)
  y = GSL::Vector.alloc(n)

# check eigenvalues
  for i in 0...n
    ei = eval[i]
    vi = evec.col(i)
    GSL::Vector.memcpy(x, vi)
    y = GSL::Blas.dgemv(GSL::Blas::NoTrans, 1.0, m, x, 0.0, y)
#    y = GSL::Blas.dgemv(GSL::Blas::NoTrans, 1.0, m, x)
    for j in 0...n
      xj = x[j]
      yj = y[j]
      desc3 = sprintf("%s, eigenvalue(%d,%d), %s", desc, i, j, desc2)
      GSL::Test::test_rel(yj, ei*xj, 1e8*DBL_EPSILON, desc3)
    end
  end
# check eigenvectors are orthonormal
  for i in 0...n
    vi = evec.col(i)
    nrm_v = GSL::Blas.dnrm2(vi)
    desc3 = sprintf("%s, normalized(%d), %s", desc, i, desc2)
    GSL::Test::test_rel(nrm_v, 1.0, n*DBL_EPSILON, desc3)
  end

  for i in 0...n
    vi = evec.col(i)
    for j in (i+1)...n
      vj = evec.col(j)
      vivj = GSL::Blas.ddot(vi, vj)
      desc3 = sprintf("%s, orthogonal(%d,%d), %s", desc, i, j, desc2)
      GSL::Test::test_abs(vivj, 0.0, n*DBL_EPSILON, desc3)
    end
  end
end

def test_eigen_complex_results(n, m, eval, evec, desc, desc2)
  x = GSL::Vector::Complex.alloc(n)
  y = GSL::Vector::Complex.alloc(n)

# check eigenvalues
  for i in 0...n
    ei = eval[i]
    vi = evec.col(i)
    GSL::Vector::Complex.memcpy(x, vi)
    c1 = GSL::Complex.alloc(1.0, 0.0)
    c0 = GSL::Complex.alloc(0.0, 0.0)
    y = GSL::Blas.zgemv(GSL::Blas::NoTrans, c1, m, x, c0, y)
    for j in 0...n
      xj = x[j]
      yj = y[j]
      desc3 = sprintf("%s, eigenvalue(%d,%d), real, %s", desc, i, j, desc2)
      GSL::Test::test_rel(yj.re, ei*xj.re, 1e8*DBL_EPSILON, desc3)
      desc3 = sprintf("%s, eigenvalue(%d,%d), imag, %s", desc, i, j, desc2)
      GSL::Test::test_rel(yj.im, ei*xj.im, 1e8*DBL_EPSILON, desc3)
    end
  end
# check eigenvectors are orthonormal
  for i in 0...n
    vi = evec.col(i)
    nrm_v = GSL::Blas.dznrm2(vi)
    desc3 = sprintf("%s, normalized(%d), %s", desc, i, desc2)
    GSL::Test::test_rel(nrm_v, 1.0, n*DBL_EPSILON, desc3)
  end

  for i in 0...n
    vi = evec.col(i)
    for j in (i+1)...n
      vj = evec.col(j)
      vivj = GSL::Blas.zdotc(vi, vj)
      desc3 = sprintf("%s, orthogonal(%d,%d), %s", desc, i, j, desc2)
      GSL::Test::test_abs(vivj.abs, 0.0, n*DBL_EPSILON, desc3)
    end
  end
end

def test_eigenvalues(n, eval, eval2, desc, desc2)
  for i in 0...n
    ei = eval[i]
    e2i = eval2[i]
    desc3 = sprintf("%s, direct eigenvalue(%d), %s", desc, i, desc2)
    GSL::Test::test_rel(ei, e2i, GSL::DBL_EPSILON, desc3)
  end
end

def chop_subnormals(x)
  # Chop any subnormal values */
  return x.abs < GSL::DBL_MIN ? 0 : x
end

def test_eigenvalues_real(eval, eval2, desc, desc2)
	n = eval.size
	emax = 0
  # check eigenvalues 
  for i in 0...n do
		ei = eval[i]
		if ei.abs > emax
			emax = ei.abs
		end
	end

  for i in 0...n do
  	ei = eval[i]
  	e2i = chop_subnormals(eval2[i])
		GSL::Test::test_abs(ei, e2i, emax * 1e8 * GSL::DBL_EPSILON,  "#{desc}, direct eigenvalue(#{i}), #{desc2}")
  end
end

def test_eigenvalues_complex(eval, eval2, desc, desc2)
	n = eval.size
	for i in 0...n do
		ei = eval[i]
		e2i = eval2[i]
		GSL::Test::test_rel(ei.real, e2i.real, 10*n*GSL::DBL_EPSILON, 
                   "#{desc}, direct eigenvalue(#{i}) real, #{desc2}")

		GSL::Test::test_rel(ei.imag, e2i.imag, 10*n*GSL::DBL_EPSILON, 
                   "#{desc}, direct eigenvalue(#{i}) imag, #{desc2}")		
	end
end

def test_eigen_schur(a, s, q, z, count, desc, desc2)
	n = a.size1

	t1 = a*z
	t2 = q*s
	for i in 0...n do
		for j in 0...n do
			x = t1[i][j]
			y = t2[i][j]
	    GSL::test_abs(x, y, 1.0e8 * GSL::DBL_EPSILON,
                       "#{desc}(N=#{n},cnt=#{count}), #{desc2}, schur(#{i},#{j})")
  	end
 	end
end
