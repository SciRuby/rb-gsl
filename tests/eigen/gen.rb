#!/usr/bin/env ruby
require("gsl")
require("../gsl_test.rb")
require("./eigen.rb")
include GSL::Test

def test_eigen_gen_results(a, b, alpha, beta, evec, count, desc, desc2)
  n = a.size1
  ma = a.to_complex
  mb = b.to_complex
  
  for i in 0...n do
    vi = evec.column(i)
    ai = alpha[i]
    bi = beta[i]
    
    x = mb*vi
    x *= ai
    
    y = ma*vi
    y *= bi
    
    for j in 0...n do
      GSL::Test::test_abs(y[j].real, x[j].real, 1e8*GSL::DBL_EPSILON, 
                          "gen(N=#{n},cnt=#{count}), #{desc}, eigenvalue(#{i},#{j}), real, #{desc2}")
      GSL::Test::test_abs(y[j].imag, x[j].imag, 1e8*GSL::DBL_EPSILON, 
                          "gen(N=#{n},cnt=#{count}), #{desc}, eigenvalue(#{i},#{j}), imag, #{desc2}")
    end
  end
end

def test_eigen_gen_pencil(a, b, count, desc, test_schur, w, wv)
  n = a.size1
  
  aa = a.clone
  bb = b.clone
  
  if test_schur == 1
    alphav, betav, evec, q, z = GSL::Eigen::genv_QZ(aa, bb, wv)
    test_eigen_schur(a, aa, q, z, count, "genv/A", desc)
    test_eigen_schur(b, bb, q, z, count, "genv/B", desc)
  else
    alphav, betav, evec = GSL::Eigen::genv(aa, bb, wv)
  end
  
  test_eigen_gen_results(a, b, alphav, betav, evec, count, desc, "unsorted")
  
  aa = a.clone
  bb = b.clone
  if test_schur == 1
    GSL::Eigen::gen_params(1, 1, 0, w)
    alpha, beta, q, z = GSL::Eigen::gen_QZ(aa, bb, w)
    test_eigen_schur(a, aa, q, z, count, "gen/A", desc)
    test_eigen_schur(b, bb, q, z, count, "gen/B", desc)
  else
    GSL::Eigen::gen_params(0, 0, 0, w)
    alpha, beta = GSL::Eigen::gen(aa, bb, w)
  end
  
  eval = GSL::Vector::Complex.alloc(n)
  evalv = GSL::Vector::Complex.alloc(n)	
  for i in 0...n do
    ai = alpha[i]
    bi = beta[i]
    z = GSL::Complex.alloc(ai.real/bi, ai.imag/bi)
    eval[i] = z
    
    ai = alphav[i]
    bi = betav[i]
    z = GSL::Complex.alloc(ai.real/bi, ai.imag/bi)
    evalv[i] = z
  end
  
  GSL::Eigen::nonsymmv_sort(eval, nil, GSL::EIGEN_SORT_ABS_ASC)
  GSL::Eigen::nonsymmv_sort(evalv, nil, GSL::EIGEN_SORT_ABS_ASC)
  test_eigenvalues_complex(evalv, eval, "gen", desc)
  
  GSL::Eigen::genv_sort(alphav, betav, evec, GSL::EIGEN_SORT_ABS_ASC)
  test_eigen_gen_results(a, b, alphav, betav, evec, count, desc, "abs/asc")
  GSL::Eigen::genv_sort(alphav, betav, evec, GSL::EIGEN_SORT_ABS_DESC)
  test_eigen_gen_results(a, b, alphav, betav, evec, count, desc, "abs/desc")
end

def test_eigen_gen()
  n_max = 20
  rng = GSL::Rng.alloc()
  for n in 1..n_max do
    w = GSL::Eigen::Gen.alloc(n)
    wv = GSL::Eigen::Genv.alloc(n)		
    for i in 0...5 do
      a = create_random_nonsymm_matrix(n, n, rng, -10, 10)
      b = create_random_nonsymm_matrix(n, n, rng, -10, 10)
      test_eigen_gen_pencil(a, b, i, "random", 0, w, wv);
      test_eigen_gen_pencil(a, b, i, "random", 1, w, wv);
    end
  end
  ma = GSL::Matrix.alloc([1, 1, 0, 0, 0, -1, 1, 0, 0], 3, 3)
  mb = GSL::Matrix.alloc([-1, 0, -1, 0, -1, 0, 0, 0, -1], 3, 3)
  w = GSL::Eigen::Gen.alloc(3)
  wv = GSL::Eigen::Genv.alloc(3)	    
  test_eigen_gen_pencil(ma, mb, 0, "integer", 0, w, wv);
  test_eigen_gen_pencil(ma, mb, 0, "integer", 1, w, wv);
end

test_eigen_gen()
