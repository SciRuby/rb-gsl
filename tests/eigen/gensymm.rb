#!/usr/bin/env ruby
require("gsl")
require("../gsl_test.rb")
require("./eigen.rb")
include GSL::Test

def test_eigen_gensymm_results(a, b, eval, evec, count, desc, desc2)
  n = a.size1
  # check A v = lambda B v 

  for i in 0...n do
    ei = eval[i]
    vi = evec.column(i)
    norm = vi.dnrm2
    # check that eigenvector is normalized 
    GSL::Test::test_rel(norm, 1.0, n * GSL::DBL_EPSILON,
                        "gensymm(N=#{n},cnt=#{count}), #{desc}, normalized(#{i}), #{desc2}")
    
    y = a*vi
    x = b*vi
    x *= ei
    
    for j in 0...n do
      GSL::Test::test_rel(y[j], x[j], 1e9 * GSL::DBL_EPSILON, 
                          "gensymm(N=#{n},cnt=#{count}), #{desc}, eigenvalue(#{i},#{j}), #{desc2}")
    end
  end
end

def test_eigen_gensymm()
  n_max = 20
  rng = GSL::Rng.alloc()
  for n in 1..n_max do
    w  = GSL::Eigen::Gensymm::Workspace.alloc(n)
    wv = GSL::Eigen::Gensymmv::Workspace.alloc(n)		
    for i in 0...5 do
      a = create_random_symm_matrix(n, n, rng, -10, 10)
      b = create_random_posdef_matrix(n, n, rng)
      
      aa = a.clone
      bb = b.clone
      evalv, evec = GSL::Eigen::gensymmv(aa, bb, wv)
      test_eigen_gensymm_results(a, b, evalv, evec, i, "random", "unsorted")
      
      aa = a.clone
      bb = b.clone
      eval = GSL::Eigen::gensymm(aa, bb, w)
      
      x = eval.sort
      y = evalv.sort
      
      test_eigenvalues_real(y, x, "gensymm, random", "unsorted");
      
      GSL::Eigen::Gensymmv::sort(evalv, evec, GSL::EIGEN_SORT_VAL_ASC);
      test_eigen_gensymm_results(a, b, evalv, evec, i, "random", "val/asc");
      
      GSL::Eigen::Gensymmv::sort(evalv, evec, GSL::EIGEN_SORT_VAL_DESC);
      test_eigen_gensymm_results(a, b, evalv, evec, i, "random", "val/desc");
      
      GSL::Eigen::Gensymmv::sort(evalv, evec, GSL::EIGEN_SORT_ABS_ASC);
      test_eigen_gensymm_results(a, b, evalv, evec, i, "random", "abs/asc");
      GSL::Eigen::Gensymmv::sort(evalv, evec, GSL::EIGEN_SORT_ABS_DESC);
      test_eigen_gensymm_results(a, b, evalv, evec, i, "random", "abs/desc");
    end
  end
end

test_eigen_gensymm()
