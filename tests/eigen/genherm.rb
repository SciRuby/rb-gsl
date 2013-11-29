#!/usr/bin/env ruby
require("gsl")
require("../gsl_test.rb")
require("./eigen.rb")
include GSL::Test

def test_eigen_genherm_results(a, b, eval, evec, count, desc, desc2)
	n = a.size1

  # check A v = lambda B v 
  for i in 0...n do
		ei = eval[i]
		vi = evec.column(i)
		norm = vi.nrm2
    # check that eigenvector is normalized 
    GSL::Test::test_rel(norm, 1.0, n * GSL::DBL_EPSILON,
                   "genherm(N=#{n},cnt=#{count}), #{desc}, normalized(#{i}), #{desc2}")

    y = a*vi
    x = b*vi
		x *= ei

		for j in 0...n do
      GSL::Test::test_rel(y[j].real, x[j].real, 1e9 * GSL::DBL_EPSILON, 
                       "genherm(N=#{n},cnt=#{count}), #{desc}, eigenvalue(#{i},#{j}), real, #{desc2}")
	    GSL::Test::test_rel(y[j].imag, x[j].imag, 1e9 * GSL::DBL_EPSILON, 
                       "genherm(N=#{n},cnt=#{count}), #{desc}, eigenvalue(#{i},#{j}), imag, #{desc2}")                       
		end
  end
end

def test_eigen_genherm()
  n_max = 20
	rng = GSL::Rng.alloc()
	for n in 1..n_max do
		w = GSL::Eigen::Genherm.alloc(n)
		wv = GSL::Eigen::Genhermv.alloc(n)		
		for i in 0...5 do
			a = create_random_herm_matrix(n, n, rng, -10, 10)
      b = create_random_complex_posdef_matrix(n, n, rng)

			evalv, evec = GSL::Eigen::genhermv(a, b, wv)
			test_eigen_genherm_results(a, b, evalv, evec, i, "random", "unsorted")
	
			eval = GSL::Eigen::genherm(a, b, w)

			x = eval.sort
			y = evalv.sort

      test_eigenvalues_real(y, x, "genherm, random", "unsorted")
			
      GSL::Eigen::genhermv_sort(evalv, evec, GSL::EIGEN_SORT_VAL_ASC);
      test_eigen_genherm_results(a, b, evalv, evec, i, "random", "val/asc");

      GSL::Eigen::genhermv_sort(evalv, evec, GSL::EIGEN_SORT_VAL_DESC);
      test_eigen_genherm_results(a, b, evalv, evec, i, "random", "val/desc");

      GSL::Eigen::genhermv_sort(evalv, evec, GSL::EIGEN_SORT_ABS_ASC);
      test_eigen_genherm_results(a, b, evalv, evec, i, "random", "abs/asc");
      GSL::Eigen::genhermv_sort(evalv, evec, GSL::EIGEN_SORT_ABS_DESC);
      test_eigen_genherm_results(a, b, evalv, evec, i, "random", "abs/desc");
		end
	end
end

test_eigen_genherm()
