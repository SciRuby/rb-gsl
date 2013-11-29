#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "eigen/test.c"
require("gsl")
require("../gsl_test.rb")
require("./eigen.rb")
include GSL::Test

def test_eigen_symm(desc, m)
  n = m.size1
  a = GSL::Matrix.alloc(n, n)

  w1 = GSL::Eigen::Symm::Workspace.alloc(n)
  w2 = GSL::Eigen::Symmv::Workspace.alloc(n)

  GSL::Matrix.memcpy(a, m)
  eval, evec = a.eigen_symmv(w2)
  test_eigen_results(n, m, eval, evec, desc, "unsorted")

  GSL::Matrix.memcpy(a, m)
  eval2 = a.eigen_symm(w1)
  test_eigenvalues(n, eval, eval2, desc, "unsorted")

  GSL::Eigen::Symmv::sort(eval, evec, GSL::Eigen::SORT_VAL_ASC)
  test_eigen_results(n, m, eval, evec, desc, "val/asc")
  GSL::Eigen::Symmv::sort(eval, evec, GSL::Eigen::SORT_VAL_DESC)
  test_eigen_results(n, m, eval, evec, desc, "val/desc")
  GSL::Eigen::Symmv::sort(eval, evec, GSL::Eigen::SORT_ABS_ASC)
  test_eigen_results(n, m, eval, evec, desc, "abs/asc")
  GSL::Eigen::Symmv::sort(eval, evec, GSL::Eigen::SORT_ABS_DESC)
  test_eigen_results(n, m, eval, evec, desc, "abs/desc")
end

def test_eigen_herm(desc, m)
  n = m.size1
  a = GSL::Matrix::Complex.alloc(n, n)

  w1 = GSL::Eigen::Herm::Workspace.alloc(n)
  w2 = GSL::Eigen::Hermv::Workspace.alloc(n)

  GSL::Matrix::Complex.memcpy(a, m)
  eval, evec = a.eigen_hermv(w2)
  test_eigen_complex_results(n, m, eval, evec, desc, "unsorted")

  GSL::Matrix::Complex.memcpy(a, m)
  eval2 = a.eigen_herm(w1)
  test_eigenvalues(n, eval, eval2, desc, "unsorted")

  GSL::Eigen::Hermv::sort(eval, evec, GSL::Eigen::SORT_VAL_ASC)
  test_eigen_complex_results(n, m, eval, evec, desc, "val/asc")
  GSL::Eigen::Hermv::sort(eval, evec, GSL::Eigen::SORT_VAL_DESC)
  test_eigen_complex_results(n, m, eval, evec, desc, "val/desc")
  GSL::Eigen::Hermv::sort(eval, evec, GSL::Eigen::SORT_ABS_ASC)
  test_eigen_complex_results(n, m, eval, evec, desc, "abs/asc")
  GSL::Eigen::Hermv::sort(eval, evec, GSL::Eigen::SORT_ABS_DESC)
  test_eigen_complex_results(n, m, eval, evec, desc, "abs/desc")
end



r = GSL::Matrix.alloc([0, 0, -1, 0], [0, 1, 0, 1], [-1, 0, 0, 0], [0, 1, 0, 0])
test_eigen_symm("symm(4)", r)

c = r.to_complex
test_eigen_herm("herm(4)", c)

r = GSL::Matrix.alloc(4, 4)
r.set_diagonal([1, 2, 3, 4])
test_eigen_symm("symm(4) diag", r)

c = r.to_complex
test_eigen_herm("herm(4) diag", c)



