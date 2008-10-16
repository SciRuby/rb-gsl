#!/usr/bin/env ruby
require("gsl")
require("../gsl_test.rb")
include GSL::Test

def test_nonsymmv2(m, eps)
  p m
  m2 = m.clone
  eval, evec = m2.eigen_nonsymmv
  evalre = eval.real
  evecre = evec.real
  a = evecre.inv*m*evecre
  GSL::Test::test_abs(a[0][0], evalre[0], eps, "GSL::Matrix::eigen_nonsymmv")
  GSL::Test::test_abs(a[1][1], evalre[1], eps, "GSL::Matrix::eigen_nonsymmv")
  GSL::Test::test_abs(a[0][1], 0, eps, "GSL::Matrix::eigen_nonsymmv")
  GSL::Test::test_abs(a[1][0], 0, eps, "GSL::Matrix::eigen_nonsymmv")
end

def test_nonsymmv3(m, eps)
  p m
  m2 = m.clone
  eval, evec = m2.eigen_nonsymmv
  evalre = eval.real
  evecre = evec.real
  a = evecre.inv*m*evecre
  GSL::Test::test_abs(a[0][0], evalre[0], eps, "GSL::Matrix::eigen_nonsymmv")
  GSL::Test::test_abs(a[1][1], evalre[1], eps, "GSL::Matrix::eigen_nonsymmv")
  GSL::Test::test_abs(a[2][2], evalre[2], eps, "GSL::Matrix::eigen_nonsymmv")
  GSL::Test::test_abs(a[0][1], 0, eps, "GSL::Matrix::eigen_nonsymmv")
  GSL::Test::test_abs(a[0][2], 0, eps, "GSL::Matrix::eigen_nonsymmv")
  GSL::Test::test_abs(a[1][0], 0, eps, "GSL::Matrix::eigen_nonsymmv")
  GSL::Test::test_abs(a[1][2], 0, eps, "GSL::Matrix::eigen_nonsymmv")
  GSL::Test::test_abs(a[2][0], 0, eps, "GSL::Matrix::eigen_nonsymmv")
  GSL::Test::test_abs(a[2][1], 0, eps, "GSL::Matrix::eigen_nonsymmv")
end

m = GSL::Matrix[[1, 2], [3, 2]]
test_nonsymmv2(m, 1e-10)

m = GSL::Matrix[[4, 2, 3, -1], 2, 2]
test_nonsymmv2(m, 1e-10)

m = GSL::Matrix[[2, 2, -5], [3, 7, -15], [1, 2, -4]]
test_nonsymmv3(m, 1e-10)

m = GSL::Matrix[[4, 1, -1], [2, 5, -2], [1, 1, 2]]
test_nonsymmv3(m, 1e-10)

m = GSL::Matrix[[-3, 1, -1], [-7, 5, -1], [-6, 6, -2]]
test_nonsymmv3(m, 1e-6)

m = GSL::Matrix[[11, -8, 4], [-8, -1, -2], [4, -2, -4]]
test_nonsymmv3(m, 1e-10)
