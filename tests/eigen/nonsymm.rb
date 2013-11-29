#!/usr/bin/env ruby
require("gsl")
require("../gsl_test.rb")
include GSL::Test

m = GSL::Matrix[[1, 2], [3, 2]]
p m
eval = m.eigen_nonsymm.real.sort
GSL::Test::test_abs(eval[0], -1, 0, "GSL::Matrix::eigen_nonsymm")
GSL::Test::test_abs(eval[1], 4, 0, "GSL::Matrix::eigen_nonsymm")

m = GSL::Matrix[[1, 4, 2, 3], 2, 2]
p m
eval = m.eigen_nonsymm.real.sort
GSL::Test::test_abs(eval[0], -1, 0, "GSL::Matrix::eigen_nonsymm")
GSL::Test::test_abs(eval[1], 5, 0, "GSL::Matrix::eigen_nonsymm")

m = GSL::Matrix[[2, 4], [3, 1]]
p m
eval = m.eigen_nonsymm.real.sort
GSL::Test::test_abs(eval[0], -2, 0, "GSL::Matrix::eigen_nonsymm")
GSL::Test::test_abs(eval[1], 5, 0, "GSL::Matrix::eigen_nonsymm")

m = GSL::Matrix[[5, 6, 3, 2], 2, 2]
p m
eval = m.eigen_nonsymm.real.sort
GSL::Test::test_abs(eval[0], -1, 0, "GSL::Matrix::eigen_nonsymm")
GSL::Test::test_abs(eval[1], 8, 0, "GSL::Matrix::eigen_nonsymm")

m = GSL::Matrix[[4, 1, -1, 2, 5, -2, 1, 1, 2], 3, 3]
p m
eval = m.eigen_nonsymm.real.sort
GSL::Test::test_abs(eval[0], 3, 1e-10, "GSL::Matrix::eigen_nonsymm")
GSL::Test::test_abs(eval[1], 3, 0, "GSL::Matrix::eigen_nonsymm")
# 2008/Oct/17 YT
# This test fails in Darwin9.5.0-gcc4.0.1
#   expected: 5
#   obtained: 4.99999999999999911
GSL::Test::test_abs(eval[2], 5, 0, "GSL::Matrix::eigen_nonsymm")

m = GSL::Matrix[[-3, 1, -1], [-7, 5, -1], [-6, 6, -2]]
p m
eval = m.eigen_nonsymm.real.sort
GSL::Test::test_abs(eval[0], -2, 1e-6, "GSL::Matrix::eigen_nonsymm")
GSL::Test::test_abs(eval[1], -2, 1e-6, "GSL::Matrix::eigen_nonsymm")
GSL::Test::test_abs(eval[2], 4, 1e-10, "GSL::Matrix::eigen_nonsymm")

m = GSL::Matrix[[11, -8, 4, -8, -1, -2, 4, -2, -4], 3, 3]
p m
eval = m.eigen_nonsymm.real.sort
GSL::Test::test_abs(eval[0], -5, 1e-10, "GSL::Matrix::eigen_nonsymm")
GSL::Test::test_abs(eval[1], -5, 1e-10, "GSL::Matrix::eigen_nonsymm")
GSL::Test::test_abs(eval[2], 16, 1e-10, "GSL::Matrix::eigen_nonsymm")
