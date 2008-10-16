#!/usr/bin/env ruby
require("gsl")
require("../gsl_test.rb")
include GSL::Test

v = GSL::Vector.alloc(0.537, 0.826)
expected = 1
k = v.idamax
GSL::Test::test_int(k, expected, "damax")

vz = GSL::Vector::Complex.alloc([0.913, -0.436], [-0.134, 0.129])
expected = 0
k = vz.izamax
GSL::Test::test_int(k, expected, "zmax")
