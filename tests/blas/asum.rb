#!/usr/bin/env ruby
require("gsl")
require("../gsl_test.rb")
include GSL::Test

dbleps = 1e-6

v = GSL::Vector.alloc(0.271, -0.012)
expected = 0.283
f = v.dasum
GSL::Test::test_rel(f, expected, dbleps, "dasum")

vz = GSL::Vector::Complex.alloc([-0.046, -0.671], [-0.323, 0.785])
expected = 1.825
f = vz.dzasum
GSL::Test::test_rel(f, expected, dbleps, "dzasum")
