#!/usr/bin/env ruby

require 'rubygems'
require 'narray'
require 'gsl'
require '../gsl_test.rb'
include GSL::Test

dbleps = 1e-6
expected = Math.sqrt((0..4).inject {|m,x| m+=x*x})

v = GSL::Vector.indgen(5)
v_dnrm2 = GSL::Blas.dnrm2(v)
GSL::Test.test_rel(v_dnrm2, expected, dbleps, "GSL::Blas.dnrm2(GSL::Vector)")

na = NArray.float(5).indgen!
na_dnrm2 = GSL::Blas.dnrm2(na)
GSL::Test.test_rel(na_dnrm2, expected, dbleps, "GSL::Blas.dnrm2(NArray)")

GSL::Test.test_rel(na_dnrm2, v_dnrm2, 0, "GSL::Blas.dnrm2(NArray) == GSL::Blas.dnrm2(GSL::Vector)")
