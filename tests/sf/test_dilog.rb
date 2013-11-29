#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
require("./gsl_test_sf.rb")
include GSL::Test
include Math
include GSL::Test::Sf
s = 0
_m = GSL::MODE_DEFAULT
  TEST_SF(s, "GSL::Sf::dilog_e", "(-3.0)",   -1.9393754207667089531,     TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(-0.5)",   -0.4484142069236462024,     TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(-0.001)", -0.0009997501110486510834,  TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(0.1)",     0.1026177910993911,        TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(0.7)",     0.8893776242860387386,     TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(1.0)",     1.6449340668482260,        TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(1.5)",     2.3743952702724802007,     TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(2.0)",     2.4674011002723397,        TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "( 5.0)",    1.7837191612666306277,     TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "( 11.0)",   0.3218540439999117111,     TEST_TOL1, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(12.59)",   0.0010060918167266208634,  TEST_TOL3, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(12.595)",  0.00003314826006436236810, TEST_TOL5, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(13.0)",   -0.07806971248458575855,    TEST_TOL2, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(20.0)",   -1.2479770861745251168,     TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(150.0)",  -9.270042702348657270,      TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::dilog_e", "(1100.0)", -21.232504073931749553,     TEST_TOL0, GSL::SUCCESS)
