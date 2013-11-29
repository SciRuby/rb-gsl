#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
require("./gsl_test_sf.rb")
include GSL::Test
include Math
include GSL::Test::Sf
s = 0
_m = GSL::MODE_DEFAULT
  TEST_SF(s, "GSL::Sf::hydrogenicR_1_e", "(3.0, 2.0)",  0.025759948256148471036,  TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::hydrogenicR_1_e", "(3.0, 10.0)", 9.724727052062819704e-13, TEST_TOL1, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::hydrogenicR_e", "(4, 0, 3.0, 2.0)", -0.03623182256981820062,  TEST_TOL2, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::hydrogenicR_e", "(4, 1, 3.0, 2.0)", -0.028065049083129581005, TEST_TOL2, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::hydrogenicR_e", "(4, 2, 3.0, 2.0)",  0.14583027278668431009,  TEST_TOL0, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::hydrogenicR_e", "(100,  0, 3.0, 2.0)", -0.00007938950980052281367, TEST_TOL3, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::hydrogenicR_e", "(100, 10, 3.0, 2.0)",  7.112823375353605977e-12,  TEST_TOL2, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::hydrogenicR_e", "(100, 90, 3.0, 2.0)",  5.845231751418131548e-245, TEST_TOL2, GSL::SUCCESS)
