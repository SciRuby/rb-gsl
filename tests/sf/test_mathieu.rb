#!/usr/bin/env ruby
require("gsl")
require("../gsl_test2.rb")
require("./gsl_test_sf.rb")
include GSL::Test
include Math
include GSL::Test::Sf
s = 0

TEST_SNGL = 1.0e-06

TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(0, 0.0, 0.0)", 0.7071067811865475, TEST_SNGL, GSL::SUCCESS)
  TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(0, 0.0, GSL::M_PI_2)",
          0.7071067811865475, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(0, 5.0, 0.0)",
          0.04480018165188902, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(0, 5.0, GSL::M_PI_2)",
          1.334848674698019, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(0, 10.0, 0.0)",
          0.007626517570935782, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(0, 10.0, GSL::M_PI_2)",
          1.468660470712856, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(0, 15.0, 0.0)",
          0.001932508315204592, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(0, 15.0, GSL::M_PI_2)",
          1.550108146686649, TEST_SNGL, GSL::SUCCESS);
 TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(0, 20.0, 0.0)",
         0.0006037438292242197, TEST_SNGL, GSL::SUCCESS);
 TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(0, 20.0, GSL::M_PI_2)",
         1.609890857395926, TEST_SNGL, GSL::SUCCESS);
 TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(0, 25.0, 0.0)",
         0.0002158630184146612, TEST_SNGL, GSL::SUCCESS);
 TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(0, 25.0, GSL::M_PI_2)",
         1.657510298323475, TEST_SNGL, GSL::SUCCESS);
 TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(1, 0.0, 0.0)",
         1.00000000, TEST_SNGL, GSL::SUCCESS);
 TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(1, 5.0, 0.0)",
         0.2565428793223637, TEST_SNGL, GSL::SUCCESS);
 TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(1, 10.0, 0.0)",
         0.05359874774717657, TEST_SNGL, GSL::SUCCESS);
 TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(1, 15.0, 0.0)",
         0.01504006645382623, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_se_e", "(5, 0.0, GSL::M_PI_2)",
          1.0000000, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_se_e", "(5, 5.0, GSL::M_PI_2)",
          0.9060779302023551, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_se_e", "(5, 10.0, GSL::M_PI_2)",
          0.8460384335355106, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_se_e", "(5, 15.0, GSL::M_PI_2)",
          0.837949340012484, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_se_e", "(5, 20.0, GSL::M_PI_2)",
          0.8635431218533667, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_se_e", "(5, 25.0, GSL::M_PI_2)",
          0.8992683245108413, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(10, 0.0, 0.0)",
          1.00000000, TEST_SNGL, GSL::SUCCESS);
  TEST_SF(s, "GSL::Sf::mathieu_ce_e", "(10, 0.0, GSL::M_PI_2)",
          -1.00000000, TEST_SNGL, GSL::SUCCESS);                   
 