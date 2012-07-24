#!/usr/bin/env ruby
# Ruby/GSL implementation of GSL "const/test.c"
require("gsl")
require("./gsl_test.rb")
include GSL::Test
include Math
include GSL::CONST

GSL::IEEE::env_setup()

c = MKSA::SPEED_OF_LIGHT
eps = MKSA::VACUUM_PERMITTIVITY
mu = MKSA::VACUUM_PERMEABILITY
GSL::Test.test_rel(c, 1.0/sqrt(eps*mu), 1e-6, "speed of light (mks)")

ly = CGSM::LIGHT_YEAR
c = CGSM::SPEED_OF_LIGHT
y = 365.2425 * CGSM::DAY
GSL::Test.test_rel(ly, c * y, 1e-6, "light year (cgs)")

micro = NUM::MICRO
mega = NUM::MEGA
kilo = NUM::KILO
GSL::Test.test_rel(mega/kilo, 1/(micro*kilo), 1e-10, "kilo (mega/kilo, 1/(micro*kilo))");
