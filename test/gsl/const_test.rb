require 'test_helper'

class ConstTest < GSL::TestCase

  def test_speed_of_light
    c   = GSL::CONST::MKSA::SPEED_OF_LIGHT
    eps = GSL::CONST::MKSA::VACUUM_PERMITTIVITY
    mu  = GSL::CONST::MKSA::VACUUM_PERMEABILITY

    assert_rel c, 1.0 / GSL.sqrt(eps * mu), 1e-6, 'speed of light (mks)'
  end

  def test_light_year
    ly = GSL::CONST::CGSM::LIGHT_YEAR
    c  = GSL::CONST::CGSM::SPEED_OF_LIGHT
    y  = 365.2425 * GSL::CONST::CGSM::DAY

    assert_rel ly, c * y, 1e-6, 'light year (cgs)'
  end

  def test_kilo
    micro = GSL::CONST::NUM::MICRO
    mega  = GSL::CONST::NUM::MEGA
    kilo  = GSL::CONST::NUM::KILO

    assert_rel mega / kilo, 1 / (micro * kilo), 1e-10, 'kilo (mega/kilo, 1/(micro*kilo))'
  end

end
