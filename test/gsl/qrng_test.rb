require 'test_helper'

class QRngTest < GSL::TestCase

  def _test_sobol(v, d, g = nil)
    status = 0

    if g
      g.init
      reinitialized = true
    else
      g = GSL::QRng.alloc(GSL::QRng::SOBOL, d)
    end

    3.times { g.get(v) }

    status += (v[0] != 0.25 || v[1] != 0.75) ? 1 : 0
    g.get(v)

    status += (v[0] != 0.375 || v[1] != 0.375) ? 1 : 0
    assert status.zero?, "Sobol d=#{d}#{' (reinitialized)' if reinitialized}"

    g
  end

  def test_sobol
    v = GSL::Vector.alloc(3)

    _test_sobol(v, 2)
    _test_sobol(v, 3, _test_sobol(v, 3))
  end

  def _test_nied2(v, d, g = nil)
    status = 0

    if g
      g.init
      reinitialized = true
    else
      g = GSL::QRng.alloc(GSL::QRng::NIEDERREITER_2, d)
    end

    3.times { g.get(v) }

    status += (v[0] != 0.75 || v[1] != 0.25) ? 1 : 0
    g.get(v)

    status += (v[0] != 0.25 || v[1] != 0.75) ? 1 : 0
    3.times { g.get(v) }

    status += (v[0] != 0.625 || v[1] != 0.125) ? 1 : 0
    assert status.zero?, "Niederreiter d=#{d}#{' (reinitialized)' if reinitialized}"

    g
  end

  def test_nied2
    v = GSL::Vector.alloc(3)

    _test_nied2(v, 2)
    _test_nied2(v, 3, _test_nied2(v, 3))
  end

  def _test_hdsobol(v, d, g = nil)
    status = 0

    if g
      g.init
      reinitialized = true
    else
      g = GSL::QRng.alloc(GSL::QRng::HDSOBOL, d)
    end

    3.times { g.get(v) }

    status += (v[0] != 0.25 || v[1] != 0.75) ? 1 : 0
    g.get(v)

    status += (v[0] != 0.375 || v[1] != 0.375) ? 1 : 0
    assert status.zero?, "HDSobol d=#{d}#{' (reinitialized)' if reinitialized}"

    g
  end

  def test_hdsobol
    return unless GSL::QRng.const_defined?(:HDSOBOL)

    v = GSL::Vector.alloc(3)

    _test_hdsobol(v, 2)
    _test_hdsobol(v, 3, _test_hdsobol(v, 3))
  end

end
