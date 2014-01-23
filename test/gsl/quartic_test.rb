require 'test_helper'

class QuarticTest < GSL::TestCase

  def _test_quartic(n, *sol)
    z = GSL::Poly.complex_solve_quartic(0.0, 0.0, 0.0, n)
    desc = "Four roots, x^4 #{n}"

    assert_rel z[0].re, sol[0], 1e-9, "#{desc}: z0.real"
    assert_rel z[0].im, sol[1], 1e-9, "#{desc}: z0.imag"
    assert_rel z[1].re, sol[2], 1e-9, "#{desc}: z1.real"
    assert_rel z[1].im, sol[3], 1e-9, "#{desc}: z1.imag"
    assert_rel z[2].re, sol[4], 1e-9, "#{desc}: z2.real"
    assert_rel z[2].im, sol[5], 1e-9, "#{desc}: z2.imag"
    assert_rel z[3].re, sol[6], 1e-9, "#{desc}: z3.real"
    assert_rel z[3].im, sol[7], 1e-9, "#{desc}: z3.imag"
  end

  def test_quartic
    return unless GSL::Poly.method_defined?(:complex_solve_quartic)

    sol = 3.0 / Math.sqrt(2.0)

    _test_quartic(-81.0, -3.0, 0.0, 0.0, -3.0, 0.0, 3.0, 3.0, 0.0)
    _test_quartic(81.0, -sol, -sol, -sol, sol, sol, -sol, sol, sol)
  end

end
