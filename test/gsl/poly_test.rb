require 'test_helper'

class PolyTest < GSL::TestCase

  EPS = 100.0 * GSL::DBL_EPSILON

  def test_poly
    c = GSL::Poly.alloc(1.0, 0.5, 0.3)
    x = 0.5
    y = c.eval(x)
    assert_rel y, 1 + 0.5 * x + 0.3 * x * x, EPS, 'gsl_poly_eval({1, 0.5, 0.3}, 0.5)'

    d = GSL::Poly.alloc( 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1)
    x = 1.0
    y = GSL::Poly.eval(d, x)
    assert_rel y, 1.0, EPS, 'gsl_poly_eval({1,-1, 1, -1, 1, -1, 1, -1, 1, -1, 1}, 1.0)'

    x0, x1 = GSL::Poly.solve_quadratic(4.0, -20.0, 25.0).to_a
    assert_rel x0, 2.5, 1e-9, 'x0, (2x - 5)^2 = 0'
    assert_rel x1, 2.5, 1e-9, 'x1, (2x - 5)^2 = 0'

    x0, x1 = GSL::Poly.solve_quadratic(4.0, 7.0, 0.0).to_a
    assert_rel x0, -1.75, 1e-9, 'x0, x(4x + 7) = 0'
    assert_rel x1, 0.0, 1e-9, 'x1, x(4x + 7) = 0'

    x0, x1 = GSL::Poly.solve_quadratic(5.0, 0.0, -20.0).to_a
    assert_rel x0, -2.0, 1e-9, 'x0, 5 x^2 = 20'
    assert_rel x1, 2.0, 1e-9, 'x1, 5 x^2 = 20'

    # Quadratic single real root (technically not a quadratic)
    x0, x1 = GSL::Poly.solve_quadratic(0.0, 1.0, 0.0).to_a
    assert_rel x0, 0.0, 0, 'x0, x = 0'
    assert x1.nil?, 'x1, x = 0 is nil'

    # Quadratic no real root
    x0, x1 = GSL::Poly.solve_quadratic(1.0, 0.0, 1.0).to_a
    assert x0.nil?, 'x0, x^2 = -1 is nil'
    assert x1.nil?, 'x1, x^2 = -1 is nil'

    x0, x1, x2 = GSL::Poly.solve_cubic(0.0, 0.0, -27.0).to_a
    assert_rel x0, 3.0, 1e-9, 'x0, x^3 = 27'

    # Cubic triple real root
    x0, x1, x2 = GSL::Poly.solve_cubic(-51.0, 867.0, -4913.0).to_a
    assert_rel x0, 17.0, 1e-9, 'x0, (x-17)^3=0'
    assert_rel x1, 17.0, 1e-9, 'x1, (x-17)^3=0'
    assert_rel x2, 17.0, 1e-9, 'x2, (x-17)^3=0'

    # Cubic double real root plus single real root
    x0, x1, x2 = GSL::Poly.solve_cubic(-57.0, 1071.0, -6647.0).to_a
    assert_rel x0, 17.0, 1e-9, 'x0, (x-17)(x-17)(x-23)=0'
    assert_rel x1, 17.0, 1e-9, 'x1, (x-17)(x-17)(x-23)=0'
    assert_rel x2, 23.0, 1e-9, 'x2, (x-17)(x-17)(x-23)=0'

    x0, x1, x2 = GSL::Poly.solve_cubic(-11.0, -493.0, +6647.0).to_a
    assert_rel x0, -23.0, 1e-9, 'x0, (x+23)(x-17)(x-17)=0'
    assert_rel x1, 17.0, 1e-9, 'x1, (x+23)(x-17)(x-17)=0'
    assert_rel x2, 17.0, 1e-9, 'x2, (x+23)(x-17)(x-17)=0'

    x0, x1, x2 = GSL::Poly.solve_cubic(-143.0, 5087.0, -50065.0).to_a
    assert_rel x0, 17.0, 1e-9, 'x0, (x-17)(x-31)(x-95)=0'
    assert_rel x1, 31.0, 1e-9, 'x1, (x-17)(x-31)(x-95)=0'
    assert_rel x2, 95.0, 1e-9, 'x2, (x-17)(x-31)(x-95)=0'

    x0, x1, x2 = GSL::Poly.solve_cubic(-109.0, 803.0, 50065.0).to_a
    assert_rel x0, -17.0, 1e-9, 'x0, (x+17)(x-31)(x-95)=0'
    assert_rel x1, 31.0, 1e-9, 'x1, (x+17)(x-31)(x-95)=0'
    assert_rel x2, 95.0, 1e-9, 'x2, (x+17)(x-31)(x-95)=0'

    # Cubic double real root only is impossible

    # Cubic single real root (and two complex roots, not returned)
    x0, x1, x2 = GSL::Poly.solve_cubic(0.0, 0.0, -1.0).to_a
    assert_rel x0, 1.0, 1e-9, 'x0, x^3 = 1'
    assert x1.nil?, 'x1, x^3 = 1 is nil'
    assert x2.nil?, 'x2, x^3 = 1 is nil'

    # Cubic no real root is impossible

    z = GSL::Poly::Complex.solve_quadratic(4.0, -20.0, 26.0)
    assert_rel z[0].re, 2.5, 1e-9, 'z0.real, (2x - 5)^2 = -1'
    assert_rel z[0].im, -0.5, 1e-9, 'z0.imag, (2x - 5)^2 = -1'
    assert_rel z[1].re, 2.5, 1e-9, 'z1.real, (2x - 5)^2 = -1'
    assert_rel z[1].im, 0.5, 1e-9, 'z1.imag, (2x - 5)^2 = -1'

    z = GSL::Poly.complex_solve_quadratic(4.0, -20.0, 25.0)
    assert_rel z[0].re, 2.5, 1e-9, 'z0.real, (2x - 5)^2 = 0'
    assert_rel z[0].im, 0.0, 1e-9, 'z0.imag (2x - 5)^2 = 0'
    assert_rel z[1].re, 2.5, 1e-9, 'z1.real, (2x - 5)^2 = 0'
    assert_rel z[1].im, 0.0, 1e-9, 'z1.imag (2x - 5)^2 = 0'
    assert_equal z[0].re, z[1].re, 'z0.real == z1.real, (2x - 5)^2 = 0'
    assert_equal z[1].im, z[1].im, 'z0.imag == z1.imag, (2x - 5)^2 = 0'

    z = GSL::Poly.complex_solve_quadratic(4.0, -20.0, 21.0)
    assert_rel z[0].re, 1.5, 1e-9, 'z0.real, (2x - 5)^2 = 4'
    assert_rel z[0].im, 0.0, 1e-9, 'z0.imag, (2x - 5)^2 = 4'
    assert_rel z[1].re, 3.5, 1e-9, 'z1.real, (2x - 5)^2 = 4'
    assert_rel z[1].im, 0.0, 1e-9, 'z1.imag, (2x - 5)^2 = 4'

    z = GSL::Poly.complex_solve_quadratic(4.0, 7.0, 0.0)
    assert_rel z[0].re, -1.75, 1e-9, 'z[0].real, x(4x + 7) = 0'
    assert_rel z[0].im, 0.0, 1e-9, 'z[0].imag, x(4x + 7) = 0'
    assert_rel z[1].re, 0.0, 1e-9, 'z[1].real, x(4x + 7) = 0'
    assert_rel z[1].im, 0.0, 1e-9, 'z[1].imag, x(4x + 7) = 0'

    z = GSL::Poly.complex_solve_quadratic(5.0, 0.0, -20.0)
    assert_rel z[0].re, -2.0, 1e-9, 'z[0].real, 5 x^2 = 20'
    assert_rel z[0].im, 0.0, 1e-9, 'z[0].imag, 5 x^2 = 20'
    assert_rel z[1].re, 2.0, 1e-9, 'z[1].real, 5 x^2 = 20'
    assert_rel z[1].im, 0.0, 1e-9, 'z[1].imag, 5 x^2 = 20'

    z = GSL::Poly.complex_solve_quadratic(5.0, 0.0, 20.0)
    assert_rel z[0].re, 0.0, 1e-9, 'z[0].real, 5 x^2 = -20'
    assert_rel z[0].im, -2.0, 1e-9, 'z[0].imag, 5 x^2 = -20'
    assert_rel z[1].re, 0.0, 1e-9, 'z[1].real, 5 x^2 = -20'
    assert_rel z[1].im, 2.0, 1e-9, 'z[1].imag, 5 x^2 = -20'

    # Quadratic single complex root (technically not quadratic and root not
    # complex since imaginary component is 0, but the data type is complex)
    z = GSL::Poly.complex_solve_quadratic(0.0, 1.0, -1.0)
    assert_rel z[0].re, 1.0, 1e-9, 'z[0].real, x = 1 (complex)'
    assert_rel z[0].im, 0.0, 0.0,  'z[0].imag, x = 1 (complex)'
    assert x1.nil?, 'z[1], x = 0 is nil'

    # Quadratic no complex root (technically not quadratic)
    z = GSL::Poly.complex_solve_quadratic(0.0, 0.0, 1.0)
    assert z[0].nil?, 'z[0], 1 = 0 is nil'
    assert z[1].nil?, 'z[1], 1 = 0 is nil'

    z = GSL::Poly.complex_solve_cubic(0.0, 0.0, -27.0)
    assert_rel z[0].re, -1.5, 1e-9, 'z[0].real, x^3 = 27'
    assert_rel z[0].im, -1.5 * Math.sqrt(3.0), 1e-9, 'z[0].imag, x^3 = 27'
    assert_rel z[1].re, -1.5, 1e-9, 'z[1].real, x^3 = 27'
    assert_rel z[1].im, 1.5 * Math.sqrt(3.0), 1e-9, 'z[1].imag, x^3 = 27'
    assert_rel z[2].re, 3.0, 1e-9, 'z[2].real, x^3 = 27'
    assert_rel z[2].im, 0.0, 1e-9, 'z[2].imag, x^3 = 27'

    z = GSL::Poly.complex_solve_cubic(-1.0, 1.0, 39.0)
    assert_rel z[0].re, -3.0, 1e-9, 'z[0].real, (x+3)(x^2+1) = 0'
    assert_rel z[0].im, 0.0, 1e-9, 'z[0].imag, (x+3)(x^2+1) = 0'
    assert_rel z[1].re, 2.0, 1e-9, 'z[1].real, (x+3)(x^2+1) = 0'
    assert_rel z[1].im, -3.0, 1e-9, 'z[1].imag, (x+3)(x^2+1) = 0'
    assert_rel z[2].re, 2.0, 1e-9, 'z[2].real, (x+3)(x^2+1) = 0'
    assert_rel z[2].im, 3.0, 1e-9, 'z[2].imag, (x+3)(x^2+1) = 0'

    z = GSL::Poly.complex_solve_cubic(-51.0, 867.0, -4913.0)
    assert_rel z[0].re, 17.0, 1e-9, 'z[0].real, (x-17)^3=0'
    assert_rel z[0].im, 0.0, 1e-9, 'z[0].imag, (x-17)^3=0'
    assert_rel z[1].re, 17.0, 1e-9, 'z[1].real, (x-17)^3=0'
    assert_rel z[1].im, 0.0, 1e-9, 'z[1].imag, (x-17)^3=0'
    assert_rel z[2].re, 17.0, 1e-9, 'z[2].real, (x-17)^3=0'
    assert_rel z[2].im, 0.0, 1e-9, 'z[2].imag, (x-17)^3=0'

    z = GSL::Poly.complex_solve_cubic(-57.0, 1071.0, -6647.0)
    assert_rel z[0].re, 17.0, 1e-9, 'z[0].real, (x-17)(x-17)(x-23)=0'
    assert_rel z[0].im, 0.0, 1e-9, 'z[0].imag, (x-17)(x-17)(x-23)=0'
    assert_rel z[1].re, 17.0, 1e-9, 'z[1].real, (x-17)(x-17)(x-23)=0'
    assert_rel z[1].im, 0.0, 1e-9, 'z[1].imag, (x-17)(x-17)(x-23)=0'
    assert_rel z[2].re, 23.0, 1e-9, 'z[2].real, (x-17)(x-17)(x-23)=0'
    assert_rel z[2].im, 0.0, 1e-9, 'z[2].imag, (x-17)(x-17)(x-23)=0'

    z = GSL::Poly.complex_solve_cubic(-11.0, -493.0, +6647.0)
    assert_rel z[0].re, -23.0, 1e-9, 'z[0].real, (x+23)(x-17)(x-17)=0'
    assert_rel z[0].im, 0.0, 1e-9, 'z[0].imag, (x+23)(x-17)(x-17)=0'
    assert_rel z[1].re, 17.0, 1e-9, 'z[1].real, (x+23)(x-17)(x-17)=0'
    assert_rel z[1].im, 0.0, 1e-9, 'z[1].imag, (x+23)(x-17)(x-17)=0'
    assert_rel z[2].re, 17.0, 1e-9, 'z[2].real, (x+23)(x-17)(x-17)=0'
    assert_rel z[2].im, 0.0, 1e-9, 'z[2].imag, (x+23)(x-17)(x-17)=0'

    z = GSL::Poly.complex_solve_cubic(-143.0, 5087.0, -50065.0)
    assert_rel z[0].re, 17.0, 1e-9, 'z[0].real, (x-17)(x-31)(x-95)=0'
    assert_rel z[0].im, 0.0, 1e-9, 'z[0].imag, (x-17)(x-31)(x-95)=0'
    assert_rel z[1].re, 31.0, 1e-9, 'z[1].real, (x-17)(x-31)(x-95)=0'
    assert_rel z[1].im, 0.0, 1e-9, 'z[1].imag, (x-17)(x-31)(x-95)=0'
    assert_rel z[2].re, 95.0, 1e-9, 'z[2].real, (x-17)(x-31)(x-95)=0'
    assert_rel z[2].im, 0.0, 1e-9, 'z[2].imag, (x-17)(x-31)(x-95)=0'

    a = GSL::Poly.alloc(-120, 274, -225, 85, -15, 1)
    w = GSL::Poly::Complex::Workspace.alloc(a.size)
    z = GSL::Poly.complex_solve(a, 6, w)
    assert_rel z[0].re, 1.0, 1e-9, 'z[0].real, 5th-order polynomial'
    assert_rel z[0].im, 0.0, 1e-9, 'z[0].imag, 5th-order polynomial'
    assert_rel z[1].re, 2.0, 1e-9, 'z[1].real, 5th-order polynomial'
    assert_rel z[1].im, 0.0, 1e-9, 'z[1].imag, 5th-order polynomial'
    assert_rel z[2].re, 3.0, 1e-9, 'z[2].real, 5th-order polynomial'
    assert_rel z[2].im, 0.0, 1e-9, 'z[2].imag, 5th-order polynomial'
    assert_rel z[3].re, 4.0, 1e-9, 'z3.real, 5th-order polynomial'
    assert_rel z[3].im, 0.0, 1e-9, 'z3.imag, 5th-order polynomial'
    assert_rel z[4].re, 5.0, 1e-9, 'z4.real, 5th-order polynomial'
    assert_rel z[4].im, 0.0, 1e-9, 'z4.imag, 5th-order polynomial'

    a = GSL::Poly.alloc(1, 0, 0, 0, 1, 0, 0, 0, 1)
    w = GSL::Poly::Complex::Workspace.alloc(a.size)
    c = 0.5
    s = Math.sqrt(3) / 2
    z = GSL::Poly.complex_solve(a, w)
    assert_rel z[0].re, -s, 1e-9, 'z[0].real, 8th-order polynomial'
    assert_rel z[0].im,  c, 1e-9, 'z[0].imag, 8th-order polynomial'
    assert_rel z[1].re, -s, 1e-9, 'z[1].real, 8th-order polynomial'
    assert_rel z[1].im, -c, 1e-9, 'z[1].imag, 8th-order polynomial'
    assert_rel z[2].re, -c, 1e-9, 'z[2].real, 8th-order polynomial'
    assert_rel z[2].im,  s, 1e-9, 'z[2].imag, 8th-order polynomial'
    assert_rel z[3].re, -c, 1e-9, 'z3.real, 8th-order polynomial'
    assert_rel z[3].im, -s, 1e-9, 'z3.imag, 8th-order polynomial'
    assert_rel z[4].re,  c, 1e-9, 'z4.real, 8th-order polynomial'
    assert_rel z[4].im,  s, 1e-9, 'z4.imag, 8th-order polynomial'
    assert_rel z[5].re,  c, 1e-9, 'z5.real, 8th-order polynomial'
    assert_rel z[5].im, -s, 1e-9, 'z5.imag, 8th-order polynomial'
    assert_rel z[6].re,  s, 1e-9, 'z6.real, 8th-order polynomial'
    assert_rel z[6].im,  c, 1e-9, 'z6.imag, 8th-order polynomial'
    assert_rel z[7].re,  s, 1e-9, 'z7.real, 8th-order polynomial'
    assert_rel z[7].im, -c, 1e-9, 'z7.imag, 8th-order polynomial'

    xa = GSL::Poly.alloc(0.16, 0.97, 1.94, 2.74, 3.58, 3.73, 4.70)
    ya = GSL::Poly.alloc(0.73, 1.11, 1.49, 1.84, 2.30, 2.41, 3.07)
    dd_expected = GSL::Vector.alloc(7.30000000000000e-01,
                                    4.69135802469136e-01,
                                   -4.34737219941284e-02,
                                    2.68681098870099e-02,
                                   -3.22937056934996e-03,
                                    6.12763259971375e-03,
                                   -6.45402453527083e-03)
    dd = GSL::Poly.dd_init(xa, ya)

    7.times { |i|
      assert_rel dd[i], dd_expected[i], 1e-10, "divided difference dd[#{i}]"
    }

    7.times { |i|
      y = dd.eval(xa, xa[i])
      assert_rel y, ya[i], 1e-10, "divided difference y[#{i}]"
    }

    coeff = dd.taylor(1.5, xa)

    7.times { |i|
      y = coeff.eval(xa[i] - 1.5)
      assert_rel y, ya[i], 1e-10, "taylor expansion about 1.5 y[#{i}]"
    }

    # Added GSL-1.12.90 (gsl-1.13)
    # gsl_poly_eval_derivs()
    return unless GSL::Poly.method_defined?(:eval_derivs)

    c = GSL::Vector.alloc([1.0, -2.0, 3.0, -4.0, 5.0, -6.0])
    x = -0.5

    dc = GSL::Poly.eval_derivs(c, x)

    assert_rel dc[0], c[0] + c[1] * x + c[2] * x * x + c[3] * x * x * x + c[4] * x * x * x * x + c[5] * x * x * x * x * x , EPS, 'eval_derivs({+1, -2, +3, -4, +5, -6}, 3.75)'
    assert_rel dc[1], c[1] + 2.0 * c[2] * x + 3.0 * c[3] * x * x + 4.0 * c[4] * x * x * x + 5.0 * c[5] * x * x * x * x , EPS, 'eval_derivs({+1, -2, +3, -4, +5, -6} deriv 1, -12.375)'
    assert_rel dc[2], 2.0 * c[2] + 3.0 * 2.0 * c[3] * x + 4.0 * 3.0 * c[4] * x * x + 5.0 * 4.0 * c[5] * x * x * x , EPS, 'eval_derivs({+1, -2, +3, -4, +5, -6} deriv 2, +48.0)'
    assert_rel dc[3], 3.0 * 2.0 * c[3] + 4.0 * 3.0 * 2.0 * c[4] * x + 5.0 * 4.0 * 3.0 * c[5] * x * x , EPS, 'eval_derivs({+1, -2, +3, -4, +5, -6} deriv 3, -174.0)'
    assert_rel dc[4], 4.0 * 3.0 * 2.0 * c[4] + 5.0 * 4.0 * 3.0 * 2.0 * c[5] * x, EPS, 'eval_derivs({+1, -2, +3, -4, +5, -6} deriv 4, +480.0)'
    assert_rel dc[5], 5.0 * 4.0 * 3.0 * 2.0 * c[5] , EPS, 'eval_derivs({+1, -2, +3, -4, +5, -6} deriv 5, -720.0)'

    # Test Poly::fit and Poly::wfit
    x = GSL::Vector[0, 2, 2]
    y = GSL::Vector[0, 1, -1]
    coef, cov, chisq, status = GSL::Poly.fit(x, y, 1)
    assert_rel coef[0], 0.0, 1e-9, 'y intercept == 0'
    assert_rel coef[1], 0.0, 1e-9, 'slope == 0'
    assert_rel chisq, 2.0, 1e-9, 'chisq == 2'

    w = GSL::Vector[1, 1, 0]
    coef, cov, chisq, status = GSL::Poly.wfit(x, w, y, 1)
    assert_rel coef[0], 0.0, 1e-9, 'y intercept == 0'
    assert_rel coef[1], 0.5, 1e-9, 'slope == 0.5'
    assert_rel chisq, 0.0, 1e-9, 'chisq == 0'
  end

  def test_special
    hermit = [
      GSL::Poly::Int[1],
      GSL::Poly::Int[0, 2],
      GSL::Poly::Int[-2, 0, 4],
      GSL::Poly::Int[0, -12, 0, 8],
      GSL::Poly::Int[12, 0, -48, 0, 16],
      GSL::Poly::Int[0, 120, 0, -160, 0, 32],
      GSL::Poly::Int[-120, 0, 720, 0, -480, 0, 64]
    ]
    hermit[0][0] = 1

    6.times { |i|
      hn = GSL::Poly.hermite(i)
      assert_equal hermit[i], hn, "Hermite polynomial, n = #{i}"
    }

    laguerre = [
      GSL::Poly::Int[1],
      GSL::Poly::Int[1, -1],
      GSL::Poly::Int[2, -4, 1],
      GSL::Poly::Int[6, -18, 9, -1],
      GSL::Poly::Int[24, -96, 72, -16, 1],
      GSL::Poly::Int[120, -600, 600, -200, 25, -1],
      GSL::Poly::Int[720, -4320, 5400, -2400, 450, -36, 1]
    ]
    laguerre[0][0] = 1

    7.times { |i|
      hn = GSL::Poly.laguerre(i)
      assert_equal laguerre[i], hn, "Laguerre polynomial, n = #{i}"
    }

    cheb = [
      GSL::Poly::Int[1],
      GSL::Poly::Int[0, 1],
      GSL::Poly::Int[-1, 0, 2],
      GSL::Poly::Int[0, -3, 0, 4],
      GSL::Poly::Int[1, 0, -8, 0, 8],
      GSL::Poly::Int[0, 5, 0, -20, 0, 16],
      GSL::Poly::Int[-1, 0, 18, 0, -48, 0, 32]
    ]
    cheb[0][0] = 1

    7.times { |i|
      hn = GSL::Poly.cheb(i)
      assert_equal cheb[i], hn, "Chebyshev polynomial, n = #{i}"
    }

    cheb_ii = [
      GSL::Poly::Int[1],
      GSL::Poly::Int[0, 2],
      GSL::Poly::Int[-1, 0, 4],
      GSL::Poly::Int[0, -4, 0, 8],
      GSL::Poly::Int[1, 0, -12, 0, 16],
      GSL::Poly::Int[0, 6, 0, -32, 0, 32],
      GSL::Poly::Int[-1, 0, 24, 0, -80, 0, 64]
    ]
    cheb_ii[0][0] = 1

    7.times { |i|
      hn = GSL::Poly.cheb_II(i)
      assert_equal cheb_ii[i], hn, "Chebyshev II polynomial, n = #{i}"
    }
  end

end
