require 'test_helper'

class ChebTest < GSL::TestCase

  def test_cheb
    tol, pi, order = 100.0 * GSL::DBL_EPSILON, GSL::M_PI, 40

    cs = GSL::Cheb.alloc(order)

    cs.init(GSL::Function.alloc { |x| 1.0 }, -1.0, 1.0)
    order.times { |i|
      assert_abs cs.c[i], i == 0 ? 2.0 : 0.0, tol, 'c[%d] for T_0(x)' % i
    }

    cs.init(GSL::Function.alloc { |x| x }, -1.0, 1.0)
    order.times { |i|
      assert_abs cs.c[i], i == 1 ? 1.0 : 0.0, tol, 'c[%d] for T_1(x)' % i
    }

    cs.init(GSL::Function.alloc { |x| 2.0 * x * x - 1.0 }, -1.0, 1.0)
    order.times { |i|
      assert_abs cs.c[i], i == 2 ? 1.0 : 0.0, tol, 'c[%d] for T_2(x)' % i
    }

    cs.init(GSL::Function.alloc { |x| Math.sin(x) }, -pi, pi)
    assert_abs cs.c[0], 0.0, tol, 'c[0] for F_sin(x)'
    assert_abs cs.c[1], 5.69230686359506e-01, tol, 'c[1] for F_sin(x)'
    assert_abs cs.c[2], 0.0, tol, 'c[2] for F_sin(x)'
    assert_abs cs.c[3], -6.66916672405979e-01, tol, 'c[3] for F_sin(x)'
    assert_abs cs.c[4], 0.0, tol, 'c[4] for F_sin(x)'
    assert_abs cs.c[5], 1.04282368734237e-01, tol, 'c[5] for F_sin(x)'

    x = -pi
    while x < pi
      assert_abs cs.eval(x), Math.sin(x), tol, 'GSL::Cheb#eval, sin(%.3g)' % x
      x += pi / 100.0
    end

    x = -pi
    while x < pi
      r, e = cs.eval_err(x)

      assert_abs r, Math.sin(x), tol, 'GSL::Cheb#eval_err, sin(%.3g)' % x
      assert_factor((r - Math.sin(x)).abs + GSL::DBL_EPSILON, e, 10.0,
        'GSL::Cheb#eval_err, error sin(%.3g)' % x)

      x += pi / 100.0
    end

    x = -pi
    while x < pi
      assert_abs cs.eval_n(25, x), Math.sin(x), tol, 'GSL::Cheb#eval_n, sin(%.3g)' % x
      x += pi / 100.0
    end

    x = -pi
    while x < pi
      r, e = cs.eval_n_err(25, x)

      assert_abs r, Math.sin(x), tol, 'GSL::Cheb#eval_n_err, sin(%.3g)' % x
      assert_factor((r - Math.sin(x)).abs + GSL::DBL_EPSILON, e, 10.0,
        'GSL::Cheb#eval_n_err, error sin(%.3g)' % x)

      x += pi / 100.0
    end

    csd, x = cs.calc_deriv, -pi
    while x < pi
      assert_abs csd.eval(x), Math.cos(x), 1600 * tol, 'GSL::Cheb#eval, deriv sin(%.3g)' % x
      x += pi / 100.0
    end

    csi, x = cs.calc_integ, -pi
    while x < pi
      assert_abs csi.eval(x), -(1 + Math.cos(x)), tol, 'GSL::Cheb#eval, integ sin(%.3g)' % x
      x += pi / 100.0
    end
  end

end
