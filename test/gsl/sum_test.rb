require 'test_helper'

class SumTest < GSL::TestCase

  N = 50

  def setup
    @t = GSL::Vector.alloc(N)
    N.times { |n| @t[n] = 1.0 / (n + 1.0) ** 2 }
  end

  def _fill_vector(x = 1.0)
    @t[0] = x
    1.upto(N - 1) { |n| @t[n] = yield n, @t[n - 1] }
  end

  def _test_sum(expected, desc)
    sum_accel, _ = GSL::Sum::Levin_utrunc.alloc(N).accel(@t)
    assert_rel sum_accel, expected, 1e-8, 'trunc result, %s' % desc

    sum_accel, err_est, = GSL::Sum::Levin_u.alloc(N).accel(@t)
    assert_rel sum_accel, expected, 1e-8, 'full result, %s' % desc

    sd_est = -Math.log10(err_est / sum_accel.abs)
    sd_actual = -Math.log10(GSL::DBL_EPSILON + ((sum_accel - expected) / expected).abs)

    refute sd_est > sd_actual + 1.0,
      'full significant digits, %s (%g vs %g)' % [desc, sd_est, sd_actual]
  end

  def test_zeta_2
    _test_sum(GSL::M_PI ** 2 / 6.0, 'zeta(2)')
  end

  def test_exp_10(x = 10.0)
    _fill_vector { |n, t| t * (x / n) }
    _test_sum(Math.exp(x), 'exp(%d)' % x)
  end

  def test_exp_neg_10
    test_exp_10(-10.0)
  end

  def test_log(x = 0.5)
    _fill_vector(x) { |n, t| t * (x * n) / (n + 1.0) }
    _test_sum(-Math.log(1 - x), "-log(#{1 - x})")
  end

  def test_log2
    test_log(-1.0)
  end

  def test_asymptotic_series
    m = GSL::M_PI ** 2
    _fill_vector(3.0 / m) { |n, t| -t * (4.0 * (n + 1.0) - 1.0) / m }
    _test_sum(0.192594048773, 'asymptotic series')
  end

  def test_eulers_constant
    _fill_vector { |n, _| 1 / (n + 1.0) + Math.log(n / (n + 1.0)) }
    _test_sum(0.5772156649015328606065120900824, "Euler's constant")
  end

  def test_eta
    N.times { |n| @t[n] = (n % 2 == 1 ? -1 : 1) * 1.0 / GSL.sqrt(n + 1.0) }
    _test_sum(0.6048986434216305, 'eta(1/2)')
  end

end
