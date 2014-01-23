require 'test_helper'

class MinTest < GSL::TestCase

  EPSABS = 0.001
  EPSREL = 0.001

  MAX_ITERATIONS = 100

  def _test_f(type, desc, f, x_lower, mid, x_upper, min)
    s = GSL::Min::FMinimizer.alloc(type)
    s.set(@f[f], mid, x_lower, x_upper)

    status = iterations = 0

    begin
      iterations += 1
      status = s.iterate

      m = s.x_minimum
      a = s.x_lower
      b = s.x_upper

      refute a > b, 'interval is invalid (%g,%g)' % [a, b]
      refute m < a || m > b, 'm lies outside interval %g (%g,%g)' % [m, a, b]

      break if status == 1

      status = GSL::Min.test_interval(a, b, EPSABS, EPSREL)
    end while status == GSL::CONTINUE && iterations < MAX_ITERATIONS

    assert status.zero?, '%s, %s (%g obs vs %g expected)' % [s.name, desc, s.x_minimum, min]

    assert_tol m, min, 'incorrect precision (%g obs vs %g expected)' % [m, min]
  end

  def _test_f_e(type, desc, f, x_lower, mid, x_upper, min)
    s = GSL::Min::FMinimizer.alloc(type)

    assert_raises(GSL::ERROR::EINVAL, '%s, %s' % [s.name, desc]) {
      s.set(@f[f], mid, x_lower, x_upper)
    }

    status = iterations = 0

    begin
      iterations += 1
      s.iterate

      _ = s.x_minimum
      a = s.x_lower
      b = s.x_upper

      status = GSL::Min.test_interval(a, b, EPSABS, EPSREL)
    rescue
    end while status == GSL::CONTINUE && iterations < MAX_ITERATIONS

    assert status.zero?, '%s, %s' % [s.name, desc]
  end

  def setup
    @f = [GSL::Function.alloc { |x| Math.cos(x) }]
    @f << GSL::Function.alloc { |x| GSL.pow(x, 4.0) - 1 }
    @f << GSL::Function.alloc { |x| Math.sqrt(x.abs) }
    @f << GSL::Function.alloc { |x| x < 1.0 ? 1 : -Math.exp(-x) }
    @f << GSL::Function.alloc { |x| x - 30.0 / (1.0 + 1e5 * GSL.pow(x - 0.8, 2.0)) }
  end

  %w[goldensection brent quad_golden].each { |type|
    {
      'cos(x) [0 (3) 6]'        => [0,  0.0,  3.0,   6.0, GSL::M_PI],
      'x^4 - 1 [-3 (-1) 17]'    => [1, -3.0, -1.0,  17.0, 0.0],
      'sqrt(|x|) [-2 (-1) 1.5]' => [2, -2.0, -1.0,   1.5, 0.0],
      'func3(x) [-2 (3) 4]'     => [3, -2.0,  3.0,   4.0, 1.0],
      'func4(x) [0 (0.782) 1]'  => [4,  0,    0.782, 1.0, 0.8]
    }.each_with_index { |(desc, args), i|
      define_method("test_f_#{type}_#{i}") { _test_f(type, desc, *args) }
    }

    {
      'invalid range check [4, 0]'  => [0,  4.0, 3.0, 0.0, GSL::M_PI],
      'invalid range check [1, 1]'  => [0,  1.0, 1.0, 1.0, GSL::M_PI],
      'invalid range check [-1, 1]' => [0, -1.0, 0.0, 1.0, GSL::M_PI]
    }.each_with_index { |(desc, args), i|
      define_method("test_f_e_#{type}_#{i}") { _test_f_e(type, desc, *args) }
    }
  }

end
