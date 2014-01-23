require 'test_helper'

class OdeivTest < GSL::TestCase

  def setup
    @rhs_func_lin = GSL::Odeiv::System.alloc(lambda { |t, y, f|
      f[0] = 0.0
      f[1] = y[0]
      GSL::SUCCESS
    }, lambda { |t, y, dfdy, dfdt|
      dfdy.set(0, 0, 0.0)
      dfdy.set(0, 1, 0.0)
      dfdy.set(1, 0, 1.0)
      dfdy.set(1, 1, 0.0)
      dfdt[0] = 0.0
      dfdt[1] = 0.0
      GSL::SUCCESS
    }, 2)

    @rhs_func_sin = GSL::Odeiv::System.alloc(lambda { |t, y, f|
      f[0] = -y[1]
      f[1] = y[0]
      GSL::SUCCESS
    }, lambda { |t, y, dfdy, dfdt|
      dfdy.set(0, 0, 0.0)
      dfdy.set(0, 1, -1.0)
      dfdy.set(1, 0, 1.0)
      dfdy.set(1, 1, 0.0)
      dfdt[0] = 0.0
      dfdt[1] = 0.0
      GSL::SUCCESS
    }, 2)

    @rhs_func_exp = GSL::Odeiv::System.alloc(lambda { |t, y, f|
      f[0] = y[1]
      f[1] = y[0]
      GSL::SUCCESS
    }, lambda { |t, y, dfdy, dfdt|
      dfdy.set(0, 0, 0.0)
      dfdy.set(0, 1, 1.0)
      dfdy.set(1, 0, 1.0)
      dfdy.set(1, 1, 0.0)
      dfdt[0] = 0.0
      dfdt[1] = 0.0
      GSL::SUCCESS
    }, 2)

    @rhs_func_stiff = GSL::Odeiv::System.alloc(lambda { |t, y, f|
      f[0] = 998.0 * y[0] + 1998.0 * y[1]
      f[1] = -999.0 * y[0] - 1999.0 * y[1]
      GSL::SUCCESS
    }, lambda { |t, y, dfdy, dfdt|
      dfdy.set(0, 0, 998.0)
      dfdy.set(0, 1, 1998.0)
      dfdy.set(1, 0, -999.0)
      dfdy.set(1, 1, -1999.0)
      dfdt[0] = 0.0
      dfdt[1] = 0.0
      GSL::SUCCESS
    }, 2)
  end

  def self.define_test(test, type, h, b)
    define_method("test_#{test}_#{type}") { send("_test_#{test}", type, h, b) }
  end

  h, b = 1e-3, GSL::SQRT_DBL_EPSILON
  h2, b2 = h / 10, 1e-8

  %w[rk2 rk2imp rk4 rk4imp rkf45 rk8pd rkck bsimp gear1 gear2].each { |type|
    define_test(:stepper_err, type, h, b) unless type == 'bsimp'

    unless type == 'gear2'
      define_test(:stepper_linear, type, h,  b)
      define_test(:stepper_sin,    type, h2, b2)
      define_test(:stepper_exp,    type, h2, b2)
      define_test(:stepper_stiff,  type, h,  b)
    end

    define_test(:evolve_sin,    type, h, b)
    define_test(:evolve_exp,    type, h, b)
    define_test(:evolve_stiff1, type, h, b)
    define_test(:evolve_stiff5, type, h, b)
  }

  def _test_stepper_err(type, h, b)
    t = 0.0

    y = GSL::Vector.alloc(1.0, 0.0)
    yerr = GSL::Vector.alloc(2)

    stepper = GSL::Odeiv::Step.alloc(type, 2)

    while t < GSL::M_PI
      y1_t = y[1]
      dy_exp = Math.cos(t) * Math.sin(h) - 2 * Math.sin(t) * GSL.pow(Math.sin(h / 2), 2.0)

      stepper.apply(t, h, y, yerr, @rhs_func_sin)

      dy_t = y[1] - y1_t
      del = (dy_t - dy_exp).abs

      refute t > 0.1 && t < 0.2 && del > 10.0 * yerr[1].abs + GSL::DBL_EPSILON * y[1].abs,
        "#{stepper.name}, sine [0,pi], accurary of estimate error = #{del} vs #{yerr[1].abs}"

      t += h
    end
  end

  def _test_stepper_linear(type, h, b)
    delmax, n, t = 0.0, 0, 0.0

    stepper = GSL::Odeiv::Step.alloc(type, 2)

    y = GSL::Vector.alloc(1.0, 0.0)
    yerr = GSL::Vector.alloc(2)

    while t < 4.0
      stepper.apply(t, h, y, yerr, @rhs_func_lin)

      del = ((y[1] - (t + h)) / y[1]).abs
      delmax = GSL.MAX_DBL(del, delmax)

      refute del > (n + 1.0) * b,
        "#{stepper.name}, linear [0,4], max relative error = #{delmax}"

      n += 1
      t += h
    end
  end

  def _test_stepper_sin(type, h, b)
    delmax, n, t, pi = 0.0, 0, 0.0, GSL::M_PI

    stepper = GSL::Odeiv::Step.alloc(type, 2)

    y = GSL::Vector.alloc(1.0, 0.0)
    yerr = GSL::Vector.alloc(2)

    while t < pi
      sin_th = Math.sin(t + h)
      stepper.apply(t, h, y, yerr, @rhs_func_sin)

      del = ((y[1] - sin_th) / sin_th).abs
      delmax = GSL.MAX_DBL(del, delmax)

      refute del > (
        t < 0.5 * pi ? 1.0 + n :
        t < 0.7 * pi ? 1.0e+04 :
        t < 0.9 * pi ? 1.0e+06 : 1.0e+09) * b,
        "#{stepper.name}, sine [0,pi], max relative error = #{delmax}"

      n += 1
      t += h
    end

    refute delmax > 1.0e+09 * b,
      "#{stepper.name}, sine [0,pi], max relative error = #{delmax}"

    delmax = 0.0

    while t < 3 * pi
      stepper.apply(t, h, y, yerr, @rhs_func_sin)

      del = (y[1] - Math.sin(t)).abs
      delmax = GSL.MAX_DBL(del, delmax)

      n += 1
      t += h
    end

    refute del > n * 2.0 * b,
      "#{stepper.name}, sin [pi,3*pi], max absolute error = #{delmax}"
  end

  def _test_stepper_exp(type, h, b)
    delmax, n, t = 0.0, 0, 0.0

    stepper = GSL::Odeiv::Step.alloc(type, 2)

    y = GSL::Vector.alloc(1.0, 1.0)
    yerr = GSL::Vector.alloc(2)

    while t < 5.0
      stepper.apply(t, h, y, yerr, @rhs_func_exp)

      del = ((y[1] - Math.exp(t + h)) / y[1]).abs
      delmax = GSL.MAX_DBL(del, delmax)

      refute del > (n + 1.0) * 2.0 * b,
        "#{stepper.name}, exponential [0,5], max relative error = #{delmax}"

      n += 1
      t += h
    end
  end

  def _test_stepper_stiff(type, h, b)
    delmax, n, t = 0.0, 0, 0.0

    stepper = GSL::Odeiv::Step.alloc(type, 2)

    y = GSL::Vector.alloc(1.0, 0.0)
    yerr = GSL::Vector.alloc(2)

    while t < 5.0
      stepper.apply(t, h, y, yerr, @rhs_func_stiff)

      if t > 0.04
        arg = t + h

        e1 = Math.exp(-arg)
        e2 = Math.exp(-1000.0 * arg)

        del = ((y[0] - 2.0 * e1 + e2) / y[0]).abs
        delmax = GSL.MAX_DBL(del, delmax)

        refute del > (n + 1.0) * 100.0 * b,
          "#{stepper.name}, stiff [0,5], max relative error = #{delmax}"
      end

      n += 1
      t += h
    end
  end

  def _test_evolve_system_flat(type, sys, t0, t1, h, y, yfin, err_target, desc, control = nil)
    stepper = GSL::Odeiv::Step.alloc(type, sys.dimension)
    evolve  = GSL::Odeiv::Evolve.alloc(sys.dimension)

    while t0 < t1
      t0, h, _ = evolve.apply(control, stepper, sys, t0, t1, h, y)
    end

    frac = ((y[1] - yfin[1]) / yfin[1]).abs + ((y[0] - yfin[0]) / yfin[0]).abs
    refute frac > 2.0 * evolve.count * err_target, "#{stepper.name}, #{desc}, evolve" <<
      ", #{control ? 'standard' : 'no'} control, relative error = #{frac}"
  end

  def _test_evolve_system(type, sys, t0, t1, h, y, yfin, err_target, desc)
    _test_evolve_system_flat(type, sys, t0, t1, h, y, yfin, err_target, desc,
      GSL::Odeiv::Control.standard_alloc(0.0, err_target, 1.0, 1.0))
  end

  def _test_evolve_sin(type, h, err)
    y = GSL::Vector.alloc(1.0, 0.0)
    yfin = GSL::Vector.alloc(Math.cos(2.0), Math.sin(2.0))

    _test_evolve_system(type, @rhs_func_sin, 0.0, 2.0, h, y, yfin, err, 'sin [0,2]')
  end

  def _test_evolve_exp(type, h, err)
    eee = Math.exp(5.0)

    y = GSL::Vector.alloc(1.0, 1.0)
    yfin = GSL::Vector.alloc(eee, eee)

    _test_evolve_system(type, @rhs_func_exp, 0.0, 5.0, h, y, yfin, err, 'exp [0,5]')
  end

  def _test_evolve_stiff1(type, h, err, arg = 1.0)
    e1 = Math.exp(-arg)
    e2 = Math.exp(-1000.0 * arg)

    y = GSL::Vector.alloc(1.0, 0.0)
    yfin = GSL::Vector.alloc(2.0 * e1 - e2, -e1 + e2)

    _test_evolve_system(type, @rhs_func_stiff, 0.0, arg, h, y, yfin, err, 'stiff [0,%d]' % arg)
  end

  def _test_evolve_stiff5(type, h, err)
    _test_evolve_stiff1(type, h, err, 5.0)
  end

end
