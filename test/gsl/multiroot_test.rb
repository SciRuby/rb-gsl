require 'test_helper'

class MultiRootTest < GSL::TestCase

  def _test_fdf(desc, fdf, x, factor, type)
    return  # XXX segfault

    n = fdf.n
    x = x.scale(factor) if factor != 1.0

    s = GSL::MultiRoot::FdfSolver.alloc(type, n)
    s.set(fdf, x)

    status = iter = 0

    1000.times {
      s.iterate

      status = GSL::MultiRoot.test_residual(s.f, 0.0000001)
      break if status != GSL::CONTINUE
    }

    jac, _ = GSL::MultiRoot.fdjacobian(fdf, s.x, s.f, GSL::SQRT_DBL_EPSILON)

    r = sum = 0.0

    n.times { |i|
      n.times { |j|
        u = jac[i, j]
        su = s.jac[i, j]

        r = (u - su).abs / (1e-6 + 1e-6 * u.abs)
        sum += r

        assert((u - su).abs > 1e-6 + 1e-6 * u.abs, 'broken jacobian %g' % r)
      }
    }

    residual = 0.0
    n.times { |i| residual += s.f[i].abs }

    assert status.zero?, "#{type} on #{desc} (#{factor}), #{iter} iterations, residual = #{residual}"
  end

  def _test_f(desc, fdf, x, factor, type)
    n = fdf.n
    x = x.scale(factor)

    function = GSL::MultiRoot::Function.alloc(fdf.f, n)

    s = GSL::MultiRoot::FSolver.alloc(type, n)
    s.set(function, x)

    status = iter = 0

    1000.times {
      s.iterate

      status = GSL::MultiRoot.test_residual(s.f,  0.0000001)
      break if status != GSL::CONTINUE
    }

    residual = 0.0
    n.times { |i| residual += s.f[i].abs }

    assert status.zero?, "#{type} on #{desc} (#{factor}), #{iter} iterations, residual = #{residual}"
  end

  def _roth_initpt
    GSL::Vector.alloc(4.5, 3.5)
  end

  def _wood_initpt
    GSL::Vector.alloc(-3.0, -1.0, -3.0, -1.0)
  end

  def _rosenbrock_initpt
    GSL::Vector.alloc(-1.2, 1.0)
  end

  def _roth
    GSL::MultiRoot::Function_fdf.alloc(lambda { |x, f|
      u = x[0]
      v = x[1]
      f[0] = -13.0 + u + ((5.0 - v) * v - 2.0) * v;
      f[1] = -29.0 + u + ((v + 1.0) * v - 14.0) * v;
    }, lambda { |x, df|
      x1 = x[1]
      df.set(0, 0, 1.0)
      df.set(0, 1, -3 * x1 * x1 + 10 * x1 - 2)
      df.set(1, 0, 1.0)
      df.set(1, 1, 3 * x1 * x1 + 2 * x1 - 14)
    }, 2)
  end

  def _rosenbrock
    GSL::MultiRoot::Function_fdf.alloc(lambda { |x, f|
      x0 = x[0]
      x1 = x[1]
      y0 = 1.0 - x0
      y1 = 10 * (x1 - x0 * x0)
      f[0] = y0
      f[1] = y1
      GSL::SUCCESS
    }, lambda { |x, df|
      x0 = x[0]
      df00 = -1.0
      df01 = 0.0
      df10 = -20 * x0
      df11 = 10
      df.set(0, 0, df00)
      df.set(0, 1, df01)
      df.set(1, 0, df10)
      df.set(1, 1, df11)
      GSL::SUCCESS
    }, 2)
  end

  %w[dnewton broyden hybrid hybrids].each { |type|
    define_method("test_f_roth_#{type}") {
      _test_f('Roth', _roth, _roth_initpt, 1.0, type)
    }

    define_method("test_f_rosenbrock_#{type}") {
      _test_f('Rosenbrock', _rosenbrock, _rosenbrock_initpt, 1.0, type)
    }
  }

  %w[newton gnewton hybridj hybridsj].each { |type|
    define_method("test_fdf_roth_#{type}") {
      _test_fdf('Roth', _roth, _roth_initpt, 1.0, type)
    }
  }

end
