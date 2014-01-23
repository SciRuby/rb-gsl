require 'test_helper'

class MultiMinTest < GSL::TestCase

  def _test_fdf(desc, f, x, type)
    s = GSL::MultiMin::FdfMinimizer.alloc(type, f.n)
    s.set(f, x, 0.1 * GSL::Blas.dnrm2(x), 0.1)

    status = iter = 0

    begin
      iter += 1
      s.iterate

      status = GSL::MultiMin.test_gradient(s.gradient, 1e-3)
    end while iter < 5000 and status == GSL::CONTINUE

    status |= s.f.abs > 1e-5 ? 1 : 0
    assert status.zero?, "#{s.name}, on #{desc}: #{iter} iterations, f(x)=#{s.f}"
  end

  def _test_f(desc, f, x)
    step_size = GSL::Vector.alloc(f.n)
    f.n.times { |i| step_size[i] = 1 }

    s = GSL::MultiMin::FMinimizer.alloc('nmsimplex', f.n)
    s.set(f, x, step_size)

    status = iter = 0

    begin
      s.iterate
      status = GSL::MultiMin.test_size(s.size, 1e-3)
    end while iter < 5000 and status == GSL::CONTINUE

    status |= s.fval.abs > 1e-5 ? 1 : 0
    assert status.zero?, "#{s.name}, on #{desc}: #{iter} iterations, f(x)=#{s.fval}"

    s = GSL::MultiMin::FMinimizer.alloc('nmsimplex2rand', f.n)
    s.set(f, x, step_size)

    status = iter = 0

    begin
      s.iterate
      status = GSL::MultiMin.test_size(s.size, 1e-3)
    end while iter < 5000 and status == GSL::CONTINUE

    status |= s.fval.abs > 1e-5 ? 1 : 0
    assert status.zero?, "#{s.name}, on #{desc}: #{iter} iterations, f(x)=#{s.fval}"
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

  def _roth_f
    lambda { |x|
      u = x[0]
      v = x[1]
      a = -13.0 + u + ((5.0 - v) * v - 2.0) * v;
      b = -29.0 + u + ((v + 1.0) * v - 14.0) * v;
      a * a + b * b;
    }
  end

  def _wood_f
    lambda { |x|
      u1 = x[0]
      u2 = x[1]
      u3 = x[2]
      u4 = x[3]
      t1 = u1 * u1 - u2
      t2 = u3 * u3 - u4
      100 * t1 * t1 + (1 - u1) * (1 - u1) + 90 * t2 * t2 + (1 - u3) * (1 - u3) + 10.1 * ((1 - u2) * (1 - u2) + (1 - u4) * (1 - u4)) + 19.8 * (1 - u2) * (1 - u4)
    }
  end

  def _rosenbrock_f
    lambda { |x|
      u = x[0]
      v = x[1]
      a = u - 1
      b = u * u - v
      a * a + 10.0 * b * b
    }
  end

  def _rothdf
    GSL::MultiMin::Function_fdf.alloc(_roth_f, lambda { |x, df|
      u = x[0]
      v = x[1]
      a = -13.0 + u + ((5.0 - v) * v - 2.0) * v
      b = -29.0 + u + ((v + 1.0) * v - 14.0) * v
      c = -2 + v * (10 - 3 * v)
      d = -14 + v * (2 + 3 * v)
      df[0] = 2 * a + 2 * b
      df[1] = 2 * a * c + 2 * b * d
    }, 2)
  end

  def _wooddf
    GSL::MultiMin::Function_fdf.alloc(_wood_f, lambda { |x, df|
      u1 = x[0]
      u2 = x[1]
      u3 = x[2]
      u4 = x[3]
      t1 = u1 * u1 - u2
      t2 = u3 * u3 - u4
      df[0] = 400 * u1 * t1 - 2 * (1 - u1)
      df[1] = -200 * t1 - 20.2 * (1 - u2) - 19.8 * (1 - u4)
      df[2] = 360 * u3 * t2 - 2 * (1 - u3)
      df[3] = -180 * t2 - 20.2 * (1 - u4) - 19.8 * (1 - u2)
    }, 4)
  end

  def _rosenbrockdf
    GSL::MultiMin::Function_fdf.alloc(_rosenbrock_f, lambda { |x, df|
      u = x[0]
      v = x[1]
      a = u - 1
      b = u * u - v
      df[0] = 2 * a + 40 * u * b
      df[1] = -20 * b
    }, 2)
  end

  fdfminimizers = %w[steepest_descent conjugate_pr conjugate_fr vector_bfgs]
  fdfminimizers << 'vector_bfgs2' if GSL::GSL_VERSION >= '1.8.90'

  fdfminimizers.each { |type|
    define_method("test_fdf_roth_#{type}") { _test_fdf('Roth', _rothdf, _roth_initpt, type) }
    define_method("test_fdf_wood_#{type}") { _test_fdf('Wood', _wooddf, _wood_initpt, type) }
    define_method("test_fdf_rosenbrock_#{type}") { _test_fdf('Rosenbrock', _rosenbrockdf, _rosenbrock_initpt, type) }
  }

  def test_f_roth
    _test_f('Roth', GSL::MultiMin::Function.alloc(_roth_f, 2), _roth_initpt)
  end

  def test_f_wood
    _test_f('Wood', GSL::MultiMin::Function.alloc(_wood_f, 4), _wood_initpt)
  end

  def test_f_rosenbrock
    _test_f('Rosenbrock', GSL::MultiMin::Function.alloc(_rosenbrock_f, 2), _rosenbrock_initpt)
  end

end
