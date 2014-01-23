require 'test_helper'

class RootsTest < GSL::TestCase

  EPSREL = 10.0 * GSL::DBL_EPSILON
  EPSABS = 10.0 * GSL::DBL_EPSILON

  MAX_ITERATIONS = 150

  def _test_f(type, desc, f, lower, upper, correct)
    s = GSL::Root::FSolver.alloc(type)
    s.set(f, lower, upper)

    status = iter = 0

    begin
      iter += 1
      s.iterate

      r = s.root
      a = s.x_lower
      b = s.x_upper

      refute a > b, "interval is invalid (#{a},#{b})"
      refute r < a || r > b, "r lies outside interval #{r} (#{a},#{b})"

      status = GSL::Root.test_interval(a, b, EPSABS, EPSREL)
    end while status == GSL::CONTINUE && iter < MAX_ITERATIONS

    assert status.zero?, "#{s.name}, #{desc} (#{s.root} obs vs #{correct} expected)"
    assert iter <= MAX_ITERATIONS, 'exceeded maximum number of iterations'

    assert_tol r, correct, "incorrect precision (#{r} obs vs #{correct} expected)"
  end

  def _test_fdf(type, desc, fdf, root, correct)
    s = GSL::Root::FdfSolver.alloc(type)
    s.set(fdf, root)

    status = iter = 0

    begin
      iter += 1
      prev = s.root

      begin
        s.iterate
      rescue GSL::ERROR::EBADFUNC
        raise unless type == 'secant'
      end

      status = GSL::Root.test_delta(s.root, prev, EPSABS, EPSREL)
    end while status == GSL::CONTINUE && iter < MAX_ITERATIONS

    assert status.zero?, "#{s.name} #{desc} (#{s.root} obs vs #{correct} expected)"
    assert iter <= MAX_ITERATIONS, 'exceeded maximum number of iterations'

    assert_tol r = s.root, correct, "incorrect precision (#{r} obs vs #{correct} expected)"
  end

  def setup
    @func = GSL::Function.alloc { |x| GSL.pow(x, 20.0) - 1 }
    @fdf  = GSL::Function_fdf.alloc(@func.f, lambda { |x| 20.0 * GSL.pow(x, 19.0) })
  end

  %w[bisection brent falsepos].each { |type|
    define_method("test_f_#{type}") {
      _test_f(type, 'x^20 - 1 [0.1, 2]', @func, 0.1, 2.0, 1.0)
    }
  }

  %w[newton secant steffenson].each { |type|
    define_method("test_fdf_#{type}") {
      _test_fdf(type, 'x^{20} - 1 {0.9}', @fdf, 0.9, 1.0)
    }
  }

end
