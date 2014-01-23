require 'test_helper'

class DerivTest < GSL::TestCase

  def setup
    @f = [GSL::Function.alloc { |x| Math.exp(x) }]
    @df = [GSL::Function.alloc { |x| Math.exp(x) }]

    @f << GSL::Function.alloc { |x|
      x >= 0.0 ? x * Math.sqrt(x) : 0.0
    }
    @df << GSL::Function.alloc { |x|
      x >= 0.0 ? 1.5 * Math.sqrt(x) : 0.0
    }

    @f << GSL::Function.alloc { |x|
      x != 0.0 ? Math.sin(1.0 / x) : 0.0
    }
    @df << GSL::Function.alloc { |x|
      x != 0.0 ? -Math.cos(1.0 / x) / (x * x) : 0.0
    }

    @f << GSL::Function.alloc { |x| Math.exp(-x * x) }
    @df << GSL::Function.alloc { |x| -2.0 * x * Math.exp(-x * x) }

    @f << GSL::Function.alloc { |x| x * x }
    @df << GSL::Function.alloc { |x| 2.0 * x }

    @f << GSL::Function.alloc { |x| 1.0 / x }
    @df << GSL::Function.alloc { |x| -1.0 / (x * x) }
  end

  {
    'exp(x)'    => 1.0,
    'x^(3/2)'   => 0.1,
    'sin(1/x)'  => 0.45,
    'exp(-x^2)' => 0.5,
    'x^2'       => 0.0,
    '1/x'       => 10.0
  }.each_with_index { |(f, x), i|
    %w[central forward backward].each { |deriv|
      define_method("test_#{deriv}_#{i}") {
        _test_deriv(deriv, @f[i], @df[i], x, "#{f}, x=#{x}, #{deriv} deriv")
      }
    }
  }

  def _test_deriv(deriv, f, df, x, desc)
    expected, h = df.eval(x), 1e-4

    result, abserr = f.send("deriv_#{deriv}", x, h)
    assert_abs result, expected, GSL.MIN(h, expected.abs) + GSL::DBL_EPSILON, desc

    if abserr < (diff = (result - expected).abs)
      assert_factor abserr, (result - expected).abs, 2, desc + ' error estimate'
    else
      zero = result == expected || expected.zero?
      assert_abs abserr, zero ? 0.0 : diff, zero ? 1e-6 : 1e6 * diff, desc + ' abserr'
    end
  end

end
